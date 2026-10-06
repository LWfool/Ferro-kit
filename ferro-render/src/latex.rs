//! Inline LaTeX → Unicode text, for formulas shown on a terminal.

/// Inline LaTeX as Unicode text, by four rules rather than a lookup table:
/// `_x` / `^x`, `\<greek>`, `\text{}`-like wrappers, and a few symbols.
///
/// A script is converted only when every character of it has a Unicode
/// sub/superscript; otherwise it stays as `_x` / `_{xy}` — `τ_c`, not `τc`,
/// which would be a different reading. A command not known here stays as
/// written, so the worst case is TeX on screen, never a wrong symbol.
pub(crate) fn latex(src: &str) -> String {
    let chars: Vec<char> = src.chars().collect();
    let mut i = 0;
    let mut out = String::new();
    while i < chars.len() {
        match chars[i] {
            c @ ('_' | '^') => {
                i += 1;
                if i < chars.len() {
                    let word = names_a_word(&chars[i..]);
                    let arg = atom(&chars, &mut i);
                    out.push_str(&script(c, &arg, word));
                } else {
                    out.push(c);
                }
            }
            _ => out.push_str(&atom(&chars, &mut i)),
        }
    }
    out
}

/// Whether a script argument is `\text{..}`-like: a word label, not variables.
fn names_a_word(rest: &[char]) -> bool {
    let name: String = rest.iter().skip(1).take_while(|c| c.is_ascii_alphabetic()).collect();
    rest.first() == Some(&'\\') && matches!(name.as_str(), "text" | "mathrm" | "operatorname")
}

/// One atom at `i`: a command with its arguments, a braced group, or a char.
fn atom(chars: &[char], i: &mut usize) -> String {
    match chars[*i] {
        '{' => latex(&braced(chars, i)),
        '\\' => command(chars, i),
        c => {
            *i += 1;
            c.to_string()
        }
    }
}

/// The raw text inside the `{...}` at `i`, leaving `i` past its closing brace.
fn braced(chars: &[char], i: &mut usize) -> String {
    let mut depth = 0;
    let start = *i + 1;
    while *i < chars.len() {
        match chars[*i] {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    *i += 1;
                    return chars[start..*i - 1].iter().collect();
                }
            }
            _ => {}
        }
        *i += 1;
    }
    chars[start.min(chars.len())..].iter().collect()
}

fn command(chars: &[char], i: &mut usize) -> String {
    *i += 1;
    let start = *i;
    while *i < chars.len() && chars[*i].is_ascii_alphabetic() {
        *i += 1;
    }
    let name: String = chars[start..*i].iter().collect();
    if name.is_empty() {
        // \, \  \; 是间距,\! 是负间距,\{ \} \_ 是转义
        let c = chars.get(*i).copied();
        *i += 1;
        return match c {
            // 后面本来就有空格时不再添一个,否则 f_\min\, m 印成两个空格
            Some(',' | ' ' | ';' | ':') if chars.get(*i) == Some(&' ') => String::new(),
            Some(',' | ' ' | ';' | ':') => " ".to_string(),
            Some('!') => String::new(),
            Some(c) => c.to_string(),
            None => "\\".to_string(),
        };
    }
    let arg = |i: &mut usize| {
        while *i < chars.len() && chars[*i] == ' ' {
            *i += 1;
        }
        if *i < chars.len() { atom(chars, i) } else { String::new() }
    };
    if let Some(g) = greek(&name) {
        // TeX 里命令名后的空格只是分隔符:\Delta r 印出来是 Δr。
        // 后面是运算符时那个空格留着,否则 \tau \leq N 会变成 τ≤ N
        if chars.get(*i) == Some(&' ') && chars.get(*i + 1).is_some_and(|c| c.is_alphanumeric()) {
            *i += 1;
        }
        return g.to_string();
    }
    let sym = match name.as_str() {
        "cdot" => "·",
        "times" => "×",
        "to" | "rightarrow" => "→",
        "infty" => "∞",
        "ell" => "ℓ",
        "in" => "∈",
        "int" => "∫",
        "propto" => "∝",
        "sim" => "∼",
        "parallel" => "∥",
        "wedge" => "∧",
        "leftarrow" => "←",
        "bmod" => "mod",
        "quad" => " ",
        "qquad" => "  ",
        "neq" | "ne" => "≠",
        "leq" | "le" => "≤",
        "geq" | "ge" => "≥",
        "approx" => "≈",
        "pm" => "±",
        "sum" => "Σ",
        "langle" => "⟨",
        "rangle" => "⟩",
        "ldots" | "dots" | "cdots" => "…",
        "lfloor" => "⌊",
        "rfloor" => "⌋",
        "lceil" => "⌈",
        "rceil" => "⌉",
        "left" | "right" | "bigl" | "bigr" | "big" | "Big" => "",
        "sin" | "cos" | "tan" | "exp" | "ln" | "log" | "det" | "min" | "max" | "arccos" | "diag" | "lim" => {
            return name;
        }
        "text" | "mathrm" | "mathbf" | "mathit" | "mathsf" | "boldsymbol" | "operatorname" | "mathrel"
        | "underbrace" => return arg(i),
        // 重音用 Unicode 组合字符叠在参数后面;多字符参数只叠在最后一个字符上
        "hat" => return format!("{}\u{302}", arg(i)),
        "bar" | "overline" => return format!("{}\u{304}", arg(i)),
        "vec" | "overrightarrow" => return format!("{}\u{20d7}", arg(i)),
        "frac" | "tfrac" | "dfrac" => {
            let (a, b) = (arg(i), arg(i));
            return format!("{}/{}", grouped(&a), grouped(&b));
        }
        "sqrt" => return format!("√{}", grouped(&arg(i))),
        _ => {
            // 认不出的命令连同参数原样保留
            let mut raw = format!("\\{name}");
            if chars.get(*i) == Some(&'{') {
                raw.push_str(&format!("{{{}}}", braced(chars, i)));
            }
            return raw;
        }
    };
    sym.to_string()
}

/// `x` as is, `ab` as `(ab)` — what a `/` needs to read unambiguously.
fn grouped(s: &str) -> String {
    if s.chars().count() <= 1 || s.chars().all(|c| c.is_ascii_digit()) {
        s.to_string()
    } else {
        format!("({s})")
    }
}

fn greek(name: &str) -> Option<char> {
    const TABLE: [(&str, char); 36] = [
        ("alpha", 'α'), ("beta", 'β'), ("gamma", 'γ'), ("delta", 'δ'), ("epsilon", 'ε'),
        ("varepsilon", 'ε'), ("zeta", 'ζ'), ("eta", 'η'), ("theta", 'θ'), ("kappa", 'κ'),
        ("lambda", 'λ'), ("mu", 'μ'), ("nu", 'ν'), ("xi", 'ξ'), ("pi", 'π'),
        ("rho", 'ρ'), ("sigma", 'σ'), ("tau", 'τ'), ("phi", 'φ'), ("varphi", 'φ'),
        ("chi", 'χ'), ("psi", 'ψ'), ("omega", 'ω'), ("iota", 'ι'),
        ("Gamma", 'Γ'), ("Delta", 'Δ'), ("Theta", 'Θ'), ("Lambda", 'Λ'), ("Xi", 'Ξ'),
        ("Pi", 'Π'), ("Sigma", 'Σ'), ("Phi", 'Φ'), ("Psi", 'Ψ'), ("Omega", 'Ω'),
        ("Upsilon", 'Υ'), ("upsilon", 'υ'),
    ];
    TABLE.iter().find(|(n, _)| *n == name).map(|(_, c)| *c)
}

/// `_` or `^` applied to already-converted text.
///
/// A `word` (from `\text{min}`) is never converted: `N_{frames}` and
/// `Nₐₜₒₘₛ` side by side — one word happens to have every letter, the other
/// not — reads worse than the TeX-ish form for both.
fn script(kind: char, arg: &str, word: bool) -> String {
    // Unicode 的下标字母不全(没有 b c d f g q w y z),上标字母只有 n 与 i
    const SUB: &str = "0₀1₁2₂3₃4₄5₅6₆7₇8₈9₉+₊-₋=₌(₍)₎aₐeₑoₒxₓhₕkₖlₗmₘnₙpₚsₛtₜiᵢjⱼrᵣuᵤvᵥ";
    const SUP: &str = "0⁰1¹2²3³4⁴5⁵6⁶7⁷8⁸9⁹+⁺-⁻=⁼(⁽)⁾nⁿiⁱ";
    let table: Vec<char> = (if kind == '_' { SUB } else { SUP }).chars().collect();
    let map = |c: char| table.chunks(2).find(|p| p[0] == c).map(|p| p[1]);
    if let Some(mapped) = arg.chars().map(map).collect::<Option<String>>() {
        if !mapped.is_empty() && !word {
            return mapped;
        }
    }
    if arg.chars().count() == 1 {
        format!("{kind}{arg}")
    } else {
        format!("{kind}{{{arg}}}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn latex_on_formulas_taken_from_the_manual() {
        let cases = [
            (r"g_{\alpha\beta}(r)", "g_{αβ}(r)"),
            (r"r_i = r_\text{min} + (i + 0.5) \Delta r", "rᵢ = r_{min} + (i + 0.5) Δr"),
            (r"\mathbf{f} = \mathbf{r} \cdot \mathbf{M}^{-1}", "f = r · M⁻¹"),
            (r"C_v(0) = \langle v^2 \rangle = 3 k_B T / m", "Cᵥ(0) = ⟨ v² ⟩ = 3 k_B T / m"),
            (r"\tau_c", "τ_c"),
            (r"\sum(\text{bridges}) \neq 2 \times |\text{O\_b}|", "Σ(bridges) ≠ 2 × |O_b|"),
            (r"r_i = r_\text{min} + (i + \tfrac{1}{2})\Delta r", "rᵢ = r_{min} + (i + 1/2)Δr"),
            (r"[k\Delta\theta,\ (k{+}1)\Delta\theta)", "[kΔθ, (k+1)Δθ)"),
            (r"\frac{4\pi\rho}{q}\sum_r r[g(r)-1]\sin(qr)\,\Delta r", "(4πρ)/qΣᵣ r[g(r)-1]sin(qr) Δr"),
            (r"r_\text{cut,AB}", "r_{cut,AB}"),
            (r"\mathbf{H} = \mathbf{U} \boldsymbol{\Sigma} \mathbf{V}^T", "H = U Σ V^T"),
            (r"Q^n_m", "Qⁿₘ"),
            (r"w_{ij}", "wᵢⱼ"),
            (r"N_\text{frames} \cdot N_\text{atoms}", "N_{frames} · N_{atoms}"),
            (r"p + \tau \leq N_\text{frames}", "p + τ ≤ N_{frames}"),
            (r"m_\mathrm{Al}=1", "m_{Al}=1"),
            (r"C_\ell(t)", "C_ℓ(t)"),
            (r"n \in \text{neighbours}", "n ∈ neighbours"),
            (r"S_{AA} \propto c_A^2", "S_{AA} ∝ c_A²"),
            (r"D_{\parallel}", "D_∥"),
            (r"\int P\,dr", "∫ P dr"),
            (r"\theta = \arccos\!\left(x\right)", "θ = arccos(x)"),
            (r"m \in [f_\min\, m_\max]", "m ∈ [fₘᵢₙ mₘₐₓ]"),
            (r"i_k = n \bmod n_k", "iₖ = n mod nₖ"),
            (r"\rho \mathrel{+}= x", "ρ += x"),
            (r"\hat{\mathbf{u}}", "u\u{302}"),
            (r"\bar\theta", "θ\u{304}"),
            (r"|\overrightarrow{BA}|", "|BA\u{20d7}|"),
        ];
        for (tex, want) in cases {
            assert_eq!(latex(tex), want, "输入: {tex}");
        }
    }

    #[test]
    fn unknown_commands_stay_as_written() {
        assert_eq!(latex(r"|\widetilde{BA}| < r"), r"|\widetilde{BA}| < r");
    }
}
