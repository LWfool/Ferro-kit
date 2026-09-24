//! Terminal rendering of the manual's markdown.
//!
//! Used only when stdout is a terminal: rendered output is escape codes and
//! box-drawing characters, which in a redirected file are noise, so
//! `ferro doc net > net.md` still gets the source as written.
//!
//! Three extension points that do not touch each other:
//!
//! - a new block type is a variant in [`split_blocks`] plus one `render_*`
//! - a new inline mark is one more flag on [`Look`], set in [`inline`]
//! - a change of style stays inside the `render_*` functions
//!
//! Syntax the splitter does not recognise comes out as [`Block::Raw`], so new
//! manual content can fail to be rendered but never be rendered wrong.

use std::fmt::Write as _;

/// What the terminal can be trusted to show.
#[derive(Clone, Copy, Debug)]
pub struct Style {
    /// SGR escape codes: bold, colour, underline.
    pub ansi: bool,
    /// Bullets and box-drawing characters outside ASCII.
    pub unicode: bool,
}

impl Style {
    pub fn for_terminal() -> Self {
        // Windows 的老 conhost 不认 ANSI,代码页也未必显示得了框线字符 ——
        // 在那里宁可朴素也要能读,不去赌控制台的配置
        if cfg!(windows) {
            return Style { ansi: false, unicode: false };
        }
        Style {
            ansi: std::env::var_os("NO_COLOR").is_none(),
            unicode: true,
        }
    }
}

/// Prose is not reflowed wider than this, however wide the terminal is.
const MAX_TEXT_WIDTH: usize = 100;

/// Columns of the terminal on stdout: the tty itself, then `$COLUMNS`, then 80.
pub fn terminal_width() -> usize {
    tty_width()
        .or_else(|| std::env::var("COLUMNS").ok()?.trim().parse().ok())
        .filter(|&w| w > 0)
        .unwrap_or(80)
}

#[cfg(unix)]
fn tty_width() -> Option<usize> {
    // SAFETY: winsize 是纯数据结构,全零是合法值;ioctl 只往里写,失败时返回 -1
    let mut ws: libc::winsize = unsafe { std::mem::zeroed() };
    let rc = unsafe { libc::ioctl(libc::STDOUT_FILENO, libc::TIOCGWINSZ, &mut ws) };
    (rc == 0 && ws.ws_col > 0).then_some(ws.ws_col as usize)
}

// Windows 不查控制台:退回 COLUMNS / 80,排版差一点但照样能读
#[cfg(not(unix))]
fn tty_width() -> Option<usize> {
    None
}

/// Renders a manual page for a terminal `width` columns wide.
pub fn render(md: &str, width: usize, style: Style) -> String {
    let text_width = width.clamp(20, MAX_TEXT_WIDTH);
    let blocks: Vec<String> = split_blocks(md)
        .iter()
        .map(|b| match b {
            Block::Heading { level, text } => render_heading(*level, text, style),
            Block::Para { indent, text } => render_para(*indent, text, text_width, style),
            Block::Code { indent, lines } => render_code(*indent, lines),
            Block::Quote(lines) => render_quote(lines, text_width, style),
            Block::List(items) => render_list(items, text_width, style),
            Block::Rule => rule(if style.unicode { '─' } else { '-' }, text_width),
            Block::Table(lines) => render_table(lines, width.max(20), style),
            Block::Raw(lines) => lines.join("\n"),
        })
        .collect();
    let mut out = blocks.join("\n\n");
    out.push('\n');
    out
}

// ─── 块切分 ─────────────────────────────────────────────────────────────

#[derive(Debug, PartialEq)]
enum Block<'a> {
    Heading { level: usize, text: &'a str },
    /// Source lines joined with single spaces; wrapping is redone at render.
    Para { indent: usize, text: String },
    /// Fence lines dropped, the fence's own indentation removed.
    Code { indent: usize, lines: Vec<&'a str> },
    Table(Vec<&'a str>),
    /// Lines with the `>` marker removed; an empty line separates paragraphs.
    Quote(Vec<&'a str>),
    List(Vec<Item>),
    Rule,
    Raw(Vec<&'a str>),
}

#[derive(Debug, PartialEq)]
struct Item {
    indent: usize,
    /// `-`, `*`, `+`, or a number with its dot.
    marker: String,
    text: String,
}

/// Cuts a page into blocks by looking only at how each line starts.
fn split_blocks(md: &str) -> Vec<Block<'_>> {
    let lines: Vec<&str> = md.lines().collect();
    let mut blocks = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let line = lines[i];
        let trimmed = line.trim_start();
        let indent = line.len() - trimmed.len();
        if trimmed.is_empty() {
            i += 1;
        } else if trimmed.starts_with("```") {
            // 列表项里的围栏整体缩进,正文去掉同样的缩进才对得齐
            let mut body = Vec::new();
            i += 1;
            while i < lines.len() && !lines[i].trim_start().starts_with("```") {
                let l = lines[i];
                let strip = l.len() - l.trim_start().len();
                body.push(&l[strip.min(indent)..]);
                i += 1;
            }
            i += 1; // 闭合围栏;文件末尾没闭合时越过末尾,循环条件兜住
            blocks.push(Block::Code { indent, lines: body });
        } else if trimmed.starts_with("$$") {
            // 块级公式不渲染(\begin{pmatrix} 在终端里无解),原样输出
            let start = i;
            let single = trimmed.len() > 2 && trimmed.trim_end().ends_with("$$");
            i += 1;
            if !single {
                while i < lines.len() {
                    let closed = lines[i].trim_end().ends_with("$$");
                    i += 1;
                    if closed {
                        break;
                    }
                }
            }
            blocks.push(Block::Raw(lines[start..i].to_vec()));
        } else if let Some(level) = heading_level(line) {
            blocks.push(Block::Heading { level, text: line[level..].trim() });
            i += 1;
        } else if is_rule(line) {
            blocks.push(Block::Rule);
            i += 1;
        } else if trimmed.starts_with('|') {
            let start = i;
            while i < lines.len() && lines[i].trim_start().starts_with('|') {
                i += 1;
            }
            blocks.push(Block::Table(lines[start..i].to_vec()));
        } else if trimmed.starts_with('>') {
            let mut body = Vec::new();
            while i < lines.len() && lines[i].trim_start().starts_with('>') {
                let rest = &lines[i].trim_start()[1..];
                body.push(rest.strip_prefix(' ').unwrap_or(rest));
                i += 1;
            }
            blocks.push(Block::Quote(body));
        } else if list_marker(trimmed).is_some() {
            let mut items: Vec<Item> = Vec::new();
            while i < lines.len() {
                let l = lines[i];
                let t = l.trim_start();
                if let Some((marker, rest)) = list_marker(t) {
                    items.push(Item {
                        indent: l.len() - t.len(),
                        marker: marker.to_string(),
                        text: rest.trim().to_string(),
                    });
                } else if t.is_empty() || starts_block(l) {
                    break;
                } else if let Some(last) = items.last_mut() {
                    // 续行:缩进的与不缩进的(惰性续行)都归上一项
                    last.text.push(' ');
                    last.text.push_str(t.trim_end());
                }
                i += 1;
            }
            blocks.push(Block::List(items));
        } else {
            let mut text = String::new();
            while i < lines.len() && !lines[i].trim().is_empty() && !starts_block(lines[i]) {
                if !text.is_empty() {
                    text.push(' ');
                }
                text.push_str(lines[i].trim());
                i += 1;
            }
            blocks.push(Block::Para { indent, text });
        }
    }
    blocks
}

/// Whether a line opens a block other than a paragraph — which also ends one.
fn starts_block(line: &str) -> bool {
    let t = line.trim_start();
    t.starts_with("```")
        || t.starts_with("$$")
        || t.starts_with('|')
        || t.starts_with('>')
        || heading_level(line).is_some()
        || is_rule(line)
        || list_marker(t).is_some()
}

fn heading_level(line: &str) -> Option<usize> {
    let level = line.chars().take_while(|&c| c == '#').count();
    ((1..=6).contains(&level) && line[level..].starts_with(' ')).then_some(level)
}

fn is_rule(line: &str) -> bool {
    line.len() >= 3 && line.chars().all(|c| c == '-')
}

/// `- x`, `* x`, `+ x` or `12. x` → the marker and the rest.
fn list_marker(t: &str) -> Option<(&str, &str)> {
    for bullet in ["- ", "* ", "+ "] {
        if let Some(rest) = t.strip_prefix(bullet) {
            return Some((&t[..1], rest));
        }
    }
    let digits = t.chars().take_while(|c| c.is_ascii_digit()).count();
    if digits > 0 && t[digits..].starts_with(". ") {
        return Some((&t[..digits + 1], &t[digits + 2..]));
    }
    None
}

// ─── 行内 ───────────────────────────────────────────────────────────────

/// The marks in effect on a run of text. Marks nest, so this is a set of
/// flags rather than one kind per span.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Look {
    bold: bool,
    italic: bool,
    code: bool,
    math: bool,
    link: bool,
    dim: bool,
}

impl Look {
    fn or(self, o: Look) -> Look {
        Look {
            bold: self.bold || o.bold,
            italic: self.italic || o.italic,
            code: self.code || o.code,
            math: self.math || o.math,
            link: self.link || o.link,
            dim: self.dim || o.dim,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Span {
    text: String,
    look: Look,
}

/// The one inline parser: a single left-to-right scan.
///
/// Not a chain of `replace` calls — marks exclude each other (`**` inside code
/// is not bold, a link's text can be bold), and a replace chain breaks the
/// moment a third mark is added. Code and math are opaque: nothing inside them
/// is parsed.
fn inline(text: &str) -> Vec<Span> {
    let chars: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut buf = String::new();
    let mut look = Look::default();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        match c {
            '`' => {
                let run = run_len(&chars, i, '`');
                match find_run(&chars, i + run, '`', run) {
                    Some(j) => {
                        flush(&mut out, &mut buf, look);
                        let inner: String = chars[i + run..j].iter().collect();
                        out.push(Span {
                            text: inner.trim().to_string(),
                            look: look.or(Look { code: true, ..Look::default() }),
                        });
                        i = j + run;
                    }
                    None => {
                        buf.extend(&chars[i..i + run]);
                        i += run;
                    }
                }
            }
            '$' if chars.get(i + 1).is_some_and(|n| !n.is_whitespace() && *n != '$') => {
                match chars[i + 1..].iter().position(|&x| x == '$') {
                    Some(k) => {
                        flush(&mut out, &mut buf, look);
                        out.push(Span {
                            text: chars[i + 1..i + 1 + k].iter().collect(),
                            look: look.or(Look { math: true, ..Look::default() }),
                        });
                        i += k + 2;
                    }
                    None => {
                        buf.push(c);
                        i += 1;
                    }
                }
            }
            '*' if chars.get(i + 1) == Some(&'*') => {
                // 只有后面还有闭合的 ** 才当粗体开头,否则是字面的星号
                if look.bold || contains(&chars[i + 2..], &['*', '*']) {
                    flush(&mut out, &mut buf, look);
                    look.bold = !look.bold;
                } else {
                    buf.push_str("**");
                }
                i += 2;
            }
            '*' => {
                let opens = !look.italic
                    && chars.get(i + 1).is_some_and(|n| !n.is_whitespace())
                    && chars[i + 1..].contains(&'*');
                let closes = look.italic && i > 0 && !chars[i - 1].is_whitespace();
                if opens || closes {
                    flush(&mut out, &mut buf, look);
                    look.italic = !look.italic;
                } else {
                    buf.push(c);
                }
                i += 1;
            }
            '[' => match link_at(&chars, i) {
                Some((label, url, end)) => {
                    flush(&mut out, &mut buf, look);
                    let base = look.or(Look { link: true, ..Look::default() });
                    for s in inline(&label) {
                        out.push(Span { text: s.text, look: base.or(s.look) });
                    }
                    out.push(Span { text: format!(" ({url})"), look: Look { dim: true, ..Look::default() } });
                    i = end;
                }
                None => {
                    buf.push(c);
                    i += 1;
                }
            },
            '\\' if chars.get(i + 1).is_some_and(|n| n.is_ascii_punctuation()) => {
                buf.push(chars[i + 1]);
                i += 2;
            }
            _ => {
                buf.push(c);
                i += 1;
            }
        }
    }
    flush(&mut out, &mut buf, look);
    out
}

fn flush(out: &mut Vec<Span>, buf: &mut String, look: Look) {
    if !buf.is_empty() {
        out.push(Span { text: std::mem::take(buf), look });
    }
}

fn run_len(chars: &[char], from: usize, c: char) -> usize {
    chars[from..].iter().take_while(|&&x| x == c).count()
}

/// Start of the next run of exactly `len` copies of `c` at or after `from`.
fn find_run(chars: &[char], from: usize, c: char, len: usize) -> Option<usize> {
    let mut j = from;
    while j < chars.len() {
        if chars[j] == c {
            let n = run_len(chars, j, c);
            if n == len {
                return Some(j);
            }
            j += n;
        } else {
            j += 1;
        }
    }
    None
}

fn contains(hay: &[char], needle: &[char]) -> bool {
    hay.windows(needle.len()).any(|w| w == needle)
}

/// `[label](url)` starting at `i` → label, url, and the index just past it.
fn link_at(chars: &[char], i: usize) -> Option<(String, String, usize)> {
    let mut depth = 0;
    let mut close = None;
    for (k, &c) in chars.iter().enumerate().skip(i) {
        match c {
            '[' => depth += 1,
            ']' => {
                depth -= 1;
                if depth == 0 {
                    close = Some(k);
                    break;
                }
            }
            _ => {}
        }
    }
    let close = close?;
    if chars.get(close + 1) != Some(&'(') {
        return None;
    }
    let end = close + 2 + chars[close + 2..].iter().position(|&c| c == ')')?;
    Some((
        chars[i + 1..close].iter().collect(),
        chars[close + 2..end].iter().collect(),
        end + 1,
    ))
}

// ─── 着色与折行 ─────────────────────────────────────────────────────────

/// The visible text of a span, with the delimiters a plain terminal needs to
/// tell code apart from prose.
fn shown(s: &Span, style: Style) -> String {
    if s.look.math {
        // 希腊字母与上下标在 ASCII 终端上会变成问号,那里宁可保留 TeX 原文
        if style.unicode { latex(&s.text) } else { format!("${}$", s.text) }
    } else if s.look.code && !style.ansi {
        format!("`{}`", s.text)
    } else {
        s.text.clone()
    }
}

fn paint(s: &Span, style: Style) -> String {
    let text = shown(s, style);
    if !style.ansi {
        return text;
    }
    let l = s.look;
    let codes: Vec<&str> = [
        (l.bold, "1"),
        (l.dim, "2"),
        (l.italic, "3"),
        (l.link, "4"),
        (l.code, "36"),
    ]
    .iter()
    .filter(|(on, _)| *on)
    .map(|(_, c)| *c)
    .collect();
    if codes.is_empty() {
        text
    } else {
        format!("\x1b[{}m{text}\x1b[0m", codes.join(";"))
    }
}

fn visible_len(spans: &[Span], style: Style) -> usize {
    spans.iter().map(|s| shown(s, style).chars().count()).sum()
}

/// Greedy word wrap of inline spans to `width` visible columns.
///
/// A word is everything between two whitespaces, across span boundaries —
/// ``` `x`, ``` keeps its comma. Code and math stay whole when they fit on a
/// line; one longer than the line breaks at its own spaces, and a single word
/// longer than the line overflows rather than being cut.
fn wrap(spans: &[Span], width: usize, style: Style) -> Vec<String> {
    let mut words: Vec<Vec<Span>> = Vec::new();
    // 上一个词是否还能继续往后粘(中间没遇到空白)
    let mut open = false;
    for s in spans {
        let opaque = s.look.code || s.look.math;
        if opaque && shown(s, style).chars().count() <= width {
            if !open {
                words.push(Vec::new());
            }
            words.last_mut().unwrap().push(s.clone());
            open = true;
            continue;
        }
        // 比整行还长的代码要在内部空格处断开。朴素样式下定界符先并进文本,
        // 否则每个碎片各带一对反引号,读起来像一串独立的代码
        let owned;
        let s = if opaque && (s.look.math || !style.ansi) {
            owned = Span { text: shown(s, style), look: Look { code: false, math: false, ..s.look } };
            &owned
        } else {
            s
        };
        let mut buf = String::new();
        for c in s.text.chars() {
            if c.is_whitespace() {
                if !buf.is_empty() {
                    if !open {
                        words.push(Vec::new());
                    }
                    words.last_mut().unwrap().push(Span { text: std::mem::take(&mut buf), look: s.look });
                }
                open = false;
            } else {
                buf.push(c);
            }
        }
        if !buf.is_empty() {
            if !open {
                words.push(Vec::new());
            }
            words.last_mut().unwrap().push(Span { text: buf, look: s.look });
            open = true;
        }
    }

    let mut lines = Vec::new();
    let mut line = String::new();
    let mut used = 0;
    for w in &words {
        let n = visible_len(w, style);
        if used > 0 && used + 1 + n > width {
            lines.push(std::mem::take(&mut line));
            used = 0;
        }
        if used > 0 {
            line.push(' ');
            used += 1;
        }
        for s in w {
            line.push_str(&paint(s, style));
        }
        used += n;
    }
    if used > 0 || lines.is_empty() {
        lines.push(line);
    }
    lines
}

// ─── 行内公式 ───────────────────────────────────────────────────────────

/// Inline LaTeX as Unicode text, by four rules rather than a lookup table:
/// `_x` / `^x`, `\<greek>`, `\text{}`-like wrappers, and a few symbols.
///
/// A script is converted only when every character of it has a Unicode
/// sub/superscript; otherwise it stays as `_x` / `_{xy}` — `τ_c`, not `τc`,
/// which would be a different reading. A command not known here stays as
/// written, so the worst case is TeX on screen, never a wrong symbol.
fn latex(src: &str) -> String {
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
        // \, \  \; 是间距,\{ \} \_ 是转义
        let c = chars.get(*i).copied();
        *i += 1;
        return match c {
            Some(',' | ' ' | ';' | ':') => " ".to_string(),
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
        "sin" | "cos" | "tan" | "exp" | "ln" | "log" | "det" | "min" | "max" | "arccos" | "diag" => {
            return name;
        }
        "text" | "mathrm" | "mathbf" | "mathit" | "boldsymbol" | "operatorname" => return arg(i),
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

// ─── 各类块 ─────────────────────────────────────────────────────────────

fn rule(c: char, width: usize) -> String {
    std::iter::repeat_n(c, width).collect()
}

fn render_heading(level: usize, text: &str, style: Style) -> String {
    let spans = inline(text);
    let bold: Vec<Span> = spans
        .iter()
        .map(|s| Span { text: s.text.clone(), look: s.look.or(Look { bold: true, ..Look::default() }) })
        .collect();
    let mut line: String = bold.iter().map(|s| paint(s, style)).collect();
    let n = visible_len(&bold, style);
    match level {
        1 => write!(line, "\n{}", rule(if style.unicode { '═' } else { '=' }, n)).unwrap(),
        2 => write!(line, "\n{}", rule(if style.unicode { '─' } else { '-' }, n)).unwrap(),
        // 没有粗体可用时,三级以下标题靠原有的 # 与正文区分
        _ if !style.ansi => line = format!("{} {line}", "#".repeat(level)),
        _ => {}
    }
    line
}

fn render_para(indent: usize, text: &str, width: usize, style: Style) -> String {
    let pad = " ".repeat(indent);
    wrap(&inline(text), width.saturating_sub(indent).max(20), style)
        .iter()
        .map(|l| format!("{pad}{l}"))
        .collect::<Vec<_>>()
        .join("\n")
}

fn render_code(indent: usize, lines: &[&str]) -> String {
    let pad = " ".repeat(indent + 4);
    lines
        .iter()
        .map(|l| if l.is_empty() { String::new() } else { format!("{pad}{l}") })
        .collect::<Vec<_>>()
        .join("\n")
}

fn render_quote(lines: &[&str], width: usize, style: Style) -> String {
    let bar = if style.unicode { "│" } else { "|" };
    let mut out: Vec<String> = Vec::new();
    for (k, para) in lines.split(|l| l.trim().is_empty()).enumerate() {
        if k > 0 {
            out.push(bar.to_string());
        }
        let text = para.iter().map(|l| l.trim()).collect::<Vec<_>>().join(" ");
        for l in wrap(&inline(&text), width - 2, style) {
            out.push(format!("{bar} {l}"));
        }
    }
    out.join("\n")
}

fn render_list(items: &[Item], width: usize, style: Style) -> String {
    let mut out: Vec<String> = Vec::new();
    for item in items {
        let bullet = if item.marker.ends_with('.') || !style.unicode {
            item.marker.as_str()
        } else {
            "•"
        };
        let head = format!("{}{bullet} ", " ".repeat(item.indent));
        let hang = " ".repeat(head.chars().count());
        let body = wrap(&inline(&item.text), width.saturating_sub(hang.len()).max(20), style);
        for (k, l) in body.iter().enumerate() {
            out.push(format!("{}{l}", if k == 0 { &head } else { &hang }));
        }
    }
    out.join("\n")
}

/// How a column's cells sit inside it, from the `|:-:|` row.
#[derive(Clone, Copy, Debug, PartialEq)]
enum Align {
    Left,
    Center,
    Right,
}

/// Splits one table row at the `|` that are not escaped as `\|`.
fn cells(row: &str) -> Vec<&str> {
    let row = row.trim();
    let row = row.strip_prefix('|').unwrap_or(row);
    let row = row.strip_suffix('|').filter(|r| !r.ends_with('\\')).unwrap_or(row);
    let mut out = Vec::new();
    let mut start = 0;
    let bytes = row.as_bytes();
    for (k, &b) in bytes.iter().enumerate() {
        if b == b'|' && (k == 0 || bytes[k - 1] != b'\\') {
            out.push(row[start..k].trim());
            start = k + 1;
        }
    }
    out.push(row[start..].trim());
    out
}

fn align_row(row: &str) -> Option<Vec<Align>> {
    let cs = cells(row);
    cs.iter()
        .all(|c| !c.is_empty() && c.chars().all(|ch| matches!(ch, '-' | ':')) && c.contains('-'))
        .then(|| {
            cs.iter()
                .map(|c| match (c.starts_with(':'), c.ends_with(':')) {
                    (true, true) => Align::Center,
                    (false, true) => Align::Right,
                    _ => Align::Left,
                })
                .collect()
        })
}

/// Columns a painted string takes on screen: escape codes take none.
fn screen_len(s: &str) -> usize {
    let mut n = 0;
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '\x1b' {
            chars.by_ref().find(|&d| d == 'm');
        } else {
            n += 1;
        }
    }
    n
}

/// A table drawn with full borders, cells word-wrapped to fit `width`.
///
/// Columns start at their natural width; while the table is too wide the
/// widest column gives up a column at a time, never below its longest word.
/// A table that still does not fit overflows — cutting a flag name in two
/// would be worse than a scroll.
fn render_table(lines: &[&str], width: usize, style: Style) -> String {
    // 第二行是对齐行才算表格;不是的话(手册里没有这种写法)原样输出而不是猜
    let Some(aligns) = lines.get(1).and_then(|r| align_row(r)) else {
        return lines.join("\n");
    };
    let bold = Look { bold: true, ..Look::default() };
    let rows: Vec<Vec<Vec<Span>>> = lines
        .iter()
        .enumerate()
        .filter(|(k, _)| *k != 1)
        .map(|(k, row)| {
            cells(row)
                .iter()
                .map(|c| {
                    let spans = inline(c);
                    if k > 0 {
                        return spans;
                    }
                    spans.into_iter().map(|s| Span { look: s.look.or(bold), ..s }).collect()
                })
                .collect()
        })
        .collect();
    let ncols = rows.iter().map(Vec::len).max().unwrap_or(0).max(aligns.len());

    let longest_word = |spans: &[Span]| {
        wrap(spans, 1, style).iter().map(|l| screen_len(l)).max().unwrap_or(0)
    };
    let mut natural = vec![1; ncols];
    let mut floor = vec![1; ncols];
    for row in &rows {
        for (j, cell) in row.iter().enumerate() {
            natural[j] = natural[j].max(visible_len(cell, style));
            floor[j] = floor[j].max(longest_word(cell));
        }
    }
    let avail = width.saturating_sub(3 * ncols + 1);
    let mut widths = natural.clone();
    while widths.iter().sum::<usize>() > avail {
        // 最宽且还能让的那一列让出一格
        let Some(j) = (0..ncols).filter(|&j| widths[j] > floor[j]).max_by_key(|&j| widths[j]) else {
            break;
        };
        widths[j] -= 1;
    }

    let (h, v, corners) = if style.unicode {
        ('─', '│', [['┌', '┬', '┐'], ['├', '┼', '┤'], ['└', '┴', '┘']])
    } else {
        ('-', '|', [['+'; 3]; 3])
    };
    let border = |[l, m, r]: [char; 3]| {
        let mut s = String::from(l);
        for (j, w) in widths.iter().enumerate() {
            s.extend(std::iter::repeat_n(h, w + 2));
            s.push(if j + 1 == ncols { r } else { m });
        }
        s
    };

    let wrapped: Vec<Vec<Vec<String>>> = rows
        .iter()
        .map(|row| (0..ncols).map(|j| row.get(j).map_or(vec![], |c| wrap(c, widths[j], style))).collect())
        .collect();
    // 有单元格折了行,行与行之间就要画线,否则分不清哪几行属于同一条记录
    let ruled = wrapped.iter().any(|row| row.iter().any(|c| c.len() > 1));

    let mut out = vec![border(corners[0])];
    for (k, row) in wrapped.iter().enumerate() {
        if k == 1 || (k > 1 && ruled) {
            out.push(border(corners[1]));
        }
        let height = row.iter().map(Vec::len).max().unwrap_or(1).max(1);
        for i in 0..height {
            let mut line = String::from(v);
            for (j, cell) in row.iter().enumerate() {
                let text = cell.get(i).map_or("", String::as_str);
                let gap = widths[j].saturating_sub(screen_len(text));
                let (l, r) = match aligns.get(j).copied().unwrap_or(Align::Left) {
                    Align::Left => (0, gap),
                    Align::Right => (gap, 0),
                    Align::Center => (gap / 2, gap - gap / 2),
                };
                write!(line, " {}{text}{} {v}", " ".repeat(l), " ".repeat(r)).unwrap();
            }
            out.push(line);
        }
    }
    out.push(border(corners[2]));
    out.join("\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    const PLAIN: Style = Style { ansi: false, unicode: false };
    const RICH: Style = Style { ansi: true, unicode: true };

    fn look(f: impl FnOnce(&mut Look)) -> Look {
        let mut l = Look::default();
        f(&mut l);
        l
    }

    #[test]
    fn inline_marks_and_their_boundaries() {
        let s = inline("a **b** `c` *d*");
        let texts: Vec<&str> = s.iter().map(|s| s.text.as_str()).collect();
        assert_eq!(texts, ["a ", "b", " ", "c", " ", "d"]);
        assert!(s[1].look.bold);
        assert!(s[3].look.code);
        assert!(s[5].look.italic);
    }

    #[test]
    fn code_is_opaque() {
        // 代码里的 ** 与 $ 不是标记
        let s = inline("`**x** $y$`");
        assert_eq!(s, [Span { text: "**x** $y$".into(), look: look(|l| l.code = true) }]);
    }

    #[test]
    fn double_backtick_code_may_contain_a_backtick() {
        let s = inline("``a`b``");
        assert_eq!(s[0].text, "a`b");
        assert!(s[0].look.code);
    }

    #[test]
    fn bold_can_contain_code() {
        let s = inline("**see `x`**");
        assert!(s[0].look.bold && !s[0].look.code);
        assert!(s[1].look.bold && s[1].look.code, "粗体里的代码应同时带两种标记");
    }

    #[test]
    fn link_text_is_parsed_and_the_url_kept() {
        let s = inline("[a **b**](x.md)");
        assert_eq!(s[0].text, "a ");
        assert!(s[0].look.link);
        assert!(s[1].look.link && s[1].look.bold);
        assert_eq!(s[2].text, " (x.md)");
        assert!(s[2].look.dim);
    }

    #[test]
    fn unmatched_marks_stay_literal() {
        let text: String = inline("a ** b * c ` d [e] f").iter().map(|s| s.text.as_str()).collect();
        assert_eq!(text, "a ** b * c ` d [e] f");
    }

    #[test]
    fn math_is_opaque_and_needs_a_closing_dollar() {
        let s = inline("$a_b * c$ and $ 5");
        assert_eq!(s[0].text, "a_b * c");
        assert!(s[0].look.math);
        assert_eq!(s[1].text, " and $ 5");
    }

    #[test]
    fn escapes_drop_the_backslash() {
        let s = inline(r"a \| b \* c");
        assert_eq!(s[0].text, "a | b * c");
    }

    #[test]
    fn wrap_respects_width_and_keeps_punctuation_on_its_word() {
        let lines = wrap(&inline("aaa `bb`, cc dd"), 8, PLAIN);
        assert_eq!(lines, ["aaa", "`bb`, cc", "dd"]);
    }

    #[test]
    fn code_that_fits_stays_whole() {
        let lines = wrap(&inline("x `a b` y"), 5, PLAIN);
        assert_eq!(lines, ["x", "`a b`", "y"]);
    }

    #[test]
    fn code_longer_than_the_line_breaks_at_its_spaces() {
        // 朴素样式:反引号只在首尾各一个
        let lines = wrap(&inline("x `a very long code span` y"), 10, PLAIN);
        assert_eq!(lines, ["x `a very", "long code", "span` y"]);
        // 着色样式:每个碎片都仍是代码
        let rich = wrap(&inline("`aaaa bbbb cccc`"), 9, RICH);
        assert_eq!(rich.len(), 2);
        assert!(rich.iter().all(|l| l.starts_with("\x1b[36m")), "{rich:?}");
    }

    #[test]
    fn wrap_counts_visible_columns_not_escape_codes() {
        // ANSI 下粗体多出的转义码不占列宽,不该提前折行
        let lines = wrap(&inline("**aaaa** bbbb"), 9, RICH);
        assert_eq!(lines.len(), 1);
    }

    #[test]
    fn blocks_are_told_apart_by_how_lines_start() {
        let md = "# T\n\npara\nline\n\n- a\n  cont\n- b\n\n> q\n\n| x |\n|---|\n\n---\n\n```rust\nfn\n```\n";
        let b = split_blocks(md);
        assert_eq!(b[0], Block::Heading { level: 1, text: "T" });
        assert_eq!(b[1], Block::Para { indent: 0, text: "para line".into() });
        let Block::List(items) = &b[2] else { panic!("应为列表: {:?}", b[2]) };
        assert_eq!(items[0].text, "a cont");
        assert_eq!(items[1].marker, "-");
        assert_eq!(b[3], Block::Quote(vec!["q"]));
        assert_eq!(b[4], Block::Table(vec!["| x |", "|---|"]));
        assert_eq!(b[5], Block::Rule);
        assert_eq!(b[6], Block::Code { indent: 0, lines: vec!["fn"] });
    }

    #[test]
    fn a_fence_inside_a_list_item_keeps_its_indent() {
        let md = "- item\n\n  ```\n  x\n    y\n  ```\n\n  after\n";
        let b = split_blocks(md);
        assert_eq!(b[1], Block::Code { indent: 2, lines: vec!["x", "  y"] });
        assert_eq!(b[2], Block::Para { indent: 2, text: "after".into() });
    }

    #[test]
    fn block_math_is_raw() {
        let md = "$$a$$\n\n$$\nb \\\\\n$$\n";
        let b = split_blocks(md);
        assert_eq!(b[0], Block::Raw(vec!["$$a$$"]));
        assert_eq!(b[1], Block::Raw(vec!["$$", "b \\\\", "$$"]));
    }

    #[test]
    fn numbered_items_keep_their_numbers() {
        let out = render("1. one\n2. two\n", 80, RICH);
        assert!(out.starts_with("1. one\n2. two"));
    }

    #[test]
    fn table_cells_split_at_unescaped_bars_only() {
        assert_eq!(cells(r"| a \| b | `c` |"), [r"a \| b", "`c`"]);
        assert_eq!(cells("|  | x |"), ["", "x"]);
    }

    #[test]
    fn table_is_boxed_and_aligned() {
        let out = render("| k | v |\n|---|:-:|\n| `a` | 1 |\n| bb | 22 |\n", 80, PLAIN);
        let want = "\
+-----+----+
| k   | v  |
+-----+----+
| `a` | 1  |
| bb  | 22 |
+-----+----+
";
        assert_eq!(out, want);
    }

    #[test]
    fn escaped_bar_renders_as_a_bar() {
        let out = render("| x |\n|---|\n| ENERGY\\| Total |\n", 80, PLAIN);
        assert!(out.contains("| ENERGY| Total |"), "{out}");
    }

    #[test]
    fn a_narrow_terminal_wraps_the_widest_column() {
        let md = "| flag | meaning |\n|---|---|\n| `-o` | one two three four five six |\n| `-s` | x |\n";
        let out = render(md, 24, PLAIN);
        for l in out.lines() {
            assert!(l.len() <= 24, "超出终端宽度: {l:?}\n{out}");
        }
        // 折了行就要有行间分隔线
        assert_eq!(out.lines().filter(|l| l.starts_with('+')).count(), 4, "{out}");
    }

    #[test]
    fn a_word_longer_than_the_room_overflows_instead_of_being_cut() {
        let out = render("| a |\n|---|\n| `--a-very-long-flag` |\n", 10, PLAIN);
        assert!(out.contains("`--a-very-long-flag`"), "{out}");
    }

    #[test]
    fn rich_table_aligns_on_screen_columns_not_bytes() {
        // 粗体表头的转义码与 Å 的多字节都不该挤歪边框
        let out = render("| Å | b |\n|---|---|\n| xx | y |\n", 80, RICH);
        let widths: Vec<usize> = out.lines().map(screen_len).collect();
        assert!(widths.windows(2).all(|w| w[0] == w[1]), "{out}");
    }

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
        ];
        for (tex, want) in cases {
            assert_eq!(latex(tex), want, "输入: {tex}");
        }
    }

    #[test]
    fn unknown_commands_stay_as_written() {
        assert_eq!(latex(r"|\overrightarrow{BA}| < r"), r"|\overrightarrow{BA}| < r");
    }

    #[test]
    fn math_is_converted_only_where_unicode_is_trusted() {
        assert_eq!(render("a $\\alpha_1$ b\n", 80, RICH), "a α₁ b\n");
        assert_eq!(render("a $\\alpha_1$ b\n", 80, PLAIN), "a $\\alpha_1$ b\n");
    }

    #[test]
    fn plain_style_is_pure_ascii() {
        let out = render("# T\n\n- **a** `b` $\\alpha$\n\n> q\n\n| a |\n|---|\n| b |\n\n---\n", 80, PLAIN);
        assert!(out.is_ascii(), "Windows 路径不能出现非 ASCII 字符:\n{out}");
        assert!(!out.contains('\x1b'));
    }
}
