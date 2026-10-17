"""Quantum ESPRESSO：pw.x 输入的解析、编辑与写回。

namelist 的解析与写出调 Rust（``ferro.parse_qe_input_parts`` / ``ferro.write_namelists``），
与 ``ferro.read`` 读 QE 结构是同一份代码，语法只实现一次。namelist 内的注释不保留；
卡片部分按行原样保留（含注释），只有被替换的卡片才重写。
"""

import re

import ferro
from ferro import config as _config
from ferro.jobs.base import Generator, InputSet, JobError, kmesh, multiplicity

RY_TIME_FS = 0.048378   # Rydberg 原子时间单位（INPUT_PW dt：1 a.u. = 4.8378e-17 s）

# pw.x 的全部卡片名（INPUT_PW）。卡片在 Rust reader 里另有按名查找（readers/qe.rs
# find_card），那里只为读结构定位三张卡片，这里为编辑切出全部卡片，用途不同，不合并
CARDS = (
    "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "ADDITIONAL_K_POINTS",
    "CELL_PARAMETERS", "CONSTRAINTS", "OCCUPATIONS", "ATOMIC_VELOCITIES",
    "ATOMIC_FORCES", "SOLVENTS", "HUBBARD",
)


class QeInputError(ValueError):
    pass


def format_value(value):
    """Python 值 → Fortran namelist 值：bool → ``.true.``，str 加单引号。"""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        if "'" in value:
            raise QeInputError(f"a namelist string cannot hold a single quote: {value!r}")
        return f"'{value}'"
    if isinstance(value, (int, float)):
        return repr(value)
    raise QeInputError(f"cannot write {value!r} as a namelist value")


def unquote(raw):
    """原文值去掉引号（``'scf'`` → ``scf``）；不带引号的原样返回。"""
    raw = raw.strip()
    if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in "'\"":
        return raw[1:-1]
    return raw


def _split_cards(text):
    """卡片部分 → (卡片之前的行, [[卡片名, 卡片行...]])。"""
    head, cards = [], []
    for line in text.splitlines():
        word = line.split()[0].upper() if line.split() else ""
        if word in CARDS:
            cards.append([word, line])
        elif cards:
            cards[-1].append(line)
        else:
            head.append(line)
    return head, cards


class QeInput:
    """pw.x 输入。namelist 名与键都不分大小写，取值是原文（``'scf'``、``.true.``）。"""

    def __init__(self, namelists, head, cards):
        self.namelists = namelists  # {名小写: {键小写: 原文值}}，保持顺序
        self.head = head            # 第一张卡片之前、namelist 之外的行（多为注释）
        self._cards = cards

    @classmethod
    def from_str(cls, text):
        try:
            parsed, rest = ferro.parse_qe_input_parts(text)
        except RuntimeError as e:
            raise QeInputError(str(e)) from None
        namelists = {}
        for name, entries in parsed:
            if name in namelists:
                raise QeInputError(f"namelist &{name.upper()} appears twice")
            namelists[name] = dict(entries)
        head, cards = _split_cards(rest)
        return cls(namelists, head, cards)

    def get(self, namelist, key, default=None):
        return self.namelists.get(namelist.lower(), {}).get(key.lower(), default)

    def set(self, namelist, key, value):
        """设一项；namelist 不存在就在末尾新建。``value`` 按 Python 类型写成 Fortran。"""
        self.namelists.setdefault(namelist.lower(), {})[key.lower()] = format_value(value)

    def delete(self, namelist, key):
        self.namelists.get(namelist.lower(), {}).pop(key.lower(), None)

    def card(self, name):
        """一张卡片的原文（首行起，含尾随空行与注释）；没有返回 None。"""
        for c in self._cards:
            if c[0] == name.upper():
                return "\n".join(c[1:])
        return None

    def set_card(self, name, text):
        """替换一张卡片（``text`` 从卡片首行起）；没有就追加到末尾。"""
        name = name.upper()
        lines = text.rstrip("\n").splitlines() + [""]
        if lines[0].split()[0].upper() != name:
            raise QeInputError(f"card text must start with {name}")
        for c in self._cards:
            if c[0] == name:
                c[1:] = lines
                return
        self._cards.append([name] + lines)

    def remove_card(self, name):
        self._cards = [c for c in self._cards if c[0] != name.upper()]

    def to_str(self):
        out = ferro.write_namelists([(n, list(e.items())) for n, e in self.namelists.items()])
        body = self.head + [line for c in self._cards for line in c[1:]]
        while body and not body[0].strip():
            body.pop(0)
        return out + "\n" + "\n".join(body).rstrip("\n") + "\n"


# ── 生成器 ───────────────────────────────────────────────────────────────────

def _fortran_float(text):
    """``6.000000000000000E+000`` / ``1.0d-8`` → float；空串为 0。"""
    return float(text.strip().lower().replace("d", "e")) if text.strip() else 0.0


def upf_header(path):
    """UPF 头部 → (元素, 价电子数, wfc 建议截断, rho 建议截断)，截断单位 Ry，未给为 0。

    UPF v2 是 ``<PP_HEADER element="O" z_valence="6.0" wfc_cutoff=… rho_cutoff=…/>``；
    v1 是 ``<PP_HEADER>`` 块里按行的「值  说明」。只读文件开头，UPF 可达数 MB。
    """
    with open(path, encoding="utf-8", errors="replace") as f:
        head = f.read(65536)
    m = re.search(r"<PP_HEADER\b([^>]*)>", head, re.S)
    if m and "=" in m.group(1):
        attrs = {k: v.strip() for k, v in re.findall(r'(\w+)\s*=\s*"([^"]*)"', m.group(1))}
        numbers = [_fortran_float(attrs.get(k, "")) for k in ("z_valence", "wfc_cutoff", "rho_cutoff")]
        return attrs.get("element", ""), *numbers
    block = head.split("<PP_HEADER>", 1)[1].split("</PP_HEADER>", 1)[0] if "<PP_HEADER>" in head else ""
    element, z, wfc, rho = "", 0.0, 0.0, 0.0
    for line in block.splitlines():
        words = line.split()
        text = line.lower()
        if text.rstrip().endswith("element") and words:
            element = words[0]
        elif "z valence" in text and words:
            z = _fortran_float(words[0])
        elif "suggested cutoff" in text and len(words) >= 2:
            wfc, rho = (_fortran_float(w) for w in words[:2])
    return element, z, wfc, rho


class QeGenerator(Generator):
    program = "qe"
    template_name = "pw.in"

    def pseudos(self, elements):
        """元素 → (文件名, 价电子数, wfc, rho)。索引 pseudo_dir 里全部 UPF，按头部的元素归类。"""
        root = _config.require(self.config, "qe", "pseudo_dir")
        if not root.is_dir():
            raise JobError(f"[qe] pseudo_dir {root} is not a directory")
        index = {}
        for path in sorted(p for p in root.iterdir() if p.suffix.lower() == ".upf"):
            element, z, wfc, rho = upf_header(path)
            element = element or re.match(r"[A-Za-z][a-z]?", path.name).group(0).capitalize()
            index.setdefault(element, []).append((path.name, z, wfc, rho))
        chosen = {}
        for e in elements:
            wanted = self.config["qe"]["pseudos"].get(e)
            found = index.get(e, [])
            if wanted is not None:
                match = [x for x in found if x[0] == wanted]
                if not match:
                    raise JobError(f"[qe.pseudos] {e} = {wanted}: no such UPF for {e} in {root}")
                chosen[e] = match[0]
            elif len(found) == 1:
                chosen[e] = found[0]
            else:
                raise JobError(
                    f"{len(found) or 'no'} UPF files for {e} in {root}"
                    + (f" ({', '.join(x[0] for x in found)}): choose one in [qe.pseudos]" if found else ""))
            if chosen[e][1] <= 0:
                raise JobError(f"{chosen[e][0]}: no z_valence in its PP_HEADER")
        return root, chosen

    def build(self, one, options, template_text):
        qe = QeInput.from_str(template_text)
        qe.head = []   # 模板头注释是写给改模板的人看的，不带进生成的 pw.in
        notes = []
        elements = one.elements()
        symbols = one.symbols(0)
        root, chosen = self.pseudos(elements)
        notes.append("UPF: " + ", ".join(f"{e} {p[0]}" for e, p in chosen.items()))

        qe.set("control", "pseudo_dir", str(root) + "/")
        qe.set("system", "nat", len(symbols))
        qe.set("system", "ntyp", len(elements))
        # 胞总由 CELL_PARAMETERS 给：ibrav ≠ 0 时 pw.x 不用这张卡片，模板的 ibrav / celldm
        # 会静默定出另一个胞；结构属于程序改写的量
        qe.set("system", "ibrav", 0)
        for key in list(qe.namelists.get("system", {})):
            if re.fullmatch(r"celldm\(\d\)|a|b|c|cosab|cosac|cosbc", key):
                qe.delete("system", key)
        rho_written = qe.get("system", "ecutrho") is not None
        for key, i in (("ecutwfc", 2), ("ecutrho", 3)):
            if qe.get("system", key) is None:
                suggested = max(p[i] for p in chosen.values())
                if suggested <= 0:
                    raise JobError(f"the UPF headers suggest no {key}: write {key} in the template")
                qe.set("system", key, suggested)
                notes.append(f"{key} = {suggested:g} Ry (largest UPF suggestion)")
        # pw.x set_cutoff（7.3.1）：ecutrho ≤ ecutwfc 才 errore，< 4 × ecutwfc 只 infomsg 照跑，
        # 但密度 |ψ|² 需要 4 倍截断。UPF 建议值不够就提到 4 倍（同 pw.x 不写 ecutrho 时的默认）；
        # 模板写的照用，只提示
        wfc, rho = (_fortran_float(qe.get("system", k)) for k in ("ecutwfc", "ecutrho"))
        if rho <= wfc:
            raise JobError(f"ecutrho must exceed ecutwfc (pw.x stops): {rho:g} <= {wfc:g} Ry")
        if rho < 4 * wfc:
            if rho_written:
                notes.append(f"ecutrho {rho:g} Ry is below 4 x ecutwfc = {4 * wfc:g} Ry (pw.x only warns)")
            else:
                qe.set("system", "ecutrho", 4 * wfc)
                notes[-1] = f"ecutrho = {4 * wfc:g} Ry (4 x ecutwfc; the largest UPF suggestion is {rho:g})"

        n_valence = sum(symbols.count(e) * chosen[e][1] for e in elements)
        mult = multiplicity(n_valence, options.charge, options.multiplicity)
        qe.set("system", "tot_charge", options.charge)
        if mult > 1:
            qe.set("system", "nspin", 2)
            qe.set("system", "tot_magnetization", mult - 1)
        notes.append(f"{n_valence - options.charge:g} electrons, multiplicity {mult}")

        cards = ferro.qe_structure_cards(one, {e: p[0] for e, p in chosen.items()})
        for block in cards.strip().split("\n\n"):
            qe.set_card(block.split()[0], block)
        if options.kspacing is not None or qe.card("K_POINTS") is None:
            spacing = options.kspacing if options.kspacing is not None else 0.5
            mesh = kmesh(one.cell(0), spacing)
            qe.set_card("K_POINTS", "K_POINTS automatic\n  " + " ".join(map(str, mesh)) + "  0 0 0")
            notes.append(f"kspacing {spacing:g} -> mesh {' x '.join(map(str, mesh))}")
        else:
            notes.append("K_POINTS: kept as written in the template")
        if options.fine_grid is not None:
            notes.append("--fine-grid ignored: pw.x has no separate fine grid (the dense grid follows ecutrho)")

        if options.task == "opt" and options.vary_cell:
            qe.set("control", "calculation", "vc-relax")
            if "cell" not in qe.namelists:
                qe.set("cell", "cell_dynamics", "bfgs")
                qe.set("cell", "press_conv_thr", 0.5)
        if options.task == "aimd":
            if options.ensemble not in (None, "nvt"):
                raise JobError(f"QE jobs are NVT only (vc-md has no thermostat); --ensemble {options.ensemble} "
                               "is a CP2K option")
            if options.temperature is not None:
                qe.set("ions", "tempw", options.temperature)
            if options.steps is not None:
                qe.set("control", "nstep", options.steps)
            if options.timestep is not None:
                qe.set("control", "dt", round(options.timestep / RY_TIME_FS, 4))
                # SVR 的时间常数 nraise 以步计，改 dt 热浴跟着变快或变慢
                if qe.get("ions", "ion_temperature", "").strip("'\"").lower() == "svr":
                    nraise = int(qe.get("ions", "nraise", "1"))
                    notes.append(f"SVR time constant nraise = {nraise} steps = {nraise * options.timestep:g} fs "
                                 "at this timestep (nraise counts steps, not time)")
        notes.append("parallel settings are not in pw.in: add your machine's to the run command (-nk, -nd)")
        return InputSet({"pw.in": qe.to_str()}, notes=notes)
