"""CP2K：输入节树的解析、编辑与写回（结构参照 ASE ``calculators/cp2k.py`` 的 ``InputSection``）。

节名、关键字都不分大小写（CP2K 如此），保存时保留原文写法。整行注释、空行、
行尾注释都挂在节树上随之写回；预处理指令（``@SET``、``@INCLUDE`` 等）原样保留、
不展开。写回统一两格缩进，所以不是逐字节还原，但再解析得到同一棵树。

后半是 CP2K 任务的生成器：&KIND 的基组与势在运行时从 CP2K 的数据文件里查。
"""

import re

import ferro
from ferro import config as _config
from ferro.jobs.base import DEFAULT_KSPACING, Generator, InputSet, JobError, kmesh, multiplicity


class Cp2kInputError(ValueError):
    pass


def _split_comment(line):
    """引号外第一个 ``#`` 或 ``!`` 起是注释；返回 (代码, 注释)，注释含起始符。

    ferro-io's CP2K reader (readers/cp2k.rs `parse_cp2k`) strips comments its own way: `!` before
    `#`, quotes ignored. Not merged, the reason is written there.
    """
    quote = None
    for i, c in enumerate(line):
        if quote:
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c in "#!":
            return line[:i].rstrip(), line[i:]
    return line.rstrip(), ""


# 与 vasp.py 的 format_value 逐字相同，按 R6 两处不合：INCAR 与 CP2K 输入的值语法各自演化
def format_value(value):
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, (list, tuple)):
        return " ".join(format_value(v) for v in value)
    return str(value)


class Keyword:
    def __init__(self, name, value="", comment=""):
        self.name, self.value, self.comment = name, value, comment

    def to_lines(self, indent):
        text = f"{indent}{self.name} {self.value}".rstrip()
        return [f"{text}  {self.comment}" if self.comment else text]


class Line:
    """整行注释、空行或预处理指令：原样写回（去掉原缩进）。"""

    def __init__(self, text):
        self.text = text

    def to_lines(self, indent):
        return [f"{indent}{self.text}" if self.text else ""]


class Section:
    """一节：``&NAME param`` …… ``&END NAME``。根节点 ``name`` 为 None。"""

    def __init__(self, name=None, param="", comment=""):
        self.name, self.param, self.comment = name, param, comment
        self.items = []

    @classmethod
    def from_str(cls, text):
        root = cls()
        stack = [root]
        for n, raw in enumerate(text.splitlines(), 1):
            code, comment = _split_comment(raw.strip())
            if not code:
                stack[-1].items.append(Line(comment))
                continue
            word, rest = (code.split(None, 1) + [""])[:2]
            if word.upper() == "&END":
                if len(stack) == 1:
                    raise Cp2kInputError(f"line {n}: &END without an open section")
                if rest and rest.upper() != stack[-1].name.upper():
                    raise Cp2kInputError(
                        f"line {n}: &END {rest} closes &{stack[-1].name}")
                stack.pop()
            elif word.startswith("&"):
                section = cls(word[1:], rest, comment)
                stack[-1].items.append(section)
                stack.append(section)
            elif word.startswith("@"):
                stack[-1].items.append(Line(raw.strip()))
            else:
                stack[-1].items.append(Keyword(word, rest, comment))
        if len(stack) > 1:
            raise Cp2kInputError(f"section &{stack[-1].name} is not closed with &END")
        return root

    # ── 查找 ───────────────────────────────────────────────────────────────
    def sections(self, name):
        return [s for s in self.items if isinstance(s, Section) and s.name.upper() == name.upper()]

    def section(self, path, create=False):
        """``"FORCE_EVAL/DFT/MGRID"`` → 那一节（同名取第一个）；没有返回 None，``create`` 时新建。"""
        node = self
        for name in path.split("/"):
            found = node.sections(name)
            if found:
                node = found[0]
            elif create:
                new = Section(name.upper())
                node.items.append(new)
                node = new
            else:
                return None
        return node

    def keywords(self, name):
        return [k for k in self.items if isinstance(k, Keyword) and k.name.upper() == name.upper()]

    def get(self, name, default=None):
        found = self.keywords(name)
        return found[0].value if found else default

    # ── 编辑 ───────────────────────────────────────────────────────────────
    def set(self, name, value):
        """设关键字（同名取第一个，其余删去）；没有就加在本节子节之前。"""
        found = self.keywords(name)
        if found:
            found[0].value = format_value(value)
            for extra in found[1:]:
                self.items.remove(extra)
            return
        position = next((i for i, x in enumerate(self.items) if isinstance(x, Section)), len(self.items))
        self.items.insert(position, Keyword(name.upper(), format_value(value)))

    def remove(self, name):
        """删去这一名字的全部关键字与子节。"""
        self.items = [x for x in self.items
                      if not (isinstance(x, (Keyword, Section)) and x.name.upper() == name.upper())]

    def add(self, item):
        self.items.append(item)
        return item

    # ── 写回 ───────────────────────────────────────────────────────────────
    def to_lines(self, indent=""):
        if self.name is None:
            return [line for x in self.items for line in x.to_lines(indent)]
        head = f"{indent}&{self.name} {self.param}".rstrip()
        lines = [f"{head}  {self.comment}" if self.comment else head]
        for x in self.items:
            lines += x.to_lines(indent + "  ")
        lines.append(f"{indent}&END {self.name}")
        return lines

    def to_str(self):
        return "\n".join(self.to_lines()) + "\n"


# ── 生成器 ───────────────────────────────────────────────────────────────────

DEFAULTS = {"basis_family": "DZVP-MOLOPT-SR-GTH", "potential_family": "GTH-PBE",
            "ae_basis_family": "pob-TZVP-rev2"}


def data_entries(path):
    """CP2K 基组 / 势文件 → {元素: [(名字与别名, 下一行的数字)]}。

    条目头是 ``O DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q6``：元素、名字、别名。势文件
    条目的下一行是各角动量的电子数（``2 4``），求和即这个势的电子数。
    """
    entries = {}
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for i, line in enumerate(lines):
        words = line.split("#")[0].split()
        if len(words) >= 2 and re.fullmatch(r"[A-Z][a-z]?", words[0]) and not re.match(r"[-+\d.]", words[1]):
            following = lines[i + 1].split("#")[0].split() if i + 1 < len(lines) else []
            numbers = [int(w) for w in following if re.fullmatch(r"\d+", w)]
            entries.setdefault(words[0], []).append((words[1:], numbers))
    return entries


class Cp2kGenerator(Generator):
    program = "cp2k"
    template_name = "cp2k.inp"

    def _files(self, dft, keyword):
        names = [k.value.split()[0] for k in dft.keywords(keyword) if k.value.split()]
        if not names:
            raise JobError(f"the template's &DFT has no {keyword}")
        root = _config.require(self.config, "cp2k", "data_dir")
        found = {}
        for name in names:
            path = root / name
            if not path.is_file():
                raise JobError(f"{keyword} {name}: {path} does not exist ([cp2k] data_dir)")
            for element, items in data_entries(path).items():
                found.setdefault(element, []).extend(items)
        return found, names

    def kinds(self, elements, dft, all_electron):
        """元素 → (基组名, 势名, 电子数)。不猜：找不到、或有歧义就报错并列出候选。"""
        cfg = self.config["cp2k"]
        family = cfg.get("ae_basis_family" if all_electron else "basis_family",
                         DEFAULTS["ae_basis_family" if all_electron else "basis_family"])
        potential_family = cfg.get("potential_family", DEFAULTS["potential_family"])
        basis_db, basis_files = self._files(dft, "BASIS_SET_FILE_NAME")
        potential_db, potential_files = self._files(dft, "POTENTIAL_FILE_NAME")
        out = {}
        for e in elements:
            override = cfg["kinds"].get(e, {})
            names = [n for names, _ in basis_db.get(e, []) for n in names]
            potentials = {n: numbers for names_, numbers in potential_db.get(e, []) for n in names_}
            basis = override.get("basis")
            if basis is None:
                if all_electron:
                    # 全电子的条目没有 -qN 别名；pob-TZVP-rev2 在 Rb 之后的条目本名相同、
                    # 别名带 -qN，配的是 ECP，不是全电子
                    entries = [n for n, _ in basis_db.get(e, []) if family in n]
                    if any(not any(re.search(r"-q\d+$", x) for x in n) for n in entries):
                        basis = family
                    elif entries:
                        raise JobError(f"{family} for {e} is an ECP basis, not all-electron")
                else:
                    candidates = list(dict.fromkeys(
                        n for n in names if re.fullmatch(rf"{re.escape(family)}-q\d+", n)))
                    if len(candidates) > 1:
                        # 同族多个 -qN（如 Rh 的 q9 与 q17）：基组与势须同 q 配对，只留势族里
                        # 有同 q 势的；仍不止一个或一个不剩就报错，不按文件顺序默取第一个
                        paired = [n for n in candidates
                                  if f"{potential_family}-q{n.rsplit('-q', 1)[1]}" in potentials]
                        if len(paired) != 1:
                            raise JobError(f"{e}: {family} comes as {', '.join(candidates)}, of which "
                                           f"{len(paired)} pair with a {potential_family}-qN potential; "
                                           f"name basis and potential in [cp2k.kinds]")
                        candidates = paired
                    basis = candidates[0] if candidates else None
                if basis is None:
                    raise JobError(f"no {family} basis for {e} in {', '.join(basis_files)}; "
                                   f"available: {', '.join(sorted(set(names))) or 'none'}")
            elif basis not in names:
                raise JobError(f"[cp2k.kinds] {e}: basis {basis} is not in {', '.join(basis_files)}")
            potential = override.get("potential")
            if potential is None:
                if all_electron:
                    potential = "ALL"
                else:
                    q = re.search(r"-q(\d+)$", basis)
                    if q is None:
                        raise JobError(f"basis {basis} of {e} has no -qN: name the potential in [cp2k.kinds]")
                    potential = f"{potential_family}-q{q.group(1)}"
            if potential not in potentials:
                raise JobError(f"no potential {potential} for {e} in {', '.join(potential_files)}; "
                               f"available: {', '.join(sorted(potentials)) or 'none'}")
            if not potentials[potential]:
                # POTENTIAL_UZH 的镧系占位条目下一行是 NA，没有电子数；求和得 0 会静默写错电荷
                raise JobError(f"potential {potential} of {e} in {', '.join(potential_files)} lists no "
                               "electron counts (a placeholder entry); name another in [cp2k.kinds]")
            out[e] =(basis, potential, sum(potentials[potential]))
        return out

    def build(self, one, options, template_text):
        root = Section.from_str(template_text)
        force_eval = root.section("FORCE_EVAL")
        dft = root.section("FORCE_EVAL/DFT")
        if force_eval is None or dft is None:
            raise JobError("the template has no &FORCE_EVAL / &DFT")
        method = (root.section("FORCE_EVAL/DFT/QS") or Section()).get("METHOD", "GPW").split()[0].upper()
        # GAPW_XC 不是全电子（CP2K 对它配 ALL 直接 CPABORT），配赝势时用户以有限差分
        # 验证过力有错：不支持
        if method == "GAPW_XC":
            raise JobError("&QS METHOD GAPW_XC is not supported (its forces fail a finite-difference "
                           "check); use GPW or GAPW")
        all_electron = method == "GAPW"
        notes = []

        elements = one.elements()
        symbols = one.symbols(0)
        kinds = self.kinds(elements, dft, all_electron)
        n_electrons = sum(symbols.count(e) * kinds[e][2] for e in elements)
        mult = multiplicity(n_electrons, options.charge, options.multiplicity)
        dft.set("CHARGE", options.charge)
        dft.set("MULTIPLICITY", mult)
        if mult > 1:
            dft.set("UKS", "T")
        notes.append(f"{method}: " + ", ".join(f"{e} {b} / {p}" for e, (b, p, _) in kinds.items()))
        notes.append(f"{n_electrons - options.charge} electrons, multiplicity {mult}")

        # 结构：&CELL / &COORD 用 Rust 写的两段，&KIND 每元素一个
        subsys = force_eval.section("SUBSYS", create=True)
        for name in ("CELL", "COORD", "KIND", "TOPOLOGY"):
            subsys.remove(name)
        for section in Section.from_str(ferro.cp2k_structure_sections(one)).items:
            subsys.add(section)
        for e, (basis, potential, _) in kinds.items():
            kind = subsys.add(Section("KIND", e))
            kind.set("BASIS_SET", basis)
            kind.set("POTENTIAL", potential)

        if options.fine_grid is not None:
            # CP2K 只有开关（倍数是它自己的 FINE_GRID_FACTOR，2026.x 起才有、默认 4），
            # --fine-grid 的数值在这里不用
            dft.section("XC/XC_GRID", create=True).set("USE_FINER_GRID", "T")
            xc = dft.section("XC/XC_FUNCTIONAL")
            names = [xc.param] + [s.name for s in xc.items if isinstance(s, Section)] if xc else []
            mgga = any(re.search(r"MGGA|SCAN|TPSS", n or "", re.I) for n in names)
            notes.append("USE_FINER_GRID on (the --fine-grid factor is not used by CP2K)"
                         + ("" if mgga else "; not a meta-GGA, the finer grid is probably unnecessary"))
        self._kpoints(one, options, dft, notes)
        self._motion(root, force_eval, options, notes)
        notes.append("parallel settings are not in cp2k.inp: set MPI ranks and OMP_NUM_THREADS in the run command")
        return InputSet({"cp2k.inp": root.to_str()}, notes=notes)

    def _kpoints(self, one, options, dft, notes):
        scf = dft.section("SCF", create=True)
        has_ot = bool(scf.sections("OT"))
        if has_ot and scf.section("SMEAR") is not None:
            raise JobError("the template smears occupations with OT; CP2K needs diagonalisation for &SMEAR")
        if scf.section("SMEAR") is not None and scf.get("ADDED_MOS") is None:
            raise JobError("the template's &SMEAR needs ADDED_MOS in &SCF (empty orbitals)")
        if options.kspacing is None and dft.sections("KPOINTS"):
            notes.append("k-points: kept as written in the template")
            return
        spacing = options.kspacing if options.kspacing is not None else DEFAULT_KSPACING
        mesh = kmesh(one.cell(0), spacing)
        dft.remove("KPOINTS")
        if mesh == [1, 1, 1]:
            notes.append(f"kspacing {spacing:g} -> Gamma only" + (", OT" if has_ot else ""))
            return
        kp = dft.add(Section("KPOINTS"))
        kp.set("SCHEME", ["MONKHORST-PACK", *mesh])
        if has_ot:
            # OT 只能 Γ 点（qs_scf_initialization.F 拒绝 OT + k 点）：换成用户脚本的对角化分支
            scf.remove("OT")
            scf.remove("OUTER_SCF")
            scf.set("MAX_SCF", 128)
            scf.add(Section("DIAGONALIZATION")).set("ALGORITHM", "STANDARD")
            mixing = scf.add(Section("MIXING"))
            mixing.set("METHOD", "BROYDEN_MIXING")
            mixing.set("ALPHA", 0.4)
            mixing.set("NBROYDEN", 10)
            notes.append(f"kspacing {spacing:g} -> mesh {' x '.join(map(str, mesh))}: "
                         "OT replaced by diagonalisation + Broyden mixing (OT is Gamma only)")
        else:
            notes.append(f"kspacing {spacing:g} -> mesh {' x '.join(map(str, mesh))}")

    def _motion(self, root, force_eval, options, notes):
        glob = root.section("GLOBAL", create=True)
        if options.task == "opt" and options.vary_cell:
            glob.set("RUN_TYPE", "CELL_OPT")
            root.section("MOTION/CELL_OPT", create=True)
            force_eval.set("STRESS_TENSOR", "ANALYTICAL")
        if options.task != "aimd":
            return
        md = root.section("MOTION/MD", create=True)
        if options.temperature is not None:
            md.set("TEMPERATURE", f"{options.temperature:g}")
        if options.steps is not None:
            md.set("STEPS", options.steps)
        if options.timestep is not None:
            md.set("TIMESTEP", f"{options.timestep:g}")
        if options.ensemble is not None:
            md.set("ENSEMBLE", options.ensemble.upper())
        ensemble = md.get("ENSEMBLE", "NVE").split()[0].upper()
        if ensemble.startswith("NPT"):
            # NPT 要解析应力（md_run.F 不给就报错，审查 D-M4）与一个压浴
            force_eval.set("STRESS_TENSOR", "ANALYTICAL")
            if md.section("BAROSTAT") is None:
                barostat = md.add(Section("BAROSTAT"))
                barostat.set("PRESSURE", 1.01325)
                barostat.set("TIMECON", 1000)
                notes.append(f"{ensemble}: &BAROSTAT added (1.01325 bar, TIMECON 1000 fs)")
