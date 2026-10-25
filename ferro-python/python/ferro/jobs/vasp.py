"""VASP：INCAR 的解析、编辑与写回，以及 VASP 任务的生成器。

模板就是用户跑通过的 INCAR，所以写回要尽量保持原文：没被改的行原样输出（含注释、
对齐、空行），改过的键只换值、行尾注释留着，新加的键追加在末尾。
"""

import math
import re

import ferro
from ferro import config as _config
from ferro.jobs.base import DEFAULT_KSPACING, Generator, InputSet, JobError, kmesh, multiplicity

# VASP 官方推荐的 PBE POTCAR（vasp.at wiki「Available pseudopotentials」，potpaw.64，
# 2026-10-16 查）。表里没有的元素用元素名本身；配置 [vasp.potcar] 逐元素覆盖
RECOMMENDED_POTCAR = {
    "Li": "Li_sv", "Na": "Na_pv", "K": "K_sv", "Ca": "Ca_sv", "Sc": "Sc_sv", "Ti": "Ti_sv",
    "V": "V_sv", "Cr": "Cr_pv", "Mn": "Mn_pv", "Fe": "Fe_pv", "Co": "Co_pv", "Ni": "Ni_pv",
    "Cu": "Cu_pv", "Ga": "Ga_d", "Ge": "Ge_d", "Rb": "Rb_sv", "Sr": "Sr_sv", "Y": "Y_sv",
    "Zr": "Zr_sv", "Nb": "Nb_sv", "Mo": "Mo_sv", "Tc": "Tc_pv", "Ru": "Ru_pv", "Rh": "Rh_pv",
    "In": "In_d", "Sn": "Sn_d", "Cs": "Cs_sv", "Ba": "Ba_sv", "Pr": "Pr_3", "Nd": "Nd_3",
    "Pm": "Pm_3", "Sm": "Sm_3", "Eu": "Eu_2", "Gd": "Gd_3", "Tb": "Tb_3", "Dy": "Dy_3",
    "Ho": "Ho_3", "Er": "Er_2", "Tm": "Tm_3", "Yb": "Yb_2", "Lu": "Lu_3", "Hf": "Hf_pv",
    "Ta": "Ta_pv", "W": "W_sv", "Tl": "Tl_d", "Pb": "Pb_d", "Bi": "Bi_d", "Po": "Po_d",
    "Fr": "Fr_sv", "Ra": "Ra_sv",
}

PARALLEL_KEYS = ("NCORE", "NPAR", "KPAR")
# Low / Medium / High are VASP 4.4 legacy settings (High also raises the default ENCUT to
# 1.3 x ENMAX, which the ENMAX-based fine grid would not see); Single / SingleN are special-purpose
PREC_VALUES = ("Normal", "Accurate")


class IncarError(ValueError):
    pass


def _split_comment(line):
    """引号外第一个 ``#`` 或 ``!`` 起是注释；返回 (代码, 注释)，注释含起始符。"""
    quoted = False
    for i, c in enumerate(line):
        if c == '"':
            quoted = not quoted
        elif c in "#!" and not quoted:
            return line[:i], line[i:]
    return line, ""


def _assignments(code):
    """代码段 → [(键大写, 值的起止下标)]；``;`` 分隔同一行的多个赋值（VASP 允许）。"""
    out, start = [], 0
    for part in re.split(r"(;)", code):
        if part != ";":
            if part.strip():
                if "=" not in part:
                    raise IncarError(f"not a KEY = value assignment: '{part.strip()}'")
                key, _, value = part.partition("=")
                if not key.strip():
                    raise IncarError(f"assignment without a key: '{part.strip()}'")
                v0 = start + len(key) + 1
                v0 += len(value) - len(value.lstrip())
                out.append((key.strip().upper(), v0, start + len(part.rstrip())))
        start += len(part)
    return out


# 与 cp2k.py 的 format_value 逐字相同，按 R6 两处不合：INCAR 与 CP2K 输入的值语法各自演化
def format_value(value):
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, (list, tuple)):
        return " ".join(format_value(v) for v in value)
    return str(value)


class Incar:
    """INCAR。键不分大小写；``incar["ENCUT"]`` 取原文值（字符串），赋值时按类型格式化。"""

    def __init__(self, lines):
        self._lines = lines  # 逻辑行；带引号的多行值并成一行（含换行符）
        self._index()

    @classmethod
    def from_str(cls, text):
        lines, pending = [], None
        for n, line in enumerate(text.splitlines(), 1):
            if pending is not None:
                pending += "\n" + line
                if line.count('"') % 2:
                    lines.append(pending)
                    pending = None
                continue
            if _split_comment(line)[0].count('"') % 2:
                pending = line
            else:
                lines.append(line)
        if pending is not None:
            raise IncarError("INCAR ends inside a quoted value")
        return cls(lines)

    def _index(self):
        self._where = {}
        for i, line in enumerate(self._lines):
            try:
                found = _assignments(_split_comment(line)[0])
            except IncarError as e:
                raise IncarError(f"INCAR line {i + 1}: {e}") from None
            for key, _, _ in found:
                if key in self._where:
                    raise IncarError(f"INCAR sets {key} twice")
                self._where[key] = i

    def keys(self):
        return list(self._where)

    def __contains__(self, key):
        return key.upper() in self._where

    def __getitem__(self, key):
        i = self._where[key.upper()]
        code = _split_comment(self._lines[i])[0]
        for k, a, b in _assignments(code):
            if k == key.upper():
                return code[a:b]

    def get(self, key, default=None):
        return self[key] if key in self else default

    def __setitem__(self, key, value):
        key, text = key.upper(), format_value(value)
        if key not in self._where:
            self._lines.append(f"{key} = {text}")
        else:
            i = self._where[key]
            code, comment = _split_comment(self._lines[i])
            for k, a, b in _assignments(code):
                if k == key:
                    code = code[:a] + text + code[b:]
            self._lines[i] = code + comment
        self._index()

    def __delitem__(self, key):
        key = key.upper()
        i = self._where[key]
        code, comment = _split_comment(self._lines[i])
        parts = [p for p in code.split(";") if p.strip() and p.split("=")[0].strip().upper() != key]
        if parts:
            self._lines[i] = ";".join(parts) + comment
        else:
            del self._lines[i]   # 行尾注释说的就是这个键，一起删
        self._index()

    def comment_out(self, key):
        """``KEY = value`` becomes ``# KEY = value`` where it stood; returns the value.

        Alone on its line, the line is commented as a whole (its trailing comment kept);
        sharing a line, the assignment leaves it and the commented copy follows the line.
        """
        key = key.upper()
        i, value = self._where[key], self[key]
        n, comment = len(self._lines), _split_comment(self._lines[i])[1]
        del self[key]
        if len(self._lines) < n:
            self._lines.insert(i, f"# {key} = {value}" + (f"  {comment}" if comment else ""))
        else:
            self._lines.insert(i + 1, f"# {key} = {value}")
        self._index()
        return value

    def to_str(self):
        return "\n".join(self._lines) + "\n"


# ħ²/2mₑ [eV·Å²] as VASP itself has it: HSQDTM = RYTOEV·AUTOA² in constant.inc, with
# RYTOEV = 13.605826 and AUTOA = 0.529177249 (older than CODATA, kept by VASP on purpose).
# CODATA gives 3.80998, which puts a grid that sits right on a size boundary one step off
# VASP's own NGX. Only that prediction below uses it, so it stays here, not in units.rs
HBAR2_2M = 13.605826 * 0.529177249 ** 2


def _good_size(m):
    """不小于 m 的最小偶数、且只含 2/3/5/7 因子（VASP 自己取的 FFT 尺寸）。"""
    while True:
        if m % 2 == 0:
            r = m
            for p in (2, 3, 5, 7):
                while r % p == 0:
                    r //= p
            if r == 1:
                return m
        m += 1


def fine_grid(cell, encut, factor, prec="Accurate"):
    """Coarse grid NGi (a good size) and fine grid NGiF = factor·NGi.

    The coarse grid covers a sphere of radius 2·G_cut (PREC = Accurate) or 3/2·G_cut (Normal),
    i.e. NGi = 4 (or 3)·G_cut·|a_i|/2π. Accurate follows the same rule as the user's
    ``make_vasp_label.py`` (which also caps the fine-grid spacing; here only the factor is used).
    """
    gcut = math.sqrt(encut / HBAR2_2M)
    k = 4 if prec.strip().upper() == "ACCURATE" else 3
    coarse = [_good_size(round(k * gcut * math.sqrt(sum(x * x for x in a)) / (2 * math.pi))) for a in cell]
    return coarse, [_good_size(math.ceil(factor * n)) for n in coarse]


def _potcar_field(text, name):
    """POTCAR 头部的 ``ZVAL = 4.000`` / ``ENMAX = 400.000``：取第一处。"""
    m = re.search(rf"\b{name}\s*=\s*([-\d.]+)", text)
    if m is None:
        raise JobError(f"no {name} in the POTCAR header")
    return float(m.group(1))


class VaspGenerator(Generator):
    program = "vasp"
    template_name = "INCAR"

    def potcar(self, element):
        """元素 → (POTCAR 文本, 所在路径)。变体：配置 > VASP 推荐 > 元素名本身。"""
        root = _config.require(self.config, "vasp", "potcar_dir")
        name = self.config["vasp"]["potcar"].get(element) or RECOMMENDED_POTCAR.get(element, element)
        path = root / name / "POTCAR"
        if not path.is_file():
            near = sorted(p.name for p in root.glob(f"{element}*") if (p / "POTCAR").is_file()
                          and re.fullmatch(rf"{element}(_.*)?", p.name))
            raise JobError(
                f"no POTCAR for {element}: {path} does not exist"
                + (f"; available: {', '.join(near)} (choose one in [vasp.potcar])" if near else "")
            )
        return path.read_text(encoding="utf-8", errors="replace"), path

    def build(self, one, options, template_text):
        incar = Incar.from_str(template_text)
        # 变胞下基组固定会带来 Pulay 应力，截断能须由收敛测试定，不替用户猜一个倍数
        # （用户定：曾经的 1.3 × ENMAX 不够）
        if options.task == "opt" and options.vary_cell and "ENCUT" not in incar:
            raise JobError("--vary-cell needs ENCUT in the template (Pulay stress): find it with "
                           "`ferro-py conv -s vasp --scan encut=...` and write it there")
        # PREC sets the coarse-grid factor of the fine grid; abbreviations (Acc, A) are not
        # guessed and only Normal / Accurate are taken (user decisions)
        prec = incar.get("PREC", "Normal").strip()
        if prec.upper() not in (p.upper() for p in PREC_VALUES):
            raise JobError(f"PREC = {prec}: write PREC in full as {' or '.join(PREC_VALUES)}")
        elements = one.elements()   # sort_atoms 之后即 POSCAR 的元素顺序
        symbols = one.symbols(0)
        counts = [symbols.count(e) for e in elements]
        notes = []

        potcars = [self.potcar(e) for e in elements]
        zval = [_potcar_field(text, "ZVAL") for text, _ in potcars]
        enmax = max(_potcar_field(text, "ENMAX") for text, _ in potcars)
        notes.append("POTCAR: " + " ".join(str(p.parent.name) for _, p in potcars))

        # meta-GGA（r2SCAN 等）用 ALGO = All + ISEARCH = 1：Normal（Davidson）对 meta-GGA
        # 收敛差，用户在 r2SCAN 标注里实测过这一组
        metagga = incar.get("METAGGA", "").strip().upper()
        if metagga and metagga not in ("NONE", "--"):
            if incar.get("ALGO", "").strip().upper() != "ALL" or incar.get("ISEARCH", "").strip() != "1":
                incar["ALGO"] = "All"
                incar["ISEARCH"] = 1
                notes.append(f"METAGGA = {metagga}: ALGO = All, ISEARCH = 1")
        # meta-GGA 的 egg-box 噪声对细网格敏感：显式写 NGXF/NGYF/NGZF（模板已写且未给
        # --fine-grid 就沿用模板）。GGA 给了 --fine-grid 也照写，只提示多半用不着
        if options.fine_grid is not None or (metagga not in ("", "NONE", "--") and "NGXF" not in incar):
            factor = options.fine_grid if options.fine_grid is not None else 2
            encut = float(incar["ENCUT"]) if "ENCUT" in incar else enmax
            coarse, fine = fine_grid(one.cell(0), encut, factor, prec)
            for key, n in zip(("NGXF", "NGYF", "NGZF"), fine):
                incar[key] = n
            notes.append(f"fine grid {factor:g} x NG {coarse} (PREC {prec}) -> NGXF/NGYF/NGZF {fine}"
                         + (" (Zn present: a factor of 3 may be needed)" if "Zn" in elements and factor < 3 else ""))
            if metagga in ("", "NONE", "--"):
                notes.append("no METAGGA: a GGA rarely needs an explicit fine grid")

        n_valence = sum(c * z for c, z in zip(counts, zval))
        mult = multiplicity(n_valence, options.charge, options.multiplicity)
        # A template's NELECT / NUPDOWN belong to the system it was tested on: left in,
        # they override the count and multiplicity reported below (review R8-1).
        # Commented out rather than deleted (user's choice), so the INCAR still shows them
        if options.charge:
            incar["NELECT"] = f"{n_valence - options.charge:g}"
        elif "NELECT" in incar:
            old = incar.comment_out("NELECT")
            notes.append(f"template NELECT = {old} commented out: neutral structure, VASP counts the electrons")
        if mult > 1:
            incar["ISPIN"] = 2
            incar["NUPDOWN"] = mult - 1
        elif options.multiplicity is not None and incar.get("ISPIN", "1").strip() == "2":
            incar["NUPDOWN"] = 0
        elif "NUPDOWN" in incar:
            old = incar.comment_out("NUPDOWN")
            notes.append(f"template NUPDOWN = {old} commented out: multiplicity 1, the moment relaxes "
                         "in the SCF (give --multiplicity to fix it)")
        notes.append(f"{n_valence - options.charge:g} electrons, multiplicity {mult}")

        if options.kspacing is not None:
            incar["KSPACING"] = options.kspacing
        elif "KSPACING" not in incar:
            incar["KSPACING"] = DEFAULT_KSPACING
        try:
            spacing = float(incar["KSPACING"])
        except ValueError:
            raise JobError(f"KSPACING = {incar['KSPACING']} is not a number") from None
        mesh = kmesh(one.cell(0), spacing)
        notes.append(f"PREC {prec}, KSPACING {spacing:g} -> mesh {' x '.join(map(str, mesh))}"
                     + (" (Gamma only: vasp_gam can run it)" if mesh == [1, 1, 1] else ""))

        if options.task == "opt" and options.vary_cell:
            incar["ISIF"] = 3
        if options.task == "aimd":
            if options.ensemble not in (None, "nvt"):
                raise JobError(f"VASP jobs are NVT only; --ensemble {options.ensemble} is a CP2K option")
            if options.temperature is not None:
                incar["TEBEG"] = f"{options.temperature:g}"
                if "TEEND" in incar:
                    incar["TEEND"] = f"{options.temperature:g}"
            if options.steps is not None:
                incar["NSW"] = options.steps
            if options.timestep is not None:
                incar["POTIM"] = f"{options.timestep:g}"

        if not any(k in incar for k in PARALLEL_KEYS):
            notes.append("no NCORE / NPAR / KPAR in the INCAR: add your machine's parallel settings")
        poscar = ferro.poscar_string(one, 0, " ".join(elements))
        return InputSet({"INCAR": incar.to_str(), "POSCAR": poscar,
                         "POTCAR": "".join(text for text, _ in potcars)},
                        shared={"POTCAR"}, notes=notes)
