"""三层结构（借鉴 pymatgen ``io.core``）：参数块类（``Incar``、``QeInput``、``Section``）
→ ``InputSet``（一个任务目录里的全部文件）→ ``Generator``（模板 + 结构 + 选项 → InputSet）。

``Generator.generate`` 做三个程序共有的检查与准备，程序相关的部分在子类的 ``build``。
加一个程序（如 Gaussian）就是加一个子类与它的模板；加一个功能（如固定原子）就是
``JobOptions`` 加一项、各子类的 ``build`` 接上。
"""

import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

import ferro

TASKS = ("sp", "opt", "aimd")
ENSEMBLES = ("nvt", "npt_i", "npt_f")
BUILTIN_TEMPLATES = Path(__file__).parent / "templates"
# k-point spacing [1/Å, 2π included] when neither --kspacing nor the template gives a mesh
DEFAULT_KSPACING = 0.5


class JobError(ValueError):
    pass


@dataclass
class JobOptions:
    """程序在模板之上改的量。``None`` 表示不改，沿用模板。"""

    task: str
    charge: int = 0
    multiplicity: int | None = None   # None：按电子数奇偶取 1 或 2
    kspacing: float | None = None     # Å⁻¹，含 2π（VASP KSPACING 的定义）；None 用程序默认
    temperature: float | None = None  # 以下四项只用于 aimd
    steps: int | None = None
    timestep: float | None = None     # fs
    ensemble: str | None = None
    vary_cell: bool = False           # 只用于 opt
    fine_grid: float | None = None    # VASP：NGXF = fine_grid × NGX；CP2K：USE_FINER_GRID

    def check(self):
        if self.task not in TASKS:
            raise JobError(f"unknown task '{self.task}'; one of {', '.join(TASKS)}")
        md = {"temperature": self.temperature, "steps": self.steps,
              "timestep": self.timestep, "ensemble": self.ensemble}
        given = [k for k, v in md.items() if v is not None]
        if given and self.task != "aimd":
            raise JobError(f"{', '.join(given)} only apply to an aimd task, not {self.task}")
        if self.vary_cell and self.task != "opt":
            raise JobError(f"vary_cell only applies to an opt task, not {self.task}")
        if self.ensemble is not None and self.ensemble not in ENSEMBLES:
            raise JobError(f"unknown ensemble '{self.ensemble}'; one of {', '.join(ENSEMBLES)}")
        for name in ("temperature", "timestep", "kspacing"):
            value = getattr(self, name)
            if value is not None and not (math.isfinite(value) and value > 0):
                raise JobError(f"{name} must be a positive number, got {value}")
        if self.fine_grid is not None and not (math.isfinite(self.fine_grid) and self.fine_grid >= 1):
            raise JobError(f"fine_grid must be at least 1, got {self.fine_grid}")
        if self.steps is not None and self.steps < 1:
            raise JobError(f"steps must be at least 1, got {self.steps}")
        if self.multiplicity is not None and self.multiplicity < 1:
            raise JobError(f"multiplicity must be at least 1, got {self.multiplicity}")


class InputSet:
    """一个任务目录：文件名 → 文本。``shared`` 里的文件批量时可在上层目录共用一份；
    ``notes`` 是生成时要告诉用户的话（取了哪个多重度、哪些量沿用模板）。"""

    def __init__(self, files, shared=(), notes=()):
        self.files = dict(files)
        self.shared = set(shared)
        self.notes = list(notes)

    def write(self, directory, links=None):
        """写进 ``directory``（须已存在）。``links`` 里的文件名写成指向给定目标的符号链接。"""
        directory = Path(directory)
        links = links or {}
        for name, text in self.files.items():
            path = directory / name
            if path.is_symlink() or path.exists():
                path.unlink()
            if name in links:
                path.symlink_to(links[name])
            else:
                path.write_text(text, encoding="utf-8")


def write_batch(jobs, out_dir):
    """``[(子目录名, InputSet)]`` 写进 ``out_dir``（须已存在），返回各子目录的路径。

    多个任务时 ``shared`` 的文件（POTCAR）在 ``out_dir`` 只放一份，子目录里是相对路径
    的符号链接（``../POTCAR``）：整个目录挪走或 ``rsync -a`` 过去链接不断。各任务的共用
    文件不一致（组成不同）就报错，让用户按组成分开生成，不在一个目录里放几份。
    """
    out_dir = Path(out_dir)
    names = [name for name, _ in jobs]
    duplicated = sorted({n for n in names if names.count(n) > 1})
    if duplicated:
        raise JobError(f"two inputs would both write {out_dir / duplicated[0]}; rename one of them")
    shared = {}
    if len(jobs) > 1:
        for name, s in jobs:
            for f in s.shared:
                first = shared.setdefault(f, (name, s.files[f]))
                if first[1] != s.files[f]:
                    raise JobError(
                        f"{name} needs a different {f} from {first[0]} (another composition or "
                        f"element order); generate jobs of each composition into their own -o"
                    )
    # 先前一批留下的共用文件不同（另一组成）就拒绝：覆盖它会让那一批的 ../POTCAR 静默改指
    for f, (_, text) in shared.items():
        old = out_dir / f
        if old.is_file() and old.read_text(encoding="utf-8") != text:
            raise JobError(f"{old} already holds the {f} of another composition; "
                           f"generate jobs of each composition into their own -o")
    for f, (_, text) in shared.items():
        (out_dir / f).write_text(text, encoding="utf-8")
    paths = []
    for name, s in jobs:
        d = out_dir / name
        d.mkdir(parents=True, exist_ok=True)
        s.write(d, links={f: Path(os.path.relpath(out_dir / f, d)) for f in shared})
        paths.append(d)
    return paths


def multiplicity(n_valence, charge, requested):
    """自旋多重度：未给时按电子数奇偶取 1 或 2；给了就校验奇偶。

    电子数用所选赝势的价电子数，不用原子序数 —— f 电子冻在核里的镧系赝势（VASP
    ``Gd_3``）核内电子数是奇数，按原子序数会判反。
    """
    n = n_valence - charge
    if abs(n - round(n)) > 1e-6:
        raise JobError(f"the pseudopotentials give {n_valence:g} valence electrons, not a whole number")
    n = round(n)
    if n < 0:
        raise JobError(f"charge {charge} removes more than the {n_valence:g} valence electrons")
    if requested is None:
        return 1 if n % 2 == 0 else 2
    if (requested - 1) % 2 != n % 2:
        raise JobError(
            f"multiplicity {requested} is impossible with {n} electrons "
            f"({'even' if n % 2 == 0 else 'odd'} count needs an {'odd' if n % 2 == 0 else 'even'} multiplicity)"
        )
    return requested


def kmesh(cell, kspacing):
    """Γ 中心网格：N_i = max(1, ⌈|b_i|·2π / kspacing⌉)，b_i 为倒格矢（VASP KSPACING 的定义）。

    |b_i| = 1 / d_i, the interplanar spacing; the geometry is Rust's (`Cell::interplanar_spacings`).
    """
    return [max(1, math.ceil(2 * math.pi / (d * kspacing) - 1e-9)) for d in ferro.interplanar_spacings(cell)]


class Generator:
    """子类给出 ``program``、``template_name``（模板文件名，也是写出的主文件名）与 ``build``。"""

    program = None
    template_name = None

    def __init__(self, config):
        self.config = config

    def template(self, task, path=None):
        """模板文本：显式 ``path`` > 用户模板目录 ``<dir>/<程序>/<任务>/<文件名>`` > 内置。"""
        if path is None:
            user = self.config["templates"].get("dir")
            candidate = user / self.program / task / self.template_name if user else None
            if candidate is not None and candidate.exists():
                path = candidate
            else:
                return self.builtin_template(task)
        path = Path(path)
        if not path.is_file():
            raise JobError(f"template {path} does not exist")
        return path.read_text(encoding="utf-8")

    @classmethod
    def builtin_template(cls, task):
        """随包的通用模板（``ferro-py job template`` 导出的就是它）。"""
        path = BUILTIN_TEMPLATES / cls.program / task / cls.template_name
        if not path.is_file():
            raise JobError(f"no built-in {cls.program} template for task '{task}'")
        return path.read_text(encoding="utf-8")

    def generate(self, traj, frame, options, template=None, template_text=None):
        """``traj`` 的第 ``frame`` 帧 → InputSet。原子按 (Z, 符号) 排序，同 ``ferro convert``。"""
        options.check()
        one = ferro.sort_atoms(traj.subset([frame]))
        if one.cell(0) is None or not all(one.pbc(0)):
            raise JobError("only periodic systems (a cell, periodic along a, b and c) are supported")
        # `X<type>` / `X` is what the Rust readers name an atom whose element they do not know
        # (ferro-io readers/lammps_dump.rs and lammps_data.rs `X<type>`, cube.rs `X<Z>`, a label
        # with no leading letter in readers/util.rs `element_from_label` and cif.rs `X`)
        placeholders = sorted(e for e in one.elements() if re.fullmatch(r"X\d*", e))
        if placeholders:
            raise JobError(
                f"atoms {', '.join(placeholders)} have no element: a dump without an element "
                "column needs --type-map (e.g. 1=Zr,2=O); a label needs a leading element symbol"
            )
        if template_text is None:
            template_text = self.template(options.task, template)
        return self.build(one, options, template_text)

    def build(self, one, options, template_text):
        raise NotImplementedError
