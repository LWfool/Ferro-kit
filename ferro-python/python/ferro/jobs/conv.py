"""收敛测试：一个参数扫一组值，每个值生成同一批结构的单点；读回输出与最高档比较。

扫哪个参数、怎么改模板，各程序一张表（``SCANS``）；目录布局与参数值记在
``ferro-conv.json`` 里，``analyze`` 只读它，不从目录名反解析。
"""

import json
import math
from pathlib import Path

import ferro
from ferro.jobs.base import JobError, JobOptions
from ferro.jobs.cp2k import Section
from ferro.jobs.vasp import Incar

MANIFEST = "ferro-conv.json"
EV_A3_TO_KBAR = 1602.1766208   # 1 eV/Å³ = 160.21766208 GPa
SWITCH = ("off", "on")


def _vasp_set(key):
    def apply(text, value, options):
        incar = Incar.from_str(text)
        incar[key] = f"{value:g}"
        return incar.to_str(), options
    return apply


def _vasp_fine_grid(text, value, options):
    options.fine_grid = value
    return text, options


def _cp2k_set(path, key, switch=False):
    def apply(text, value, options):
        root = Section.from_str(text)
        if switch:
            section = root.section(path)
            if value == "on":
                root.section(path, create=True).set(key, "T")
            elif section is not None:
                section.remove(key)
        else:
            section = root.section(path)
            if section is None:
                raise JobError(f"the template has no &{path.replace('/', ' &')} for {key}")
            section.set(key, f"{value:g}")
        return root.to_str(), options
    return apply


def _cp2k_needs_gapw(apply):
    def check(text, value, options):
        qs = Section.from_str(text).section("FORCE_EVAL/DFT/QS")
        method = (qs.get("METHOD", "GPW") if qs else "GPW").split()[0].upper()
        if method != "GAPW":
            raise JobError("GAPW_ACCURATE_XCINT only acts with &QS METHOD GAPW (all-electron)")
        return apply(text, value, options)
    return check


# 程序 → 参数 → (改模板的函数, 取值类型, 单位)。取值类型 float 为数值档，SWITCH 为关 / 开
SCANS = {
    "vasp": {
        "encut": (_vasp_set("ENCUT"), float, "eV"),
        "fine-grid": (_vasp_fine_grid, float, "x NG"),
    },
    "cp2k": {
        "cutoff": (_cp2k_set("FORCE_EVAL/DFT/MGRID", "CUTOFF"), float, "Ry"),
        "rel-cutoff": (_cp2k_set("FORCE_EVAL/DFT/MGRID", "REL_CUTOFF"), float, "Ry"),
        "finer-grid": (_cp2k_set("FORCE_EVAL/DFT/XC/XC_GRID", "USE_FINER_GRID", switch=True), SWITCH, ""),
        "accurate-xcint": (_cp2k_needs_gapw(_cp2k_set("FORCE_EVAL/DFT/QS", "GAPW_ACCURATE_XCINT", switch=True)),
                           SWITCH, ""),
    },
}


def parse_scan(software, text):
    """``encut=400:900:50`` / ``cutoff=300,400,600`` / ``finer-grid`` → (参数, [值])。"""
    if software not in SCANS:
        raise JobError(f"convergence tests are offered for {', '.join(SCANS)}, not {software}")
    name, sep, values = text.partition("=")
    if name not in SCANS[software]:
        raise JobError(f"{software} scans {', '.join(SCANS[software])}, not '{name}'")
    kind = SCANS[software][name][1]
    if kind is SWITCH:
        if sep:
            raise JobError(f"{name} is a switch: give --scan {name} (off and on), no values")
        return name, list(SWITCH)
    if not values:
        raise JobError(f"--scan {name}=VALUES: e.g. {name}=400:900:50 or {name}=400,500,700")
    try:
        if ":" in values:
            start, stop, step = (float(x) for x in values.split(":"))
            if step <= 0 or stop < start:
                raise ValueError
            out = [start + i * step for i in range(int(math.floor((stop - start) / step + 1e-9)) + 1)]
        else:
            out = [float(x) for x in values.split(",")]
    except ValueError:
        raise JobError(f"--scan {text}: values are START:STOP:STEP or a comma list of numbers") from None
    if len(out) < 2 or len(set(out)) != len(out):
        raise JobError(f"--scan {text}: needs at least two distinct values")
    return name, sorted(out)


def label(value):
    return value if isinstance(value, str) else f"{value:g}"


def build(generator, software, parameter, values, structures, options, template=None):
    """结构 [(名, Trajectory, 帧号)] × 参数值 → [(子目录, InputSet)] 与清单。"""
    if options.task != "sp":
        raise JobError("convergence tests are single points")
    apply = SCANS[software][parameter][0]
    base = generator.template("sp", template)
    jobs, manifest = [], {"software": software, "parameter": parameter,
                          "unit": SCANS[software][parameter][2], "values": values, "jobs": []}
    for value in values:
        text, opts = apply(base, value, JobOptions(**vars(options)))
        for name, traj, frame in structures:
            directory = f"{parameter}_{label(value)}/{name}"
            jobs.append((directory, generator.generate(traj, frame, opts, template_text=text)))
            manifest["jobs"].append({"dir": directory, "value": value, "structure": name, "frame": frame})
    return jobs, manifest


# ── 分析 ─────────────────────────────────────────────────────────────────────

OUTPUT_NAMES = {"vasp": ("OUTCAR", "vasprun.xml"), "cp2k": ("cp2k.out", "*.out", "*.log")}
# 提示用的阈值（2026-10-16 查惯例：能量 1 meV/atom 为常用判据，SSSP 结合能 2 meV/atom；
# 压力 1 kbar 取自 SSSP 对 vc-relax 的提醒；力无权威标准，取常见势函数力误差的约 1/10）
THRESHOLDS = {"de": 1.0, "dp": 1.0, "df_rms": 5.0, "df_max": 10.0}


def read_result(directory, software):
    """任务目录 → (能量 eV, 力 [(x,y,z)], 压力 kbar, 状态)；读不到时数值为 None、状态说明原因。"""
    directory = Path(directory)
    candidates = []
    for pattern in OUTPUT_NAMES[software]:
        candidates += sorted(directory.glob(pattern))
    if not candidates:
        return None, None, None, "no output"
    last_error = None
    for path in candidates:
        try:
            traj, stats = ferro.read_output(str(path))
        except RuntimeError as e:
            last_error = str(e)
            continue
        if len(traj) == 0:
            return None, None, None, "SCF not converged" if stats.get("scf_failed") else "no complete step"
        f = len(traj) - 1
        stress = traj.stress(f)
        pressure = None if stress is None else sum(stress[i][i] for i in range(3)) / 3 * EV_A3_TO_KBAR
        return traj.energy(f), traj.forces(f), pressure, "ok"
    return None, None, None, f"unreadable ({last_error})"


def compare(manifest, out_dir, ref=None, thresholds=None):
    """每帧每档一行：与参考档（默认最高档，开关为 on）的偏差。返回 (行, 参考值)。"""
    thresholds = {**THRESHOLDS, **(thresholds or {})}
    values = manifest["values"]
    if ref is None:
        ref = values[-1]
    elif ref not in values:
        raise JobError(f"--ref {label(ref)} is not one of the scanned values ({', '.join(map(label, values))})")
    results = {}
    for job in manifest["jobs"]:
        results[(job["structure"], job["value"])] = read_result(Path(out_dir) / job["dir"], manifest["software"])
    rows = []
    for job in manifest["jobs"]:
        name, value = job["structure"], job["value"]
        energy, forces, pressure, status = results[(name, value)]
        r_energy, r_forces, r_pressure, r_status = results[(name, ref)]
        row = {"value": label(value), "structure": name, "frame": str(job["frame"]), "status": status,
               "n_atoms": "", "energy_eV": math.nan, "dE_meV_atom": math.nan, "pressure_kbar": math.nan,
               "dP_kbar": math.nan, "dF_max_meV_A": math.nan, "dF_rms_meV_A": math.nan, "within": ""}
        if status == "ok":
            row["energy_eV"] = energy
            row["n_atoms"] = str(len(forces)) if forces else ""
            if pressure is not None:
                row["pressure_kbar"] = pressure
        if status == "ok" and r_status == "ok":
            n = len(forces)
            if len(r_forces) != n:
                raise JobError(f"{name}: {n} atoms at {label(value)}, {len(r_forces)} at the reference")
            row["dE_meV_atom"] = (energy - r_energy) / n * 1000
            if pressure is not None and r_pressure is not None:
                row["dP_kbar"] = pressure - r_pressure
            d = [math.dist(a, b) * 1000 for a, b in zip(forces, r_forces)]
            row["dF_max_meV_A"] = max(d)
            row["dF_rms_meV_A"] = math.sqrt(sum(x * x for x in d) / n)
            checks = [abs(row["dE_meV_atom"]) <= thresholds["de"], row["dF_rms_meV_A"] <= thresholds["df_rms"],
                      row["dF_max_meV_A"] <= thresholds["df_max"]]
            if not math.isnan(row["dP_kbar"]):
                checks.append(abs(row["dP_kbar"]) <= thresholds["dp"])
            row["within"] = "yes" if all(checks) else "no"
        elif status == "ok":
            row["status"] = f"reference {r_status}"
        rows.append(row)
    return rows, ref


def summary(rows, values):
    """每档取各帧的最大偏差（绝对值），并找出第一档「全部帧都在阈值内」的值。"""
    out = []
    for value in map(label, values):
        mine = [r for r in rows if r["value"] == value]
        def worst(key):
            xs = [abs(r[key]) for r in mine if not math.isnan(r[key])]
            return max(xs) if xs else math.nan
        ok = [r for r in mine if r["status"] == "ok"]
        out.append({"value": value, "frames": f"{len(ok)}/{len(mine)}", "dE": worst("dE_meV_atom"),
                    "dP": worst("dP_kbar"), "dF_max": worst("dF_max_meV_A"), "dF_rms": worst("dF_rms_meV_A"),
                    "within": all(r["within"] == "yes" for r in mine) if mine else False})
    return out
