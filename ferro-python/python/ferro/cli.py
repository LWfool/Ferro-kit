"""``ferro-py``：ferro 的 Python 命令行。子命令 ``config``、``job``，以后加 ``plot``。

帮助页与 Rust 的 ``ferro`` 同为五段（用途 / Parameters / Output / Examples /
Full documentation），手写在各子命令旁；``tests/test_cli.py`` 核对帮助页与
argparse 的参数一致（同 ``ferro-cli`` 的 ``help_sync``）。
"""

import argparse
import glob
import math
import os
import sys
from pathlib import Path

import ferro
from ferro import config

OVERVIEW = """\
ferro-py — the Python side of ferro

Usage:
  ferro-py <COMMAND> [OPTIONS]

Commands:
  job           Generate VASP / QE / CP2K inputs (sp, opt, aimd) from structures
  job template  Write a built-in template to edit and pass back with --template
  conv          Generate a convergence test (VASP fine grid / ENCUT, CP2K cutoffs)
  conv analyze  Compare its single points: energy, pressure, forces
  config init   Write a commented configuration template

Run `ferro-py <COMMAND> -h` for the parameters of one command.

Full documentation:  ferro doc python"""

CONFIG_INIT_HELP = """\
ferro-py config init — Write a commented configuration template

  Writes the global configuration that ferro-py reads: where the POTCARs, the QE
  pseudopotentials, CP2K's data files and your own templates are. Every entry is
  commented out; uncomment the ones you use. An existing file is never overwritten
  without --force.

Parameters:
  -o, --output PATH   File to write   [$FERRO_CONFIG, else ~/.config/ferro/config.toml]
      --force         Overwrite an existing file
      --mkdir         Create -o's directory without asking
  -h, --help          This page

Output:
  The configuration file. The default location's directory is created if missing;
  a missing directory of an explicit -o is asked about (--mkdir when not interactive).

Examples:
  ferro-py config init
  ferro-py config init -o ~/hpc.toml --mkdir     # then: export FERRO_CONFIG=~/hpc.toml

Full documentation:  ferro doc python config"""


JOB_HELP = """\
ferro-py job — Generate calculation inputs from structures

  Fills a template (a native input file of the program) with each structure and a
  few quantities: charge, multiplicity, k-point mesh, the element list of the
  pseudopotentials / basis sets, and for aimd the temperature, steps, time step and
  ensemble. Everything else stays as the template has it. Periodic systems only.

Parameters:
  -s, --software S      vasp | qe | cp2k  (required)
  -t, --task T          sp | opt | aimd  (required)
  -i, --input FILE      Structure file(s), any format ferro reads; quoted globs allowed
  -o, --output DIR      Directory receiving one sub-directory per job  (required)
      --template FILE   Template instead of the configured / built-in one
      --charge N        Total charge  [0]
      --multiplicity M  2S+1, checked against the electron count  [1 if even, 2 if odd]
      --kspacing K      k-point spacing in 1/Å, 2π included (as VASP KSPACING); a
                        mesh written in the template is kept  [per program, see manual]
      --temperature T   aimd: temperature in K  [the template's]
      --steps N         aimd: number of MD steps  [the template's]
      --timestep FS     aimd: time step in fs  [the template's]
      --ensemble E      aimd: nvt | npt_i | npt_f (npt_* CP2K only)  [the template's]
      --vary-cell       opt: relax the cell too; vasp needs ENCUT in the template  [cell fixed]
      --fine-grid F     vasp: NGXF/NGYF/NGZF = F x NGX/NGY/NGZ  [2 with METAGGA, else
                        unset]; cp2k: USE_FINER_GRID on; qe: ignored (no separate grid)
      --start N         First frame to take  (0-based, inclusive)  [0]
      --end N           Last frame to take  (0-based, INCLUSIVE)  [last]
      --stride N        Take every Nth frame within [start, end]  [1]
      --number N        Take this many frames, spread evenly, both ends kept
      --atom-style S    LAMMPS data input: atomic | charge | full  (required)
      --units U         Dump vx/fx units: real | metal  (required if present)
      --type-map M      Dump without element column: 1=Zr,2=O,...
      --mkdir           Create -o without asking
  -h, --help            This page

Output:
  DIR/<input name>/ for a single-structure file, DIR/<input name>_<frame>/ (frame
  0-based, 4 digits) for each frame taken from a trajectory. File names are fixed:
  INCAR KPOINTS POSCAR POTCAR (vasp), pw.in (qe), cp2k.inp (cp2k). With several
  jobs the POTCAR is written once in DIR and linked as ../POTCAR; jobs that need
  different POTCARs are an error. All frames by default.

Examples:
  ferro-py job -s vasp -t sp -i POSCAR -o sp/
  ferro-py job -s cp2k -t aimd -i traj.lammpstrj --number 20 --temperature 1500 -o md/
  ferro-py job -s qe -t opt -i 'cifs/*.cif' --vary-cell --template my_relax.in -o relax/

Full documentation:  ferro doc python job"""


JOB_TEMPLATE_HELP = """\
ferro-py job template — Write the built-in template of one job type, to edit

  The template is a native input file of the program. Its header lists what
  ferro-py rewrites in it and what it deliberately leaves out (parallel settings).
  Edit it, then pass it back with `ferro-py job --template FILE`.

Parameters:
  -s, --software S      vasp | qe | cp2k  (required)
  -t, --task T          sp | opt | aimd  (required)
  -o, --output DIR      Directory to write the template into  (required)
      --mkdir           Create -o without asking
  -h, --help            This page

Output:
  DIR/INCAR (vasp), DIR/pw.in (qe) or DIR/cp2k.inp (cp2k). An existing file is
  never overwritten.

Examples:
  ferro-py job template -s vasp -t aimd -o my_aimd/
  ferro-py job -s vasp -t aimd -i POSCAR --template my_aimd/INCAR -o md/

Full documentation:  ferro doc python job"""


class UsageError(Exception):
    pass


def parse(parser, help_text, argv):
    """``-h`` 打印手写的五段帮助页；argparse 自带的帮助页关掉，免得两份并存。"""
    if "-h" in argv or "--help" in argv:
        print(help_text)
        raise SystemExit(0)
    return parser.parse_args(argv)


def new_parser(prog):
    parser = argparse.ArgumentParser(prog=prog, add_help=False)
    parser.add_argument("-h", "--help", action="store_true")
    return parser


def ensure_dir(directory, mkdir):
    """同 ``ferro-cli/src/outpath.rs`` 的 ``ensure_dir``：缺目录先问，非交互须 ``--mkdir``。

    敲错的路径若被静默建出来，看上去就像一次成功的运行，所以要问；默认答「否」。
    """
    directory = Path(directory)
    if directory.exists():
        if not directory.is_dir():
            raise UsageError(f"-o '{directory}' exists and is not a directory")
        return
    if not mkdir:
        if not sys.stdin.isatty():
            raise UsageError(
                f"output directory '{directory}' does not exist, and there is no terminal "
                "to ask. Pass --mkdir to create it."
            )
        print(f"Directory '{directory}' does not exist. Create it? [y/N] ", end="", file=sys.stderr)
        sys.stderr.flush()
        if sys.stdin.readline().strip() not in ("y", "Y", "yes", "Yes"):
            raise UsageError(f"no output directory: '{directory}' was not created")
    directory.mkdir(parents=True)
    print(f"Created: {directory}")


def config_parser():
    p = new_parser("ferro-py config init")
    p.add_argument("-o", "--output")
    p.add_argument("--force", action="store_true")
    p.add_argument("--mkdir", action="store_true")
    return p


def run_config(argv):
    if not argv or argv[0] != "init":
        if argv and argv[0] not in ("-h", "--help"):
            raise UsageError(f"unknown config action '{argv[0]}'; the only one is `init`")
        print(CONFIG_INIT_HELP)
        return 0
    args = parse(config_parser(), CONFIG_INIT_HELP, argv[1:])
    if args.output is None:
        # 默认位置由程序定、不是用户敲的，没有敲错的可能，目录直接建
        path = config.init(force=args.force)
    else:
        ensure_dir(Path(args.output).parent, args.mkdir)
        path = config.init(args.output, force=args.force)
    print(f"Wrote: {path}")
    return 0


def job_parser():
    p = new_parser("ferro-py job")
    p.add_argument("-s", "--software", required=True)
    p.add_argument("-t", "--task", required=True)
    p.add_argument("-i", "--input", nargs="+", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--template")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--multiplicity", type=int)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--temperature", type=float)
    p.add_argument("--steps", type=int)
    p.add_argument("--timestep", type=float)
    p.add_argument("--ensemble")
    p.add_argument("--vary-cell", action="store_true")
    p.add_argument("--fine-grid", type=float)
    p.add_argument("--start", type=int)
    p.add_argument("--end", type=int)
    p.add_argument("--stride", type=int)
    p.add_argument("--number", type=int)
    p.add_argument("--atom-style")
    p.add_argument("--units")
    p.add_argument("--type-map")
    p.add_argument("--mkdir", action="store_true")
    return p


def parse_type_map(text):
    """``1=Zr,2=P`` → ``{1: "Zr", 2: "P"}``，语法同 CLI ``ferro`` 的 ``--type-map``。"""
    mapping = {}
    for item in text.split(","):
        t, sep, element = item.partition("=")
        if not sep:
            raise UsageError(f"--type-map: '{item}' is not TYPE=ELEMENT (e.g. 1=Zr,2=O)")
        try:
            t = int(t)
        except ValueError:
            raise UsageError(f"--type-map: type '{t}' is not a number") from None
        if t in mapping:
            raise UsageError(f"--type-map: type {t} is given twice")
        mapping[t] = element.strip()
    return mapping


def expand_inputs(patterns):
    """带通配符的 -i 自己展开（加了引号，shell 没展开）；不带的原样，不存在留给读取时报错。"""
    files = []
    for pattern in patterns:
        if any(c in pattern for c in "*?["):
            found = sorted(glob.glob(pattern))
            if not found:
                raise UsageError(f"-i '{pattern}' matches no file")
            files += found
        else:
            files.append(pattern)
    return files


def job_template_parser():
    p = new_parser("ferro-py job template")
    p.add_argument("-s", "--software", required=True)
    p.add_argument("-t", "--task", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--mkdir", action="store_true")
    return p


def find_generator(software):
    from ferro import jobs

    generator = jobs.GENERATORS.get(software)
    if generator is None:
        raise UsageError(f"unknown software '{software}'; one of {', '.join(jobs.GENERATORS)}")
    return generator


def run_job_template(argv):
    from ferro.jobs.base import JobOptions

    args = parse(job_template_parser(), JOB_TEMPLATE_HELP, argv)
    generator = find_generator(args.software)
    JobOptions(task=args.task).check()
    # 导出内置模板，不取用户模板目录里的：导出的是「通用设置」，用户在它上面改
    text = generator.builtin_template(args.task)
    ensure_dir(args.output, args.mkdir)
    path = Path(args.output) / generator.template_name
    if path.exists():
        raise UsageError(f"{path} already exists; it is never overwritten")
    path.write_text(text, encoding="utf-8")
    print(f"Wrote: {path}")
    return 0


def run_job(argv):
    from ferro.jobs.base import JobOptions, write_batch

    if not argv:
        print(JOB_HELP)
        return 0
    if argv[0] == "template":
        return run_job_template(argv[1:])
    args = parse(job_parser(), JOB_HELP, argv)
    # 参数错误在读第一个文件之前报
    generator = find_generator(args.software)
    options = JobOptions(
        task=args.task, charge=args.charge, multiplicity=args.multiplicity,
        kspacing=args.kspacing, temperature=args.temperature, steps=args.steps,
        timestep=args.timestep, ensemble=args.ensemble, vary_cell=args.vary_cell,
        fine_grid=args.fine_grid,
    )
    options.check()
    type_map = parse_type_map(args.type_map) if args.type_map else None
    files = expand_inputs(args.input)
    generator = generator(config.load())
    ensure_dir(args.output, args.mkdir)

    built = []
    for name, traj, i, f in read_structures(files, args, type_map):
        try:
            built.append((name, generator.generate(traj, i, options, args.template)))
        except (ValueError, RuntimeError) as e:
            raise UsageError(f"{f} frame {i}: {e}") from None
    # 先全部生成、再写：一个失败就什么都不写，不留半截的批次
    report_batch(built, write_batch(built, args.output))
    return 0


def read_structures(files, args, type_map):
    """-i 的文件 × 选帧 → [(任务名, Trajectory, 帧号, 文件)]。单结构文件的任务名不带帧号。"""
    out = []
    for f in files:
        traj = ferro.read(f, units=args.units, atom_style=args.atom_style, type_map=type_map)
        for i in traj.frame_indices(args.start, args.end, args.stride, args.number):
            out.append((Path(f).stem if len(traj) == 1 else f"{Path(f).stem}_{i:04d}", traj, i, f))
    return out


def report_batch(built, paths):
    for (name, s), path in zip(built, paths):
        listing = " ".join(f"{f}->{os.readlink(path / f)}" if (path / f).is_symlink() else f for f in s.files)
        print(f"Wrote: {path}/  {listing}")
        for note in s.notes:
            print(f"  {note}")


CONV_HELP = """\
ferro-py conv — Generate a convergence test: one parameter, several values

  Writes the same single points for every value of one parameter, the rest taken
  from the sp template; `ferro-py conv analyze DIR` then compares them. For a
  meta-GGA, test VASP's fine grid first and ENCUT after.

Parameters:
  -s, --software S      vasp | cp2k  (required)
      --scan P[=VALUES] The parameter and its values, START:STOP:STEP or a,b,c  (required)
                        vasp: encut [eV], fine-grid [x NG]
                        cp2k: cutoff [Ry], rel-cutoff [Ry], finer-grid, accurate-xcint
                        (the last two are switches: off and on, no values)
  -i, --input FILE      Structure file(s), any format ferro reads; quoted globs allowed
  -o, --output DIR      Directory receiving the test  (required)
      --template FILE   sp template instead of the configured / built-in one
      --charge N        Total charge  [0]
      --multiplicity M  2S+1, checked against the electron count  [1 if even, 2 if odd]
      --kspacing K      k-point spacing in 1/Å, 2π included  [per program]
      --fine-grid F     vasp: NGXF = F x NGX at every value  [2 with METAGGA]
      --start N         First frame to take  (0-based, inclusive)  [0]
      --end N           Last frame to take  (0-based, INCLUSIVE)  [last]
      --stride N        Take every Nth frame within [start, end]  [1]
      --number N        Take this many frames, spread evenly, both ends kept
      --atom-style S    LAMMPS data input: atomic | charge | full  (required)
      --units U         Dump vx/fx units: real | metal  (required if present)
      --type-map M      Dump without element column: 1=Zr,2=O,...
      --mkdir           Create -o without asking
  -h, --help            This page

Output:
  DIR/<parameter>_<value>/<structure>/ per value and structure, laid out as
  `ferro-py job` does (one POTCAR in DIR, linked), and DIR/ferro-conv.json listing
  them for `conv analyze`.

Examples:
  ferro-py conv -s vasp --scan fine-grid=2,2.5,3,4 -i POSCAR --template r2scan/INCAR -o ng/
  ferro-py conv -s vasp --scan encut=500:900:100 -i traj.lammpstrj --number 3 -o encut/
  ferro-py conv -s cp2k --scan finer-grid -i POSCAR -o fg/

Full documentation:  ferro doc python conv"""

CONV_ANALYZE_HELP = """\
ferro-py conv analyze — Compare the single points of a convergence test

  Reads every job's output (VASP OUTCAR / vasprun.xml; CP2K cp2k.out, *.out or
  *.log) and compares each value with the reference, frame by frame: energy per
  atom, pressure, the largest per-atom force difference and the RMS over atoms of
  |F_i - F_i(ref)|. The thresholds only mark rows; they decide nothing.

Parameters:
  DIR                   Directory written by `ferro-py conv`  (required)
  -o, --output FILE     CSV to write  [DIR/conv_<parameter>.csv]
      --ref VALUE       Reference value  [the highest; on for a switch]
      --de X            Threshold |dE| in meV/atom  [1]
      --dp X            Threshold |dP| in kbar  [1]
      --df-rms X        Threshold force RMS in meV/Å  [5]
      --df-max X        Threshold largest force difference in meV/Å  [10]
  -h, --help            This page

Output:
  The CSV: one row per value and structure — value, structure, frame, status,
  n_atoms, energy_eV, dE_meV_atom, pressure_kbar, dP_kbar, dF_max_meV_A,
  dF_rms_meV_A, within (yes / no against the thresholds). The terminal shows the
  worst frame of each value and the first value within every threshold.

Examples:
  ferro-py conv analyze encut/
  ferro-py conv analyze ng/ --ref 4 --df-rms 2

Full documentation:  ferro doc python conv"""


def conv_parser():
    p = new_parser("ferro-py conv")
    p.add_argument("-s", "--software", required=True)
    p.add_argument("--scan", required=True)
    p.add_argument("-i", "--input", nargs="+", required=True)
    p.add_argument("-o", "--output", required=True)
    p.add_argument("--template")
    p.add_argument("--charge", type=int, default=0)
    p.add_argument("--multiplicity", type=int)
    p.add_argument("--kspacing", type=float)
    p.add_argument("--fine-grid", type=float)
    for name in ("--start", "--end", "--stride", "--number"):
        p.add_argument(name, type=int)
    p.add_argument("--atom-style")
    p.add_argument("--units")
    p.add_argument("--type-map")
    p.add_argument("--mkdir", action="store_true")
    return p


def conv_analyze_parser():
    p = new_parser("ferro-py conv analyze")
    p.add_argument("dir")
    p.add_argument("-o", "--output")
    p.add_argument("--ref")
    for name in ("--de", "--dp", "--df-rms", "--df-max"):
        p.add_argument(name, type=float)
    return p


def run_conv(argv):
    import json

    from ferro.jobs import conv
    from ferro.jobs.base import JobOptions, write_batch

    if not argv:
        print(CONV_HELP)
        return 0
    if argv[0] == "analyze":
        return run_conv_analyze(argv[1:])
    args = parse(conv_parser(), CONV_HELP, argv)
    parameter, values = conv.parse_scan(args.software, args.scan)
    generator = find_generator(args.software)
    options = JobOptions(task="sp", charge=args.charge, multiplicity=args.multiplicity,
                         kspacing=args.kspacing, fine_grid=args.fine_grid)
    options.check()
    if parameter == "fine-grid" and args.fine_grid is not None:
        raise UsageError("--fine-grid and --scan fine-grid both set the fine grid; give one")
    type_map = parse_type_map(args.type_map) if args.type_map else None
    files = expand_inputs(args.input)
    generator = generator(config.load())
    ensure_dir(args.output, args.mkdir)
    structures = [(name, traj, i) for name, traj, i, _ in read_structures(files, args, type_map)]
    try:
        built, manifest = conv.build(generator, args.software, parameter, values, structures, options,
                                     args.template)
    except (ValueError, RuntimeError) as e:
        raise UsageError(str(e)) from None
    report_batch(built, write_batch(built, args.output))
    path = Path(args.output) / conv.MANIFEST
    path.write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")
    print(f"Wrote: {path}  ({len(values)} values x {len(structures)} structures)")
    return 0


def _fixed(x, digits):
    return "" if math.isnan(x) else f"{x:.{digits}f}"


def run_conv_analyze(argv):
    import json

    from ferro.jobs import conv

    args = parse(conv_analyze_parser(), CONV_ANALYZE_HELP, argv)
    path = Path(args.dir) / conv.MANIFEST
    if not path.is_file():
        raise UsageError(f"{path} does not exist: DIR must be a directory written by `ferro-py conv`")
    manifest = json.loads(path.read_text(encoding="utf-8"))
    ref = args.ref
    if ref is not None and conv.SCANS[manifest["software"]][manifest["parameter"]][1] is float:
        try:
            ref = float(ref)
        except ValueError:
            raise UsageError(f"--ref {ref} is not a number") from None
    thresholds = {k: v for k, v in (("de", args.de), ("dp", args.dp), ("df_rms", args.df_rms),
                                    ("df_max", args.df_max)) if v is not None}
    rows, ref = conv.compare(manifest, args.dir, ref, thresholds)
    limits = {**conv.THRESHOLDS, **thresholds}
    out = Path(args.output) if args.output else Path(args.dir) / f"conv_{manifest['parameter']}.csv"
    columns = list(rows[0])
    meta = [f"ferro-py conv analyze: {manifest['software']} {manifest['parameter']}"
            + (f" [{manifest['unit']}]" if manifest["unit"] else ""),
            f"reference: {conv.label(ref)}",
            f"thresholds (hints): |dE| <= {limits['de']:g} meV/atom, |dP| <= {limits['dp']:g} kbar, "
            f"dF_rms <= {limits['df_rms']:g} meV/A, dF_max <= {limits['df_max']:g} meV/A"]
    ferro.write_table(str(out), [(c, [r[c] for r in rows]) if isinstance(rows[0][c], str)
                                 else (c, [float(r[c]) for r in rows]) for c in columns], meta)

    print(f"{manifest['parameter']} (reference {conv.label(ref)}), worst frame of each value:")
    print(f"  {'value':>10} {'frames':>7} {'dE meV/at':>10} {'dP kbar':>9} {'dF_max':>8} {'dF_rms':>8}  within")
    table = conv.summary(rows, manifest["values"])
    for s in table:
        print(f"  {s['value']:>10} {s['frames']:>7} {_fixed(s['dE'], 3):>10} {_fixed(s['dP'], 2):>9} "
              f"{_fixed(s['dF_max'], 1):>8} {_fixed(s['dF_rms'], 1):>8}  {'yes' if s['within'] else 'no'}")
    first = next((s["value"] for s in table if s["within"] and s["value"] != conv.label(ref)), None)
    print(f"First value within every threshold: {first}" if first else
          "No value below the reference is within every threshold")
    failed = [r for r in rows if r["status"] != "ok"]
    for r in failed:
        print(f"  {r['value']} {r['structure']}: {r['status']}")
    print(f"Wrote: {out}")
    return 1 if failed else 0


COMMANDS = {"config": run_config, "job": run_job, "conv": run_conv}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ("-h", "--help"):
        print(OVERVIEW)
        return 0
    command = COMMANDS.get(argv[0])
    try:
        if command is None:
            raise UsageError(f"unknown command '{argv[0]}'; run `ferro-py -h` for the list")
        return command(argv[1:])
    except (UsageError, ValueError, RuntimeError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
