"""ferro-py 命令行：帮助页与 argparse 一致（同 ferro-cli 的 help_sync），以及各命令的行为。"""

import re
from pathlib import Path

import pytest

import ferro
from ferro import cli, jobs
from ferro.jobs.base import Generator, InputSet


def help_options(text):
    """帮助页 Parameters 段里列出的全部选项名。"""
    section = text.split("Parameters:\n", 1)[1].split("\n\n", 1)[0]
    options = set()
    for line in section.splitlines():
        first_column = re.split(r"\s{2,}", line.strip())[0]   # "-o, --output PATH"
        options.update(re.findall(r"(?<![\w-])--?\w[\w-]*", first_column))
    return options


def parser_options(parser):
    return {s for action in parser._actions for s in action.option_strings}


# 每个子命令：(argparse 构造函数, 帮助页)
PAGES = [(cli.config_parser, cli.CONFIG_INIT_HELP), (cli.job_parser, cli.JOB_HELP),
         (cli.job_template_parser, cli.JOB_TEMPLATE_HELP), (cli.conv_parser, cli.CONV_HELP),
         (cli.conv_analyze_parser, cli.CONV_ANALYZE_HELP), (cli.doc_parser, cli.DOC_HELP)] + \
        [(lambda figure=figure: cli.plot_parser(figure), text) for figure, text in cli.PLOT_HELPS.items()]


@pytest.mark.parametrize("make_parser, text", PAGES)
def test_help_pages_match_parsers(make_parser, text):
    sections = [line for line in text.splitlines() if line and not line.startswith(" ")]
    assert sections[1:] == ["Parameters:", "Output:", "Examples:",
                            "Full documentation:  " + sections[-1].split("  ", 1)[1]], \
        "帮助页五段，且只有五段"
    assert help_options(text) == parser_options(make_parser()), "帮助页与 argparse 的选项应一致"


def test_overview_and_unknown_command(capsys):
    assert cli.main([]) == 0
    assert "config init" in capsys.readouterr().out
    assert cli.main(["plott"]) == 1
    assert "ferro-py: FerroError: unknown command 'plott'" in capsys.readouterr().err


# ── ferro-py doc ─────────────────────────────────────────────────────────────

MANUAL = Path(__file__).resolve().parents[2] / "docs" / "src"


def test_every_help_page_points_at_a_manual_page():
    # 同 ferro doc 的 PAGES 登记：帮助页末行的 `ferro-py doc X` 必须真有那一页
    topics = {topic for topic, _, _ in ferro._ferro._doc_pages()}
    for text in [cli.OVERVIEW, cli.PLOT_HELP] + [text for _, text in PAGES]:
        asked = text.rsplit("Full documentation:  ferro-py doc", 1)[1].strip()
        assert asked == "" or asked in topics, f"帮助页指向不存在的手册页 `{asked}`"


def test_manual_pages_match_the_book():
    pages = ferro._ferro._doc_pages()
    assert len({t for t, _, _ in pages}) == len(pages), "主题不得重名"
    summary = (MANUAL / "SUMMARY.md").read_text(encoding="utf-8")
    for topic, source, text in pages:
        assert text == (MANUAL / source).read_text(encoding="utf-8"), f"{topic}：编进去的不是当前源文件"
        assert f"]({source})" in summary, f"{source} 不在 SUMMARY.md 里"
    registered = {source for _, source, _ in pages}
    for path in (MANUAL / "python").glob("*.md"):
        assert f"python/{path.name}" in registered, f"{path.name} 没登记进 ferro-python/src/doc.rs 的 PAGES"


def test_doc_list_page_and_unknown(capsys):
    assert cli.main(["doc"]) == 0
    out = capsys.readouterr().out
    assert "job vasp" in out and "python/job-vasp.md" in out
    assert cli.main(["doc", "Job", "VASP"]) == 0, "主题不分大小写"
    assert capsys.readouterr().out == (MANUAL / "python" / "job-vasp.md").read_text(encoding="utf-8"), \
        "非终端输出 markdown 原文"
    assert cli.main(["doc", "jobs"]) == 1
    err = capsys.readouterr().err
    assert "no manual page for `jobs`" in err and "job cp2k" in err, "猜错时把列表给出来"


def test_config_init(tmp_path, monkeypatch, capsys):
    monkeypatch.setenv("FERRO_CONFIG", str(tmp_path / "a" / "c.toml"))
    assert cli.main(["config", "init"]) == 0, "默认位置的目录直接建"
    assert (tmp_path / "a" / "c.toml").exists()
    assert cli.main(["config", "init"]) == 1
    assert "--force" in capsys.readouterr().err
    assert cli.main(["config", "init", "--force"]) == 0

    target = tmp_path / "b" / "c.toml"
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    assert cli.main(["config", "init", "-o", str(target)]) == 1, "显式 -o 的缺失目录在非交互时要 --mkdir"
    assert "--mkdir" in capsys.readouterr().err
    assert cli.main(["config", "init", "-o", str(target), "--mkdir"]) == 0
    assert target.exists()


# ── ferro-py job ─────────────────────────────────────────────────────────────

from conftest import NVT


class Fake(Generator):
    """POTCAR 只取决于元素顺序，同组成的帧共用一份；in.txt 带上帧的第一个坐标。"""
    program = "fake"
    template_name = "in.txt"

    def build(self, one, options, template_text):
        return InputSet({"in.txt": template_text + f"{one.positions(0)[0]}\n",
                         "POTCAR": " ".join(one.elements()) + "\n"},
                        shared={"POTCAR"}, notes=[f"task {options.task}"])


@pytest.fixture
def job(tmp_path, monkeypatch):
    monkeypatch.setitem(jobs.GENERATORS, "fake", Fake)
    monkeypatch.setenv("FERRO_CONFIG", str(tmp_path / "none.toml"))
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    (tmp_path / "t.txt").write_text("template\n")
    return lambda *a: cli.main(["job", "-s", "fake", "-t", "sp", "--template", str(tmp_path / "t.txt"), *a])


def test_job_trajectory_shares_potcar(tmp_path, job, capsys):
    out = tmp_path / "out"
    assert job("-i", str(NVT), "--number", "3", "-o", str(out), "--mkdir") == 0
    dirs = sorted(p.name for p in out.iterdir() if p.is_dir())
    assert dirs == ["70Z30P00A_NVT_5_0000", "70Z30P00A_NVT_5_0002", "70Z30P00A_NVT_5_0004"], \
        "每帧一个子目录，帧号 0 起、4 位"
    link = out / dirs[0] / "POTCAR"
    assert link.is_symlink() and str(link.readlink()) == "../POTCAR", "子目录里是相对路径的链接"
    assert (out / "POTCAR").read_text() == link.read_text()
    assert (out / dirs[0] / "in.txt").read_text() != (out / dirs[1] / "in.txt").read_text()
    printed = capsys.readouterr()
    assert "POTCAR->../POTCAR" in printed.out
    # A generator's notes are messages, in the CLI's layout on stderr
    assert f"{out / dirs[0]}: FerroNote: task sp" in printed.err.splitlines()


def test_job_single_structure(tmp_path, job):
    ferro.write(ferro.read(str(NVT)).tail(1), str(tmp_path / "one.vasp"), format="vasp")
    out = tmp_path / "out"
    out.mkdir()
    assert job("-i", str(tmp_path / "*.vasp"), "-o", str(out)) == 0, "加引号的通配符由程序展开"
    assert str((out / "one" / "POTCAR").readlink()) == "../POTCAR" and (out / "POTCAR").is_file(), \
        "a single job has the same layout as several (N=1 is a case of N)"


def test_job_mixed_compositions_write_nothing(tmp_path, job, capsys):
    t = ferro.read(str(NVT)).tail(1)
    ferro.write(t, str(tmp_path / "a.vasp"), format="vasp")
    ferro.write(t, str(tmp_path / "c.vasp"), format="vasp")
    lines = (tmp_path / "c.vasp").read_text().splitlines()
    lines[5] = lines[5].replace(t.elements()[0], "S", 1)   # 第 6 行元素行：换一个元素
    (tmp_path / "c.vasp").write_text("\n".join(lines) + "\n")
    out = tmp_path / "out"
    out.mkdir()
    assert job("-i", str(tmp_path / "a.vasp"), str(tmp_path / "c.vasp"), "-o", str(out)) == 1
    assert "different POTCAR" in capsys.readouterr().err
    assert list(out.iterdir()) == [], "一个失败就什么都不写"


@pytest.mark.parametrize("extra, message", [
    (["--temperature", "300"], "only apply to an aimd task"),
    (["--type-map", "1:Zr"], "is not TYPE=ELEMENT"),
    (["-s", "gaussian"], "unknown software 'gaussian'"),
    (["-t", "aimd", "--ensemble", "npt_i"], "NVT only"),
    (["--stride", "2", "--number", "3"], "give stride or number"),
    (["--start", "3", "--end", "1"], "start 3 is past end 1"),
    (["--units", "lj"], "unknown LAMMPS units 'lj'"),
    (["--atom-style", "bond"], "unknown atom style 'bond'"),
])
def test_job_parameter_errors_come_first(tmp_path, job, capsys, extra, message):
    # 输入文件不存在，也应先报参数错
    assert job("-i", str(tmp_path / "missing.vasp"), "-o", str(tmp_path), *extra) == 1
    assert message in capsys.readouterr().err


@pytest.mark.parametrize("command", ["job", "conv"])
@pytest.mark.parametrize("option", ["--start", "--end", "--stride", "--number"])
def test_negative_frame_selection_is_a_usage_error(tmp_path, monkeypatch, capsys, command, option):
    monkeypatch.setenv("FERRO_CONFIG", str(tmp_path / "none.toml"))
    head = ["-t", "sp"] if command == "job" else ["--scan", "encut=500,600"]
    with pytest.raises(SystemExit) as e:
        cli.main([command, "-s", "vasp", *head, "-i", str(NVT), "-o", str(tmp_path), option, "-1"])
    assert e.value.code == 2
    assert f"argument {option}: -1 is negative" in capsys.readouterr().err


def test_job_needs_mkdir_without_terminal(tmp_path, job, capsys):
    assert job("-i", str(NVT), "-o", str(tmp_path / "new")) == 1
    assert "--mkdir" in capsys.readouterr().err


def test_failed_job_leaves_no_output_directory(tmp_path, job, capsys):
    out = tmp_path / "new"
    assert job("-i", str(tmp_path / "missing.vasp"), "-o", str(out), "--mkdir") == 1
    assert "missing.vasp" in capsys.readouterr().err
    assert not out.exists(), "-o is made only once every input has been read and generated"


def test_job_errors_name_the_frame(tmp_path, job, capsys):
    ferro.write(ferro.read(str(NVT)).tail(1), str(tmp_path / "mol.xyz"), format="xyz")
    assert job("-i", str(tmp_path / "mol.xyz"), "-o", str(tmp_path)) == 1
    assert "mol.xyz frame 0: only periodic systems" in capsys.readouterr().err


# -- ferro-py plot -----------------------------------------------------------

def test_plot_overview_and_unknown_figure(capsys):
    assert cli.main(["plot"]) == 0
    assert "net" in capsys.readouterr().out
    assert cli.main(["plot", "rdf", "-i", "x.csv"]) == 1
    assert "unknown figure 'rdf'; one of gr, sq, angle, msd, net" in capsys.readouterr().err


GR_CSV = """\
# ferro gr
file,r,center,neighbor,gr,cn
a,1.0,P,O,0.0,0.0
a,1.5,P,O,2.0,1.0
b,1.0,P,O,0.0,0.0
b,1.5,P,O,3.0,2.0
"""


def can_plot():
    import importlib.util
    import shutil
    return shutil.which("latex") is not None and all(
        importlib.util.find_spec(m) for m in ("matplotlib", "pandas", "scienceplots"))


@pytest.mark.skipif(not can_plot(), reason="needs matplotlib, pandas, scienceplots and LaTeX")
def test_plot_shared_legend_covers_every_panel(tmp_path):
    # the second csv has a trajectory the first lacks: it still needs a legend entry
    from ferro.plot import common, gr
    (tmp_path / "gr_P-O.csv").write_text(GR_CSV)
    (tmp_path / "gr_Si-O.csv").write_text(GR_CSV.replace("b,", "c,").replace("P,O", "Si,O"))
    with common.style():
        frames = common.read_all([str(tmp_path / "gr_P-O.csv"), str(tmp_path / "gr_Si-O.csv")])
        fig = gr.build_figure(frames, common.context(frames, tmp_path))   # labels.csv goes to tmp_path
    legends = [ax.get_legend() for ax in fig.axes if ax.get_legend() is not None]
    assert len(legends) == 1 and [t.get_text() for t in legends[0].get_texts()] == ["a", "b", "c"]


@pytest.mark.skipif(not can_plot(), reason="needs matplotlib, pandas, scienceplots and LaTeX")
def test_plot_gr_writes_figure_data_and_labels(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("sys.stdin.isatty", lambda: False)
    (tmp_path / "gr_P-O.csv").write_text(GR_CSV)
    out = tmp_path / "figs"
    assert cli.main(["plot", "gr", "-i", str(tmp_path / "gr_*.csv"), "-o", str(out)]) == 1, \
        "a missing -o needs --mkdir without a terminal"
    capsys.readouterr()
    assert cli.main(["plot", "gr", "-i", str(tmp_path / "gr_*.csv"), "-o", str(out), "--mkdir",
                     "-s", "run1"]) == 0
    printed = capsys.readouterr()
    assert sorted(p.name for p in out.iterdir()) == \
        ["gr_run1.pdf", "gr_run1.png", "gr_run1_data.csv", "labels.csv"], "-s is the last segment"
    assert f"Wrote: {out / 'gr_run1_data.csv'}" in printed.out
    assert f"{out / 'labels.csv'}: FerroNote: written" in printed.err, "notes in the CLI's layout"
    assert (out / "labels.csv").read_text() == "key,label,show\na,a,1\nb,b,1\n"
    data = (out / "gr_run1_data.csv").read_text().splitlines()
    assert data[0] == "panel,file,label,quantity,x,y" and len(data) == 1 + 2 * 2 * 2, \
        "two inputs x (gr, cn) x two points"


@pytest.mark.skipif(not can_plot(), reason="needs matplotlib, pandas, scienceplots and LaTeX")
def test_plot_net_names_the_wrong_table(tmp_path, capsys):
    (tmp_path / "gr.csv").write_text(GR_CSV)
    assert cli.main(["plot", "net", "qn", "-i", str(tmp_path / "gr.csv"), "-o", str(tmp_path)]) == 1
    assert "'qn' needs column(s) count, former, qn" in capsys.readouterr().err
