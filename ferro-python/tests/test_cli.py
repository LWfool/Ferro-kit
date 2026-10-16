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
         (cli.conv_analyze_parser, cli.CONV_ANALYZE_HELP)]


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
    assert "unknown command 'plott'" in capsys.readouterr().err


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

NVT = Path(__file__).resolve().parents[2] / "tests" / "70Z30P00A_NVT_5.lammpstrj"


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
    printed = capsys.readouterr().out
    assert "POTCAR->../POTCAR" in printed and "task sp" in printed


def test_job_single_structure(tmp_path, job):
    ferro.write(ferro.read(str(NVT)).tail(1), str(tmp_path / "one.vasp"))
    out = tmp_path / "out"
    out.mkdir()
    assert job("-i", str(tmp_path / "*.vasp"), "-o", str(out)) == 0, "加引号的通配符由程序展开"
    assert not (out / "one" / "POTCAR").is_symlink(), "单个任务不拆出共用文件"
    assert not (out / "POTCAR").exists()


def test_job_mixed_compositions_write_nothing(tmp_path, job, capsys):
    t = ferro.read(str(NVT)).tail(1)
    ferro.write(t, str(tmp_path / "a.vasp"))
    ferro.write(t, str(tmp_path / "c.vasp"))
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
])
def test_job_parameter_errors_come_first(tmp_path, job, capsys, extra, message):
    # 输入文件不存在，也应先报参数错
    assert job("-i", str(tmp_path / "missing.vasp"), "-o", str(tmp_path), *extra) == 1
    assert message in capsys.readouterr().err


def test_job_needs_mkdir_without_terminal(tmp_path, job, capsys):
    assert job("-i", str(NVT), "-o", str(tmp_path / "new")) == 1
    assert "--mkdir" in capsys.readouterr().err


def test_job_errors_name_the_frame(tmp_path, job, capsys):
    ferro.write(ferro.read(str(NVT)).tail(1), str(tmp_path / "mol.xyz"))
    assert job("-i", str(tmp_path / "mol.xyz"), "-o", str(tmp_path)) == 1
    assert "mol.xyz frame 0: only periodic systems" in capsys.readouterr().err
