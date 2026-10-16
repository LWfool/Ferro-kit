"""ferro-py 命令行：帮助页与 argparse 一致（同 ferro-cli 的 help_sync），以及各命令的行为。"""

import re

import pytest

from ferro import cli


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
PAGES = [(cli.config_parser, cli.CONFIG_INIT_HELP)]


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
