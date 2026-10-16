"""全局配置：样板可读回、路径解析、写错的键报错。"""

import pytest

from ferro import config


def test_missing_file_is_empty_config(tmp_path):
    c = config.load(tmp_path / "none.toml")
    assert c["vasp"] == {"potcar": {}} and c["templates"] == {}, "没有配置文件时只有空的元素表"
    with pytest.raises(config.ConfigError, match=r"\[vasp\] potcar_dir is not set.*config init"):
        config.require(c, "vasp", "potcar_dir")


def test_template_loads_as_empty_config(tmp_path):
    path = config.init(tmp_path / "sub" / "c.toml")
    assert config.load(path) == config.load(tmp_path / "none.toml"), "样板全是注释，读回应与无文件相同"
    with pytest.raises(config.ConfigError, match="already exists"):
        config.init(path)
    config.init(path, force=True)


def test_template_entries_are_all_known(tmp_path):
    # 把样板里的注释行全部取消注释，每一项都应通过校验 —— 样板与 SCHEMA 不漂
    text = "\n".join(
        line[2:] if line.startswith("# ") and "=" in line else line
        for line in config.TEMPLATE.splitlines()
    )
    path = tmp_path / "c.toml"
    path.write_text(text)
    c = config.load(path)
    assert c["vasp"]["potcar"] == {"Li": "Li_sv"}
    assert c["cp2k"]["kinds"]["Li"]["potential"] == "GTH-PBE-q3"
    for section, keys in config.SCHEMA.items():
        assert set(c[section]) == set(keys), f"样板应覆盖 [{section}] 的全部键"


def test_paths_resolve(tmp_path, monkeypatch):
    monkeypatch.setenv("PP_ROOT", "/data/pp")
    path = tmp_path / "c.toml"
    path.write_text('[vasp]\npotcar_dir = "$PP_ROOT/pbe"\n[qe]\npseudo_dir = "upf"\n'
                    '[templates]\ndir = "~/tpl"\n')
    c = config.load(path)
    assert str(c["vasp"]["potcar_dir"]) == "/data/pp/pbe", "展开环境变量"
    assert c["qe"]["pseudo_dir"] == tmp_path / "upf", "相对路径相对配置文件所在目录"
    assert not str(c["templates"]["dir"]).startswith("~"), "展开 ~"


@pytest.mark.parametrize("text, message", [
    ('[vasp]\npotcar_dri = "x"\n', "unknown key 'potcar_dri'"),
    ('[gaussian]\nx = 1\n', r"unknown section \[gaussian\]"),
    ('[vasp]\npotcar_dir = 3\n', "must be a path"),
    ('[vasp.potcar]\nLi = 1\n', "vasp.potcar.Li must be a str"),
    ('[cp2k.kinds]\nLi = { basis = "x", charge = "3" }\n', "takes string values"),
    ('[vasp\n', r"c\.toml"),
])
def test_bad_entries_are_errors(tmp_path, text, message):
    path = tmp_path / "c.toml"
    path.write_text(text)
    with pytest.raises(config.ConfigError, match=message):
        config.load(path)


def test_ferro_config_overrides_location(tmp_path, monkeypatch):
    monkeypatch.setenv("FERRO_CONFIG", str(tmp_path / "x.toml"))
    assert config.config_path() == tmp_path / "x.toml"
