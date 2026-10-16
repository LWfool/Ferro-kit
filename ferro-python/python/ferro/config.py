"""全局配置：赝势、基组数据与模板目录的位置（`~/.config/ferro/config.toml`）。

位置可用环境变量 ``FERRO_CONFIG`` 覆盖。文件不存在时是空配置：只有真用到某一项
（比如生成 VASP 任务要 POTCAR）才报错，报错里指向 ``ferro-py config init``。
写错的键（``potcar_dri``）直接报错，不静默忽略 —— 否则用户以为配了、程序却在用默认。
"""

import os
import tomllib
from pathlib import Path

# 每节允许的键："path" 是路径（展开 ~ 与 $VAR，相对路径相对配置文件所在目录）；
# "elements" 是 元素 → 值 的表，值的类型在第三列
SCHEMA = {
    "vasp": {"potcar_dir": "path", "potcar": ("elements", str)},
    "qe": {"pseudo_dir": "path", "pseudo_table": "path", "pseudos": ("elements", str)},
    "cp2k": {
        "data_dir": "path",
        "basis_family": str,
        "potential_family": str,
        "kinds": ("elements", dict),
    },
    "templates": {"dir": "path"},
}

# CP2K 逐元素覆盖里允许的键
KIND_KEYS = ("basis", "potential")

TEMPLATE = """\
# ferro-py global configuration. Every entry is optional: an entry is only needed
# when a job uses it, and the error then names the entry to add.
# Paths may use ~ and $VARS; a relative path is relative to this file.

[vasp]
# One sub-directory per POTCAR: <potcar_dir>/<name>/POTCAR
# potcar_dir = "~/pp/vasp/potpaw_PBE"

[vasp.potcar]
# Element -> sub-directory name, for the elements that should not use <element>/
# Li = "Li_sv"

[qe]
# Directory holding the UPF files; written into pw.x's pseudo_dir
# pseudo_dir = "~/pp/qe/SSSP_1.3.0_PBE_precision"
# SSSP-style JSON {element: {filename, cutoff_wfc, cutoff_rho}}: file names and cutoffs
# pseudo_table = "~/pp/qe/SSSP_1.3.0_PBE_precision.json"

[qe.pseudos]
# Element -> UPF file name, overriding pseudo_table
# O = "O.pbe-n-kjpaw_psl.0.1.UPF"

[cp2k]
# Directory holding BASIS_MOLOPT, POTENTIAL, ... (CP2K's data/)
# data_dir = "~/cp2k/data"
# basis_family = "DZVP-MOLOPT-SR-GTH"
# potential_family = "GTH-PBE"

[cp2k.kinds]
# Per-element override of the basis set and/or the pseudopotential
# Li = { basis = "DZVP-MOLOPT-SR-GTH-q3", potential = "GTH-PBE-q3" }

[templates]
# Your own templates; a file here replaces the built-in one of the same name
# dir = "~/.config/ferro/templates"
"""


class ConfigError(ValueError):
    pass


def config_path():
    env = os.environ.get("FERRO_CONFIG")
    return Path(env) if env else Path.home() / ".config" / "ferro" / "config.toml"


def load(path=None):
    """读配置并校验，返回 ``{节: {键: 值}}``；没有文件时只有空的元素表。"""
    path = Path(path) if path is not None else config_path()
    # 元素表恒在（缺省为空表），用的地方直接查，不用先判有没有这一节
    config = {
        section: {key: {} for key, kind in keys.items() if isinstance(kind, tuple)}
        for section, keys in SCHEMA.items()
    }
    if not path.exists():
        return config
    try:
        raw = tomllib.loads(path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as e:
        raise ConfigError(f"{path}: {e}") from None
    for section, table in raw.items():
        if section not in SCHEMA:
            raise ConfigError(f"{path}: unknown section [{section}]; known: {', '.join(SCHEMA)}")
        if not isinstance(table, dict):
            raise ConfigError(f"{path}: '{section}' must be a [{section}] section")
        for key, value in table.items():
            kind = SCHEMA[section].get(key)
            if kind is None:
                raise ConfigError(
                    f"{path}: unknown key '{key}' in [{section}]; known: {', '.join(SCHEMA[section])}"
                )
            config[section][key] = _check(path, f"{section}.{key}", kind, value)
    return config


def _check(path, name, kind, value):
    if kind == "path":
        if not isinstance(value, str):
            raise ConfigError(f"{path}: {name} must be a path string")
        p = Path(os.path.expandvars(os.path.expanduser(value)))
        return p if p.is_absolute() else (path.parent / p)
    if isinstance(kind, tuple):
        _, value_type = kind
        if not isinstance(value, dict):
            raise ConfigError(f"{path}: [{name}] must be a table of element = value")
        for element, v in value.items():
            if not isinstance(v, value_type):
                raise ConfigError(f"{path}: {name}.{element} must be a {value_type.__name__}")
            if value_type is dict:
                unknown = set(v) - set(KIND_KEYS)
                if unknown or not all(isinstance(x, str) for x in v.values()):
                    raise ConfigError(
                        f"{path}: {name}.{element} takes string values for {' / '.join(KIND_KEYS)}"
                    )
        return dict(value)
    if not isinstance(value, kind):
        raise ConfigError(f"{path}: {name} must be a {kind.__name__}")
    return value


def require(config, section, key):
    """取一项必须有的配置；缺了就报错，指明要在哪个文件里加哪一行。"""
    value = config[section].get(key)
    if value is None:
        raise ConfigError(
            f"[{section}] {key} is not set in {config_path()} "
            f"(create the file with `ferro-py config init`)"
        )
    return value


def init(path=None, force=False):
    """写一份带注释的样板；已存在时不覆盖，除非 ``force``。返回写到的路径。"""
    path = Path(path) if path is not None else config_path()
    if path.exists() and not force:
        raise ConfigError(f"{path} already exists; pass --force to overwrite it")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(TEMPLATE, encoding="utf-8")
    return path
