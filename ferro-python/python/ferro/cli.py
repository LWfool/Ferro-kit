"""``ferro-py``：ferro 的 Python 命令行。子命令 ``config``、``job``，以后加 ``plot``。

帮助页与 Rust 的 ``ferro`` 同为五段（用途 / Parameters / Output / Examples /
Full documentation），手写在各子命令旁；``tests/test_cli.py`` 核对帮助页与
argparse 的参数一致（同 ``ferro-cli`` 的 ``help_sync``）。
"""

import argparse
import sys
from pathlib import Path

from ferro import config

OVERVIEW = """\
ferro-py — the Python side of ferro

Usage:
  ferro-py <COMMAND> [OPTIONS]

Commands:
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


COMMANDS = {"config": run_config}


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
    except (UsageError, config.ConfigError, RuntimeError, OSError) as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
