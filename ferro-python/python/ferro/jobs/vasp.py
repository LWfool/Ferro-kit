"""VASP：INCAR 的解析、编辑与写回。

模板就是用户跑通过的 INCAR，所以写回要尽量保持原文：没被改的行原样输出（含注释、
对齐、空行），改过的键只换值、行尾注释留着，新加的键追加在末尾。
"""

import re


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

    def to_str(self):
        return "\n".join(self._lines) + "\n"
