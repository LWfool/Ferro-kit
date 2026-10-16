"""CP2K：输入节树的解析、编辑与写回（结构参照 ASE ``calculators/cp2k.py`` 的 ``InputSection``）。

节名、关键字都不分大小写（CP2K 如此），保存时保留原文写法。整行注释、空行、
行尾注释都挂在节树上随之写回；预处理指令（``@SET``、``@INCLUDE`` 等）原样保留、
不展开。写回统一两格缩进，所以不是逐字节还原，但再解析得到同一棵树。
"""


class Cp2kInputError(ValueError):
    pass


def _split_comment(line):
    """引号外第一个 ``#`` 或 ``!`` 起是注释；返回 (代码, 注释)，注释含起始符。"""
    quote = None
    for i, c in enumerate(line):
        if quote:
            if c == quote:
                quote = None
        elif c in "'\"":
            quote = c
        elif c in "#!":
            return line[:i].rstrip(), line[i:]
    return line.rstrip(), ""


def format_value(value):
    if isinstance(value, bool):
        return ".TRUE." if value else ".FALSE."
    if isinstance(value, (list, tuple)):
        return " ".join(format_value(v) for v in value)
    return str(value)


class Keyword:
    def __init__(self, name, value="", comment=""):
        self.name, self.value, self.comment = name, value, comment

    def to_lines(self, indent):
        text = f"{indent}{self.name} {self.value}".rstrip()
        return [f"{text}  {self.comment}" if self.comment else text]


class Line:
    """整行注释、空行或预处理指令：原样写回（去掉原缩进）。"""

    def __init__(self, text):
        self.text = text

    def to_lines(self, indent):
        return [f"{indent}{self.text}" if self.text else ""]


class Section:
    """一节：``&NAME param`` …… ``&END NAME``。根节点 ``name`` 为 None。"""

    def __init__(self, name=None, param="", comment=""):
        self.name, self.param, self.comment = name, param, comment
        self.items = []

    @classmethod
    def from_str(cls, text):
        root = cls()
        stack = [root]
        for n, raw in enumerate(text.splitlines(), 1):
            code, comment = _split_comment(raw.strip())
            if not code:
                stack[-1].items.append(Line(comment))
                continue
            word, rest = (code.split(None, 1) + [""])[:2]
            if word.upper() == "&END":
                if len(stack) == 1:
                    raise Cp2kInputError(f"line {n}: &END without an open section")
                if rest and rest.upper() != stack[-1].name.upper():
                    raise Cp2kInputError(
                        f"line {n}: &END {rest} closes &{stack[-1].name}")
                stack.pop()
            elif word.startswith("&"):
                section = cls(word[1:], rest, comment)
                stack[-1].items.append(section)
                stack.append(section)
            elif word.startswith("@"):
                stack[-1].items.append(Line(raw.strip()))
            else:
                stack[-1].items.append(Keyword(word, rest, comment))
        if len(stack) > 1:
            raise Cp2kInputError(f"section &{stack[-1].name} is not closed with &END")
        return root

    # ── 查找 ───────────────────────────────────────────────────────────────
    def sections(self, name):
        return [s for s in self.items if isinstance(s, Section) and s.name.upper() == name.upper()]

    def section(self, path, create=False):
        """``"FORCE_EVAL/DFT/MGRID"`` → 那一节（同名取第一个）；没有返回 None，``create`` 时新建。"""
        node = self
        for name in path.split("/"):
            found = node.sections(name)
            if found:
                node = found[0]
            elif create:
                new = Section(name.upper())
                node.items.append(new)
                node = new
            else:
                return None
        return node

    def keywords(self, name):
        return [k for k in self.items if isinstance(k, Keyword) and k.name.upper() == name.upper()]

    def get(self, name, default=None):
        found = self.keywords(name)
        return found[0].value if found else default

    # ── 编辑 ───────────────────────────────────────────────────────────────
    def set(self, name, value):
        """设关键字（同名取第一个，其余删去）；没有就加在本节子节之前。"""
        found = self.keywords(name)
        if found:
            found[0].value = format_value(value)
            for extra in found[1:]:
                self.items.remove(extra)
            return
        position = next((i for i, x in enumerate(self.items) if isinstance(x, Section)), len(self.items))
        self.items.insert(position, Keyword(name.upper(), format_value(value)))

    def remove(self, name):
        """删去这一名字的全部关键字与子节。"""
        self.items = [x for x in self.items
                      if not (isinstance(x, (Keyword, Section)) and x.name.upper() == name.upper())]

    def add(self, item):
        self.items.append(item)
        return item

    # ── 写回 ───────────────────────────────────────────────────────────────
    def to_lines(self, indent=""):
        if self.name is None:
            return [line for x in self.items for line in x.to_lines(indent)]
        head = f"{indent}&{self.name} {self.param}".rstrip()
        lines = [f"{head}  {self.comment}" if self.comment else head]
        for x in self.items:
            lines += x.to_lines(indent + "  ")
        lines.append(f"{indent}&END {self.name}")
        return lines

    def to_str(self):
        return "\n".join(self.to_lines()) + "\n"
