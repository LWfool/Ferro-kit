"""Quantum ESPRESSO：pw.x 输入的解析、编辑与写回。

namelist 的解析与写出调 Rust（``ferro.parse_qe_input_parts`` / ``ferro.write_namelists``），
与 ``ferro.read`` 读 QE 结构是同一份代码，语法只实现一次。namelist 内的注释不保留；
卡片部分按行原样保留（含注释），只有被替换的卡片才重写。
"""

import ferro

# pw.x 的全部卡片名（INPUT_PW）。卡片在 Rust reader 里另有按名查找（readers/qe.rs
# find_card），那里只为读结构定位三张卡片，这里为编辑切出全部卡片，用途不同，不合并
CARDS = (
    "ATOMIC_SPECIES", "ATOMIC_POSITIONS", "K_POINTS", "ADDITIONAL_K_POINTS",
    "CELL_PARAMETERS", "CONSTRAINTS", "OCCUPATIONS", "ATOMIC_VELOCITIES",
    "ATOMIC_FORCES", "SOLVENTS", "HUBBARD",
)


class QeInputError(ValueError):
    pass


def format_value(value):
    """Python 值 → Fortran namelist 值：bool → ``.true.``，str 加单引号。"""
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, str):
        if "'" in value:
            raise QeInputError(f"a namelist string cannot hold a single quote: {value!r}")
        return f"'{value}'"
    if isinstance(value, (int, float)):
        return repr(value)
    raise QeInputError(f"cannot write {value!r} as a namelist value")


def unquote(raw):
    """原文值去掉引号（``'scf'`` → ``scf``）；不带引号的原样返回。"""
    raw = raw.strip()
    if len(raw) >= 2 and raw[0] == raw[-1] and raw[0] in "'\"":
        return raw[1:-1]
    return raw


def _split_cards(text):
    """卡片部分 → (卡片之前的行, [[卡片名, 卡片行...]])。"""
    head, cards = [], []
    for line in text.splitlines():
        word = line.split()[0].upper() if line.split() else ""
        if word in CARDS:
            cards.append([word, line])
        elif cards:
            cards[-1].append(line)
        else:
            head.append(line)
    return head, cards


class QeInput:
    """pw.x 输入。namelist 名与键都不分大小写，取值是原文（``'scf'``、``.true.``）。"""

    def __init__(self, namelists, head, cards):
        self.namelists = namelists  # {名小写: {键小写: 原文值}}，保持顺序
        self._head = head
        self._cards = cards

    @classmethod
    def from_str(cls, text):
        try:
            parsed, rest = ferro.parse_qe_input_parts(text)
        except RuntimeError as e:
            raise QeInputError(str(e)) from None
        namelists = {}
        for name, entries in parsed:
            if name in namelists:
                raise QeInputError(f"namelist &{name.upper()} appears twice")
            namelists[name] = dict(entries)
        head, cards = _split_cards(rest)
        return cls(namelists, head, cards)

    def get(self, namelist, key, default=None):
        return self.namelists.get(namelist.lower(), {}).get(key.lower(), default)

    def set(self, namelist, key, value):
        """设一项；namelist 不存在就在末尾新建。``value`` 按 Python 类型写成 Fortran。"""
        self.namelists.setdefault(namelist.lower(), {})[key.lower()] = format_value(value)

    def delete(self, namelist, key):
        self.namelists.get(namelist.lower(), {}).pop(key.lower(), None)

    def card(self, name):
        """一张卡片的原文（首行起，含尾随空行与注释）；没有返回 None。"""
        for c in self._cards:
            if c[0] == name.upper():
                return "\n".join(c[1:])
        return None

    def set_card(self, name, text):
        """替换一张卡片（``text`` 从卡片首行起）；没有就追加到末尾。"""
        name = name.upper()
        lines = text.rstrip("\n").splitlines() + [""]
        if lines[0].split()[0].upper() != name:
            raise QeInputError(f"card text must start with {name}")
        for c in self._cards:
            if c[0] == name:
                c[1:] = lines
                return
        self._cards.append([name] + lines)

    def remove_card(self, name):
        self._cards = [c for c in self._cards if c[0] != name.upper()]

    def to_str(self):
        out = ferro.write_namelists([(n, list(e.items())) for n, e in self.namelists.items()])
        body = self._head + [line for c in self._cards for line in c[1:]]
        while body and not body[0].strip():
            body.pop(0)
        return out + "\n" + "\n".join(body).rstrip("\n") + "\n"
