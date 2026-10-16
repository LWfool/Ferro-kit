"""ferro — Python 接口：Rust 扩展（结构读写、结构操作、分析）+ 纯 Python 子模块。

Rust 部分编译为 ``ferro._ferro``，这里整体再导出，``import ferro`` 的用法不变。
"""

from ferro._ferro import *  # noqa: F401,F403
from ferro._ferro import __version__  # noqa: F401  （下划线开头，星号导入不带）
