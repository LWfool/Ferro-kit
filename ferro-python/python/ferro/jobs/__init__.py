"""计算输入生成：VASP / QE / CP2K × SP / opt / AIMD。

每个程序一个模块（``vasp``、``qe``、``cp2k``），以后加 Gaussian 就加一个模块。
模板就是原生输入文件，程序只改结构与少数几个量，其余原样保留。
"""

# 程序名 → Generator 子类；ferro-py job 的 -s 取值即此表的键
GENERATORS = {}
