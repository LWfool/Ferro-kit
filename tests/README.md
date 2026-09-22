# 测试 fixture 的出处

逐字节保留原始输出，**不要手工编辑**：这些文件的价值正在于它们是真实程序
写出来的，改一个空格就可能把 reader 真正要对付的排版磨平。需要一个变体时
在测试里 `replace` 出来，不要动盘上的文件。

| 文件 | 来源 | 覆盖什么 |
|---|---|---|
| `43Z43P15A_NPT_5.lammpstrj` | LAMMPS NPT 运行 | 变胞轨迹 |
| `70Z30P00A_NVT_5.lammpstrj` | LAMMPS NVT 运行 | 定胞轨迹 |
| `triclinic_2frames.lammpstrj` | 手工构造 | dump 的三斜 `BOX BOUNDS` 换算 |
| `CHGCAR_2atoms` | VASP | 电荷密度网格 |
| `vasp_OUTCAR_2frames` | VASP AIMD | OUTCAR 逐帧量 |
| `vasp_vasprun_2frames.xml` | VASP AIMD | vasprun.xml 逐帧量 |
| `cp2k_md_3frames.out` | CP2K 2025.2 AIMD | MD 布局：xyz 块、`MD\| Step number` |
| `5Al_0003_1500K_f394.out` | CP2K 2025.2 单点 | 新布局：`FORCES\|`、`STRESS\| … [bar]`、`CELL_TOP\|` 干扰项 |
| `cp2k_sp_v61.out` | [cp2kdata][] 的 `tests/test_energy_force/v6.1/normal/output`（MIT） | 老布局：`energy (a.u.):`、`ATOMIC FORCES in [a.u.]`、无前缀的 ` STRESS TENSOR [GPa]`、`Fe1`/`Fe2` 两个 kind 名 |
| `Experiment.xlsx` | 实验数据 | 表格读取 |

[cp2kdata]: https://github.com/robinzyb/cp2kdata
