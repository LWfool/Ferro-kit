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
| `vasp651_sp/frame_{0011,0161,0194}/` | VASP 6.5.1 单点（r2SCAN-D3，302 原子），200 帧中挑 3 帧，`OUTCAR` 与 `vasprun.xml` 原样 | 单点布局；一目录一次运行（collect 的分组单位）。0011 的 xml 应力非对称最大，0161 力最大（13.76 eV/Å），0194 能量最低 |
| `cp2k_md_3frames.out` | CP2K 2025.2 AIMD | MD 布局：xyz 块、`MD\| Step number` |
| `5Al_0003_1500K_f394.out` | CP2K 2025.2 单点 | 新布局：`FORCES\|`、`STRESS\| … [bar]`、`CELL_TOP\|` 干扰项 |
| `cp2k_sp_v61.out` | [cp2kdata][] 的 `tests/test_energy_force/v6.1/normal/output`（MIT） | 老布局：`energy (a.u.):`、`ATOMIC FORCES in [a.u.]`、无前缀的 ` STRESS TENSOR [GPa]`、`Fe1`/`Fe2` 两个 kind 名 |
| `qe75_{sp,sp_spin,sp_noconv,md,vcmd}.out` | QE 7.5 pw.x，6 原子 AlPO₄（`exp_603719.cif`），PBE-D3，`verbosity = 'high'`；由 `ferro-python/scripts/qe_fixtures_hpc.sh` 在超算算出（2026-10-17），输入与作业脚本在 `examples/qe/pwx/` | 单点：总力后跟七组力分量；`nspin = 2`；`electron_maxstep = 3` 不收敛；SVR 5 步 md 与 5 步 vc-md 的「结构印在上一步末尾」 |
| `Experiment.xlsx` | 实验数据 | 不是测试 fixture：只被 `scripts/compare_sq_experiment.py` 当实验 S(q) 读来对照 |

[cp2kdata]: https://github.com/robinzyb/cp2kdata
