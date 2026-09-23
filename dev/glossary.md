# 术语表（翻译用）

> 手册英文化的约束。翻译时**一律照此表**；表里没有的词停下来问，不要自己发挥。
> 每条给出现次数与两条原句出处，按语境判断。确认后本表即生效。

确认方式：直接改右侧的英文，或在 `备注` 里写要求。

## A · 领域术语（有文献定论）

### 形成子 ×48 → **network former**

- `analysis/network.md:46` | **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
- `analysis/network.md:9` | Qn 分布 | `qn` | Qn 形成子的**同核**连接数 $n$（P–O–P），文献 $Q^n_m$ 的 $n$ |

### 修饰子 ×19 → **modifier**

- `cli-reference.md:679` | `network_coordination.csv` | 配位数分布（形成子 + 修饰子） |
- `analysis/network.md:48` | **修饰子** (modifier) | `--modifier` 点名，截断同样用 `--<M>-<L>=<Å>` | **只**参与配位数 |

### 配体 ×50 → **ligand**

- `analysis/network.md:46` | **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
- `cli-reference.md:678` | `network_ligand_type.csv` | 配体分类，`label` 读作 `Al-O_b-P` |

### 桥氧 ×17 → **bridging oxygen (BO)**

- `analysis/network.md:344` 所以自由氧就是 `O_f`、非桥氧是 `P-O_n`、三簇氧是 `O_t`。`former_a` / `former_b`
- `analysis/network.md:20` > （文献写作 Al[4] / Al[5] / Al[6]），它们仍然参与桥氧判定、`ligand_type` 与

### 非桥氧 ×5 → **non-bridging oxygen (NBO)**

- `analysis/network.md:344` 所以自由氧就是 `O_f`、非桥氧是 `P-O_n`、三簇氧是 `O_t`。`former_a` / `former_b`
- `analysis/network.md:83` > 表中没有任何 `O_n, Al` 行），于是两张分布表逐档相同。换一个 Al 带非桥氧的体系，

### 桥联 ×14 → **linkage**  · 表名就叫 linkage

- `analysis/network.md:46` | **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
- `analysis/network.md:433` | **原子** | `P_2`、`Al_4` | `linkage`、导出轨迹 | 桥联描述的是**原子之间**的连接，而一个 Qn 单元含多个原子 |

### 桥接数 ×6 → **number of bridges**

- `analysis/network.md:46` | **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
- `analysis/network.md:79` $\sum(\text{桥接数}) \neq 2 \times |\text{O\_b}|$——一个 `O_t` 被三边各数一次。

### 连接数 ×8 → **number of connections**  · 与桥接数的区别是三簇氧

- `analysis/network.md:9` | Qn 分布 | `qn` | Qn 形成子的**同核**连接数 $n$（P–O–P），文献 $Q^n_m$ 的 $n$ |
- `analysis/network.md:417` | Qn 形成子 | `P_0` `P_1` … `Si_4` | **n**（同元素连接数，P–O–P；异核连接不计入） |

### 配位数 ×30 → **coordination number**

- `analysis/network.md:46` | **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
- `analysis/network.md:94` | 其他形成子（Al, …） | **否** | **配位数** | `coordination` 表 |

### 配位 ×43 → **coordination**

- `analysis/network.md:46` | **形成子** (former) | `--<F>-<L>=<Å>` 的左侧 | 桥接数、配位数、配体分类、桥联统计 |
- `analysis/network.md:94` | 其他形成子（Al, …） | **否** | **配位数** | `coordination` 表 |

### 同核 ×8 → **homonuclear**

- `analysis/network.md:9` | Qn 分布 | `qn` | Qn 形成子的**同核**连接数 $n$（P–O–P），文献 $Q^n_m$ 的 $n$ |
- `analysis/network.md:93` | Qn 形成子（默认 B, P, Si） | 是 | $n$（同核连接数） | Qn 分布 |

### 异核 ×5 → **heteronuclear**

- `analysis/network.md:10` | 异核桥分解 | `qn_partner` | 上表再按伙伴元素拆一维，即 $Q^n(m\mathrm{Al})$ |
- `analysis/network.md:417` | Qn 形成子 | `P_0` `P_1` … `Si_4` | **n**（同元素连接数，P–O–P；异核连接不计入） |

### 位点 ×15 → **site**

- `cli-reference.md:103` | `-x` / `-y` / `-z` | 按 `Atom::label` 选（位点标签） |
- `cli-reference.md:680` | `network_linkage.csv` | 桥的连接情况：配体元素 + 两端位点状态 |

### 位点标签 ×7 → **site label**

- `cli-reference.md:103` | `-x` / `-y` / `-z` | 按 `Atom::label` 选（位点标签） |
- `analysis/sq.md:123` 位点标签对应的原子数往往不足以让它的 partial 显出信号。**库层的 `GroupBy::Label` 不动**

### 截断 ×30 → **cutoff**

- `cli-reference.md:452` | `--r-cut-ab` | 2.3 | 端 A 到中心 B 的截断 [Å] —— A 是 `-a`/`-x` 给的那个 |
- `cli-reference.md:453` | `--r-cut-bc` | 2.3 | 端 C 到中心 B 的截断 [Å] —— C 是 `-c`/`-z` 给的那个 |

### 晶胞 ×11 → **cell**  · 不是 unit cell,MD 盒子

- `cli-reference.md:238` | `Density` | **g/cm³** = Σ(原子质量) / 晶胞体积。质量优先取文件里的显式值，否则查元素表 |
- `cli-reference.md:324` | `--pbc` | (auto) | `xyz`, `z`, `none`；省略时从晶胞自动判断 |

### 玻璃 ×10 → **glass**

- `analysis/network.md:501` | ZnO–P₂O₅ 玻璃 | P | Zn | P-O: 2.4 Å，Zn-O: 2.6 Å |
- `analysis/network.md:63` | ≥3 | `O_t` | Tricluster，三配位配体（含 Al 玻璃中的常见结构） |

### 磷酸盐 ×5 → **phosphate**

- `analysis/network.md:104` 文献的扩展记号 $Q^n_m$（铝磷酸盐写 $Q^n(m\mathrm{Al})$，硼磷酸盐写 $Q^n(m\mathrm{B})$，
- `analysis/network.md:441` 带元素前缀（`P-Q2` 而非 `Q2`）是因为双 Qn 形成子体系（硼磷酸盐 B+P、铝硅酸盐

### 伙伴 ×12 → **partner**  · 表名 qn_partner

- `analysis/network.md:10` | 异核桥分解 | `qn_partner` | 上表再按伙伴元素拆一维，即 $Q^n(m\mathrm{Al})$ |
- `cli-reference.md:677` | `network_qn_partner.csv` | 同上按伙伴元素拆开，即 $Q^n(m\mathrm{Al})$ |

### 团簇 ×2 → **cluster**

- `cli-reference.md:600` 从一组 QE `pp.x` 电荷密度 cube 计算 Qn 团簇周围取向平均的电子密度。**不用 `-i`**，
- `cli-reference.md:576` ### `sdf` — 团簇 SDF

### 单点 ×10 → **single-point**  · CP2K ENERGY/ENERGY_FORCE

- `cli-reference.md:769` | `-i <FILE>...` | CP2K MD 日志 / CP2K 单点输出 / VASP OUTCAR / vasprun.xml，支持 glob |
- `data-model.md:177` | CP2K 单点 out reader | `&KIND` 的名字（`ATOMIC KIND INFORMATION` 块），**仅在它与元素不同时**才填。同元素多 kin

### 逐帧 ×10 → **per-frame**

- `cli-reference.md:877` | `-f, --f-max <EV_PER_A>` | 20.0 | 逐帧最大力**矢量模长**超过则删；0 关闭 |
- `cli-reference.md:878` | `-s, --s-max <GPA>` | 10.0 | 逐帧 9 个应力分量绝对值的最大值超过则删；0 关闭 |

### 成分 ×8 → **composition**  · 表名 composition

- `cli-reference.md:947` | `shuffle` | 同成分全部拼接 → 按 seed 打乱 → 按 `--set-size` 切 |
- `cli-reference.md:930` | `-o <DIR>` | | 输出根目录，一个成分一个子目录 |

### 组成 ×6 → **composition**  · 与「成分」同义?需区分

- `analysis/network.md:8` | **结构组成** | `composition` | 一物种一行：`P-Q2` `Al_4` `O_b` `Zn_4`，各占其元素的比例 |
- `cli-reference.md:675` | `network_composition.csv` | **结构组成一览**：`P-Q2` `Al_4` `O_b` `Zn_4`，各占其元素的比例（每元素求和为 1） |

### 应力 ×6 → **stress**

- `cli-reference.md:878` | `-s, --s-max <GPA>` | 10.0 | 逐帧 9 个应力分量绝对值的最大值超过则删；0 关闭 |
- `cli-reference.md:827` **MD** 要求 CP2K 把坐标、力、应力全部打到 `__STD_OUT__`，这样一个 out 文件

### 电荷 ×8 → **charge**

- `cli-reference.md:736` | `<输入stem>_BCF.dat` | Bader Charge File —— 逐 Bader 体积的电荷、体积、坐标 |
- `cli-reference.md:298` | `--charge` | (from file) | 覆盖体系总电荷（在自旋推断之前生效） |

### 三簇氧 ×1 → **tricluster oxygen**

- `analysis/network.md:344` 所以自由氧就是 `O_f`、非桥氧是 `P-O_n`、三簇氧是 `O_t`。`former_a` / `former_b`


## B · 写作词（无唯一译法，必须先定）

### 口径 ×11 → **statement**  · 用户定。「旧口径（总桥）vs 文献口径」= 同一个量的两种说法

- `analysis/network.md:124` | | 旧口径（总桥） | 文献口径 |
- `analysis/network.md:303` > P–O–Al 与配位数一侧（看 `Al_*` 与 `Zn_*` 的 `sd`）。旧口径把刚性骨架与涨落混在

### 约定 ×8 → **regulation**  · 用户定。与「口径」(statement) 分开：约定是规矩，口径是说法

- `data-model.md:176` | CIF / CP2K inp / QE reader | 各自的位点名（`O1`、`Fe1`）——注意这些**不合** `<元素>_<后缀>` 约定 |
- `analysis/network.md:172` - **`linkage` 展示列**：`Al_4-O-P_2` 这样的人可读形式，数字按各自约定（Qn 形成子取

### 判据 ×7 → **criterion**

- `cli-reference.md:885` | `--shuffle` | 关闭 | 写出前打乱，**在所有判据与抽帧之后** |
- `cli-reference.md:901` 报告三张表：`[funnel]` 逐步剩余、`[criteria]` 每条判据判坏多少及**独占**多少、

### 产物 ×43 → **output**  · 或 product/artifact

- `cli-reference.md:44` | `-s <SUFFIX>` | 批次标记，追加在产物名末尾：`<命令>[_<表>][_<label>]_<后缀>.csv` |
- `cli-reference.md:770` | `-o <DIR>` | 输出根目录，**可选**；不给则产物落在各输入目录的同级 |

### 主产物 ×5 → **primary output**

- `cli-reference.md:425` （`_sq` / `_xrd` / `_neutron`，只出规范半边）。主产物是两条 total（一行一个 $q$），
- `analysis/sq.md:107` 主产物是两条 total（一行一个 $q$），加权 partial $w_{ij}(q)\,S_{ij}(q)$ 是能逐点求和

### 堆叠 ×15 → **stack**

- `cli-reference.md:530` 3-D 空间分布图（Gaussian cube 格式）。**逐输入一个 `.cube`**，没有可堆叠的表也没有图，
- `cli-reference.md:918` 经与其余产物同一个 writer，自带 `#` 头与 `[inputs]` 清单。多 system 堆叠成一份，

### 长表 ×7 → **long table**

- `cli-reference.md:522` `--plot` **冻结在自查质量**，不会去追 matplotlib：数据是长表 csv，一行 seaborn 就是
- `cli-reference.md:402` **长表**：`file, r, center, neighbor, gr, cn`。类型进数据列，故元素集不同的轨迹可直接

### 清单 ×12 → **list**

- `analysis/sq.md:113` `[inputs]` 清单（`pandas.read_csv(comment="#")` 会丢掉）。`-o` 给目录，
- `cli-reference.md:918` 经与其余产物同一个 writer，自带 `#` 头与 `[inputs]` 清单。多 system 堆叠成一份，

### 名单 ×5 → **list**  · 与「清单」撞车

- `analysis/network.md:249` | `--qn E,E` | `B,P,Si` | 报 Qn 的形成子，逗号分隔。**替换**默认名单而非叠加 |
- `cli-reference.md:665` | `--qn E,E` | `B,P,Si` | 报 Qn 的形成子。**替换**默认名单而非叠加；点名非形成子或已被 `--modifier` 占用的元素会报错 |

### 批次 ×12 → **batch**  · -s/--suffix

- `cli-reference.md:44` | `-s <SUFFIX>` | 批次标记，追加在产物名末尾：`<命令>[_<表>][_<label>]_<后缀>.csv` |
- `analysis/network.md:243` | `-s <SUFFIX>` | — | 批次标记：`network_<表>_<后缀>.csv` |

### 落盘 ×2 → **write to disk**

- `cli-reference.md:912` | | 屏幕 | 落盘 |
- `cli-reference.md:910` 打印与落盘分开：

### 丢帧 ×2 → **dropped frames**

- `cli-reference.md:852` 单点缺应力**不算**丢帧 —— 没开 `STRESS_TENSOR` 的单点照样是好数据，只是没有
- `cli-reference.md:851` 丢帧三类，**始终计数**：SCF 未收敛 / 块截断（含力数与原子数不等）/ 组成不符。

### 对拍 ×3 → **cross-check**  · 跟参考实现逐点比

- `cli-reference.md:846` **与 dpdata 的对拍**：单点与 AIMD 逐项比过 cp2kdata 0.7.4（CP2K 2025.2）。
- `analysis/angle.md:113` 两列：整数直方图是与 `dump2analysis` 逐 bin 对拍的依据，只留 `p` 就对不了。

### 静默 ×7 → **silently**

- `cli-reference.md:838` 两代同形。2025–2026 静默读；**2024 及之前每个文件多打一行 `NOTE:` 点名版本，
- `cli-reference.md:5` `fe-traj` 会让旧脚本「跑成功」却吐出自己解析不了的 csv，静默坏数据比命令消失难查。

### 告警 ×7 → **warning**

- `cli-reference.md:618` | `--rmsd-warn` | `0.5` | 对齐 RMSD 告警阈值 [Å] |
- `cli-reference.md:594` | `--rmsd-warn` | 0.5 | RMSD 告警阈值 [Å] |

### 报错 ×17 → **error out**

- `cli-reference.md:45` | `--mkdir` | 不询问直接创建 `-o` 的目录。**非交互环境（脚本、CI）下必须给**，否则报错退出 |
- `cli-reference.md:886` | `--seed <N>` | 666 | `--shuffle` 的种子；不带 `--shuffle` 给它会报错 |

### 点名 ×6 → **declare**（`--modifier` 处）/ **name**（其余）  · 用户定。
「`--modifier` 点名」= 用 `--modifier` 声明修饰子元素 → declare；
「点名版本」「未点名配对」→ name

- `analysis/network.md:48` | **修饰子** (modifier) | `--modifier` 点名，截断同样用 `--<M>-<L>=<Å>` | **只**参与配位数 |
- `data-model.md:176` | CIF / CP2K inp / QE reader | 各自的位点名（`O1`、`Fe1`）——注意这些**不合** `<元素>_<后缀>` 约定 |

### 位点名 ×1 → **site name**  · 用户指出：n-gram 把它切成了「点名」，是断词错误

- `data-model.md:176` | CIF / CP2K inp / QE reader | 各自的位点名（`O1`、`Fe1`）——注意这些**不合** `<元素>_<后缀>` 约定 |

### 补零 ×5 → **pad with zeros**

- `cli-reference.md:816` `T = 2·E_kin/(3N·k_B)` 反算，`#` 头会注明。缺失值渲染成**空字段**，不补零。
- `cli-reference.md:60` **缺的留空（NaN），不补零、不插值**。失败的输入被跳过、在输出的 `[inputs]` 块里留下

### 逐输入 ×5 → **per-input**

- `cli-reference.md:64` 产物是逐输入的命令（`ferro map` 的 cube、`ferro net --export-traj` 的轨迹）例外：
- `cli-reference.md:530` 3-D 空间分布图（Gaussian cube 格式）。**逐输入一个 `.cube`**，没有可堆叠的表也没有图，

### 逐原子 ×5 → **per-atom**

- `cli-reference.md:735` | `<输入stem>_ACF.dat` | Atomic Charges File —— 逐原子 Bader 电荷、体积、到表面的最小距离 |
- `cli-reference.md:596` 产物：逐原子类型一个 `<stem>_<label>.cube`（多族时 `<stem>_fam<N>_<label>.cube`）。

### 抽帧 ×5 → **frame sampling**

- `cli-reference.md:885` | `--shuffle` | 关闭 | 写出前打乱，**在所有判据与抽帧之后** |
- `cli-reference.md:893` 全部帧 → |F|max → |σ|max → min d(O-O) → Al6 → [区间/抽帧] → shuffle

### 打乱 ×4 → **shuffle**

- `cli-reference.md:947` | `shuffle` | 同成分全部拼接 → 按 seed 打乱 → 按 `--set-size` 切 |
- `cli-reference.md:885` | `--shuffle` | 关闭 | 写出前打乱，**在所有判据与抽帧之后** |

### 去重 ×5 → **deduplicate**

- `cli-reference.md:80` | `traj msd` / `vacf` / `vanhove` | `--elements`，**排序去重** | `msd_O-P.csv`、无筛选 `msd_all.c
- `cli-reference.md:818` 文件按首个 `MD| Step number` 排序，文件内保持原序。**重叠帧不去重**（重启只

### 变号 ×6 → **sign flip**

- `data-model.md:81` | DeePMD 的 `virial.npy` | eV，正 = 压缩 | `stress × V`，不变号 |
- `data-model.md:79` | extxyz 的 `stress=` | eV/Å³，正 = 拉伸 | 读写两侧**变号** |

### 折进 ×5 → **fold into**  · dump 的 element 列

- `analysis/network.md:456` | 标签存放 | 折进 `element` 列 | 独立的 `label:S:1` 列 |
- `data-model.md:182` **writer 从不自作主张把 `label` 折进元素列。** `ferro convert` 无论标签如何都写干净的

### 折叠 ×5 → **fold**

- `data-model.md:183` 元素符号。只有 `ferro net --export-traj` 折叠，且只对 LAMMPS dump——那个格式只有一列
- `analysis/network.md:462` **只有 `ferro net --export-traj` 会折叠**——`ferro convert` 无论标签如何都写干净的

### 选型 ×5 → **type selection**

- `analysis/network.md:437` > **下游选型认原子词汇。** 导出轨迹里那个 P 叫 `P_2`，不叫 `P-Q2`——轨迹标签必须能
- `analysis/network.md:479` > **按标签选型只在单帧成立。** `g(r)` 要求逐类型的粒子数守恒，而标签是动态的——

### 诊断 ×7 → **diagnostics**

- `cli-reference.md:771` | `--type <WHAT>` | `deepmd`（DeePMD system）\| `inspect`（只出诊断，不出数据集）[deepmd] |
- `cli-reference.md:108` partial 是能加回 total 的诊断分解（$\sum w_{ij}S_{ij} = \mathrm{total}$），只留一对

### 规范序 ×1 → **canonical order**  · 用户释义：规范化之后的序列/排序。
实现是按 `(Z, 元素符号)` 排（`merge.rs`），dpdata 用字母序 —— 两边都靠 `type_map.raw` 自描述

- `cli-reference.md:941` 组内各 system 的原子排列可以不同：合并时统一到规范序 `(Z, 符号)`，**逐原子

### 规范半边 ×3 → **canonical half**  · ⚠ 待用户裁决，见下方证据

定义句在 `analysis/network.md:177`：「桥联无方向，两端按 `(元素, 同核连接数, CN)`
排序后小的在前，**每对只存一次**」—— 指**元素对只输出一半**（P-O 出，O-P 不出），
与 `_sq`/`_xrd`/`_neutron` 三种加权是同一句话里的两件事。
`sq.md:105` 的下半句写着理由：「$S(q)$ 没有有向的对应物」（g(r) 有 CN(r) 那个有向量，
S(q) 没有，所以不必写满 n² 个有序对）

- `analysis/sq.md:105` （`<pair>_sq` / `_xrd` / `_neutron`，只出规范半边——$S(q)$ 没有有向的对应物）。
- `cli-reference.md:425` （`_sq` / `_xrd` / `_neutron`，只出规范半边）。主产物是两条 total（一行一个 $q$），

### 实测 ×5 → **measured**

- `analysis/network.md:488` - `traj gr -x Al_5 -y O_b --last-n 1` 实测参考轨迹得 4.933 而非 5.000：差的 0.067
- `analysis/network.md:445` 按非桥配体数 0/1/2/≥3 分档，而修饰子的实际配位数在 3–6，实测参考轨迹 97% 落进

### 词汇 ×9 → **vocabulary**  · 两套标签词汇

- `analysis/network.md:430` | 词汇 | 形如 | 用在 | 为什么 |
- `cli-reference.md:686` 标签有**两套词汇**：分布表（`composition` / `qn` / `qn_partner`）用**单元**词汇

### 兜底 ×2 → **fallback**

- `analysis/network.md:446` 兜底桶，没有分辨力。
- `cli-reference.md:834` 按 a.u. 兜底。

### 守卫 ×4 → **guard**

- `analysis/network.md:482` > 「per-type atom counts change」错误拒绝，这是守卫而非缺陷。多帧请按元素选。
- `analysis/network.md:463` 元素符号。折叠还带守卫：标签不形如 `<元素>_…` 时不折并计数告警（CIF 的 `O1`、CP2K

### 退化 ×3 → **degenerate**

- `analysis/network.md:24` > `ferro net`；导出退化为开关 `--export-traj`。旧标签（`P0` / `Of` / `On_P` /
- `cli-reference.md:247` 行过短退化出的 `X`）在 `effective_mass()` 里回退成 1 amu，除了数值变小之外

### 闭合 ×4 → **consistent**  · 两张表 count 闭合

- `analysis/network.md:135` - `qn` 就是 `qn_partner` 对 `m_` 各列的**边际**（count 与 fraction 精确闭合）
- `cli-reference.md:109` 恰好把这条闭合藏起来；要看某一对在 pandas 里选列即可。按 label 分辨的 partial 一并

### 剥掉 ×4 → **strip**

- `cli-reference.md:754` 产物带 `.db`（原始收集数据），它**不是**划分后缀，被 `filter` / `merge` 剥掉后
- `cli-reference.md:799` `.db` 进去会被一起拒掉。`filter` / `merge` 在命名自己的产物前剥掉它，于是

