# Pe 控制的碳酸钙溶蚀模式与震电响应机制分析

## 一针见血结论

三组 Pe 溶蚀模拟的最终孔隙率、最终渗透率和最终曲折度几乎收敛到同一个末态，但震电响应差别达到数十到上千倍。这说明震电响应不是由“最终溶掉多少”控制，而是由溶蚀路径控制。

最核心的机制是：

> Pe 改变酸输运和通道化路径，通道化一方面提前增强水力连通性，另一方面使 H+ 提前突破并显著提高孔隙流体电导率；H+ 突破使 zeta 电位接近零甚至变号，同时动态电导率暴涨，导致有效电动耦合效率 `L_abs/sigma_abs` 塌陷，最终削弱或翻转 interface EM response。

换句话说，高 Pe 并不是简单“通道越强，震电越强”。在这些数据里，高 Pe 的通道化更像是一次 electrokinetic short-circuit：水力上更通了，但电动源项被酸化和高电导屏蔽掉了。

## 论文阶段诊断

Current stage: Stage 1-2，Research Story Discovery / Manuscript Architecture Design。

Core contradiction: 结果不能写成“把 RT 输出接到震电正演模型”。顶级期刊需要一个机制主角。

Recommended story: `dissolution-regime-controlled electrokinetic current imbalance`。

Currently forbidden: 不能只用孔隙率、渗透率或峰值幅度单变量解释；不能把 `phi > 0.95` 后的孔弹模型外推结果当作可靠主证据；不能把 H+ 到电导率/zeta 的半经验映射说成实测。

This round's goal: 从 Pe=0.1、1、10 的数据中识别真正控制震电响应差异的参数组合和机制链条。

## 数据来源与对应关系

溶蚀源数据：

- Pe=0.1: `dissolution_results-Da_0.0369_Pe_0.1000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random`
- Pe=1: `dissolution_results-Da_0.0369_Pe_1.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random`
- Pe=10: `dissolution_results-Da_0.0369_Pe_10.0000_L_0.0010_lengthXAxis_0.060000_lengthYAxis_0.040000_random`

震电结果：

- Pe=0.1: `se_results_offset0p1`
- Pe=1: `se_results_offsetPep1`
- Pe=10: `se_results_offset10`

关键输入列来自 `global_evolution_log.csv`：

`timestep`, `time_s`, `porosity`, `permeability_mD`, `k_k0`, `avg_dissolution_rate`, `surface_area_cm2`, `grain_volume_cm3`, `injected_pv`, `outlet_H_conc`, `tortuosity`。

关键震电列来自 `seismoelectric_timeseries_results.csv`：

`Time_s`, `Porosity_raw`, `valid_poroelastic`, `Permeability_mD`, `Tortuosity`, `cH_molL`, `pH`, `zeta`, `omega_t`, `L_abs`, `sigma_abs`, `RE_abs`, `TTM_abs`, `Amax_waveform_spectral`。

## 现象：同一末态，不同路径，不同震电响应

三组最终状态相近：

| Pe | 终止孔隙率 | 终止 `k/k0` | 终止曲折度 | 终止 `outlet_H_conc` |
| --- | ---: | ---: | ---: | ---: |
| 0.1 | 1.0000 | 247.90 | 1.0037 | 5.73e-5 |
| 1 | 1.0000 | 247.90 | 1.0037 | 9.27e-5 |
| 10 | 1.0000 | 247.90 | 1.0037 | 9.62e-5 |

但路径差别非常大：

| 指标 | Pe=0.1 | Pe=1 | Pe=10 |
| --- | ---: | ---: | ---: |
| 全局时间步 | 44 | 62 | 87 |
| 完全溶蚀时间 | 14319.2 s | 1350.6 s | 406.4 s |
| `Hout >= 1e-5` | 13419.2 s | 500.6 s | 11.4 s |
| `k/k0 >= 100` | 12519.2 s, phi=0.977 | 600.6 s, phi=0.882 | 101.4 s, phi=0.894 |
| 最大 ICN | 0.635 | 0.996 | 0.991 |

这已经说明：Pe 控制的是溶蚀路径和突破时序，而不是最终孔隙率。

## 溶蚀模式解释

### Pe=0.1：扩散/反应前缘推进，弱通道化

Pe=0.1 中，H+ 主要在入口附近推进，出口 H+ 长时间接近背景值。到 `phi≈0.8` 时，`k/k0≈2.93`，`outlet_H_conc≈7.0e-14`，说明样品内部已经溶蚀不少，但酸还没有形成有效的贯穿突破。

机制含义：

- 溶蚀更接近 face/front dissolution。
- 水力连通性增强较慢。
- 孔隙流体仍保持低电导、zeta 强负值。
- 动电耦合仍强，因此 interface EM response 保持大幅值。

### Pe=1：通道化开始主导，水力突破早于完全溶蚀

Pe=1 中，ICN 接近 1，说明已经出现强连通通道。到 `phi≈0.8` 时，`k/k0≈23.6`，比 Pe=0.1 同孔隙率高约 8 倍；同时 H+ 已经升到 `5.6e-6`，开始明显酸化。

机制含义：

- 水力通道已经形成。
- H+ 已开始突破但 zeta 还没有完全翻转。
- 震电响应明显降低，但尚未像 Pe=10 那样完全进入强屏蔽状态。

### Pe=10：快速 channeling / acid breakthrough，电动源项短路

Pe=10 中，`Hout >= 1e-5` 在 11.4 s 就发生。这个时刻孔隙率只有约 0.605，但 H+ 已经足以让 zeta 接近零：

| 同一孔隙率附近 | Pe=0.1 | Pe=10 |
| --- | ---: | ---: |
| phi | 0.6039 | 0.6055 |
| `k/k0` | 1.30 | 2.14 |
| `outlet_H_conc` | ~0 | 1.01e-5 |
| `zeta` | -0.065 | 2.33e-5 |
| `L_abs` | 1.31e-9 | 3.82e-13 |
| `sigma_abs` | 6.36e-3 | 7.18e-2 |
| `RE_abs` | 1.40e9 | 6.30e4 |
| `Amax_waveform_spectral` | 9.37e15 | 3.88e11 |

这是整套结果里最强的机制证据：孔隙率几乎一样，渗透率只增加到 1.65 倍，但 `RE_abs` 降低约 2.2e4 倍，波形峰值降低约 2.4e4 倍。原因不是结构变化本身，而是酸突破导致 `zeta` 归零和 `sigma` 暴涨。

## 真正的主控参数

### 第一主控：H+ / pH / zeta / conductivity 组合

同一孔隙率下，Pe=10 的 H+ 远高于 Pe=0.1，使 `zeta` 从强负值变成接近零或正值，并使 `sigma_abs` 增加数十倍。

在 `phi≈0.75` 的空间快照中：

| Pe | `cH_molL` | `zeta` | `L_abs/sigma_abs` | `RE_abs` | `TTM_abs` | `Amax` | 空间峰值 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0.1 | 1.00e-7 | -6.50e-2 | 1.78e-7 | 9.57e8 | 1.27e-6 | 7.81e15 | 3.65e15 |
| 1 | 1.76e-3 | -8.15e-3 | 4.03e-9 | 3.17e7 | 1.29e-9 | 2.23e14 | 1.08e14 |
| 10 | 1.74e-2 | +1.60e-3 | 1.30e-10 | 9.23e5 | 5.55e-11 | 9.30e12 | 4.39e12 |

Pe=0.1 / Pe=10 在相近孔隙率快照下的差异：

- `RE_abs`: 约 1.0e3 倍。
- `TTM_abs`: 约 2.3e4 倍。
- 波形峰值: 约 8.4e2 倍。
- 空间峰值: 约 8.3e2 倍。
- `L_abs/sigma_abs`: 约 1.37e3 倍。

因此，最适合写成论文主控指标的不是 `L_abs` 或 `sigma_abs` 单独，而是：

> dynamic electrokinetic efficiency: `L_abs/sigma_abs`

它表达了“单位导电泄漏/屏蔽背景下还剩多少可转化的电动耦合”。

### 第二主控：连通性增强的时机

渗透率不是无关变量，而是上游过程的重要中介。它控制酸前缘是否能快速贯穿样品，并决定声波诱导的流固相对运动和动态渗透率。

但在这批结果里，渗透率不是造成幅值差异的第一原因。证据是同一孔隙率的 one-at-a-time 替换：

| 对 Pe=0.1 基准，只换 Pe=10 参数 | `phi≈0.60` 时 `RE_abs` 比值 |
| --- | ---: |
| 只换孔隙率 | 1.00 |
| 只换渗透率 | 0.78 |
| 只换曲折度 | 1.02 |
| 只换 H+/电化学 | 5.70e-5 |
| 结构参数全换 | 0.80 |
| 结构 + H+/电化学全换 | 4.50e-5 |

结论非常明确：渗透率和连通性解释了酸为什么能早突破；但震电幅值崩塌本身主要由电化学分支控制。

### 第三主控：zeta 变号导致极性翻转

空间峰值诊断显示：

| Pe | R_E 侧峰值 | T_E 侧峰值 | 极性 |
| --- | ---: | ---: | --- |
| 0.1 | -3.65e15 at -37 mm | +2.74e15 at +37 mm | R_E 负、T_E 正 |
| 1 | -1.08e14 at -37 mm | +4.21e13 at +37 mm | R_E 负、T_E 正 |
| 10 | +4.39e12 at -37 mm | -5.24e11 at +36 mm | 相对前两组翻转 |

Pe=10 的 zeta 在快照时已经为正，空间极性也翻转。这是非常强的机制点：interface EM response 不只是幅值指标，还可以作为溶蚀诱发电动源项变号的诊断。

## 可投稿的核心机制图景

可以把结果组织为一个三阶段模型：

### Stage A: Electrokinetically active front

低 Pe 或早期阶段，H+ 尚未突破，流体电导率低，zeta 保持强负值。虽然渗透率增强有限，但 `L_abs/sigma_abs` 高，界面电流不平衡强，interface EM response 大。

### Stage B: Hydraulic channeling

中等 Pe 或中期阶段，通道形成，渗透率和连通性快速增强。这个阶段对震电响应有双重作用：增强流固相对运动，但也给酸突破提供高速通道。

### Stage C: Electrochemical short-circuit / coupling collapse

高 Pe 或酸突破后，H+ 进入主流通道并到达出口。流体电导率暴涨，zeta 接近零或变号，动态电动耦合效率 `L_abs/sigma_abs` 塌陷。最终即使介质更高渗、更连通，interface EM response 反而减弱甚至极性翻转。

这条机制链条比“渗透率影响震电响应”更有创新性，因为它解释了一个反直觉现象：

> 溶蚀增强水力连通性并不必然增强震电响应；当连通性增强伴随酸突破和高电导屏蔽时，震电响应会被削弱。

## 创新点提炼

### 创新点 1：路径依赖，而不是终态依赖

三组最终孔隙率和渗透率接近，但震电响应差别巨大。说明 interface EM response 对溶蚀路径、酸突破时机和电化学状态敏感，而不是只对最终孔隙率敏感。

可写成：

> Seismoelectric interface responses retain a memory of dissolution pathway rather than simply tracking final porosity or permeability.

### 创新点 2：高渗通道的双重角色

通道化同时是水力增强机制和电动耦合破坏机制。它增强流体连通，但也把高 H+ / 高电导流体提前带到界面和出口，导致电动源项塌陷。

可写成：

> Channeling amplifies hydraulic connectivity but can suppress electrokinetic conversion by accelerating acid breakthrough and electrical screening.

### 创新点 3：`L_abs/sigma_abs` 作为机制指标

`L_abs` 表示可转化的电动耦合，`sigma_abs` 表示导电泄漏/屏蔽背景。二者比值比单独的渗透率或孔隙率更接近震电响应的物理本质。

可写成：

> The ratio between dynamic electrokinetic coupling and dynamic conductivity provides a compact diagnostic of dissolution-induced conversion efficiency.

### 创新点 4：zeta 零点是极性开关

Pe=10 的 zeta 变号和空间 signed peak 翻转相互对应。这个结果可以把“幅值变化”升级为“电动源项极性重构”。

可写成：

> A zeta-potential zero crossing acts as a polarity switch for the finite-offset interface EM response.

### 创新点 5：反应输运慢时间与震电波形快时间的耦合

论文可以强调两个时间尺度：

- dissolution time: 秒到小时，控制 `phi, k0, tau, H+`。
- waveform time: 微秒，控制 interface EM 到时、空间峰值和极性。

你的工作把慢时间上的反应输运路径映射到快时间上的震电界面波形，这正是 gap。

## 推荐图件序列

下面这些图已经生成到：

`paper_figure_sequence_analysis/pe_mode_mechanism/`

同时提供 PNG 和 PDF，其中快照拼图只有 PNG。

### 已生成 Figure A: 机制概念图

文件：

- `fig1_mechanism_cartoon.png`
- `fig1_mechanism_cartoon.pdf`

图上讨论：

这张图的主线不是“Pe 越大响应越小”，而是一个竞争机制。左侧三种 Pe 首先进入两个分支：水力分支和电化学分支。水力分支中，通道化提高 `k/k0`、降低曲折度，理论上可能增强流固相对运动；但电化学分支中，高 Pe 使 H+ 提前到达主通道和出口，导致 zeta 接近零或变号，并使 `sigma` 增加。右侧结果表明，在本数据中电化学抑制强于水力放大，所以最终得到弱化甚至极性翻转的 interface EM response。

这张图适合作为论文的机制总览图，放在 Results 开头或 Discussion 开头。

### 已生成 Figure B: 同孔隙率参数桥梁图

文件：

- `fig2_same_porosity_parameter_bridge.png`
- `fig2_same_porosity_parameter_bridge.pdf`

图上讨论：

这张图直接回答“根本参数是什么”。六个面板显示，在相同孔隙率坐标下，三组 Pe 的终态相近，但路径完全不同。

最重要的三点：

1. `k/k0` 面板显示 Pe=1 和 Pe=10 更早进入高渗状态，说明通道化确实先改变了水力连通性。
2. `H+`, `zeta`, `sigma` 面板显示 Pe=10 在很低孔隙率阶段就发生酸突破，zeta 很快从强负值走向零/正值，`sigma` 同时增加近两个数量级。
3. `|L|/|sigma|` 和 `Amax` 面板几乎同步塌陷，说明震电峰值不是直接追随 `k/k0` 增强，而是追随电动耦合效率的塌陷。

因此，这张图建议作为主结果图：它把 RT 输出、动态电动参数和震电响应放在同一个物理链条里。

### 已生成 Figure C: one-at-a-time 参数归因图

文件：

- `fig3_same_phi_parameter_attribution.png`
- `fig3_same_phi_parameter_attribution.pdf`
- `same_phi_pe10_vs_pe0p1_attribution.csv`

图上讨论：

这张图最适合用来反驳“是不是渗透率导致的”这个审稿人问题。以 Pe=0.1 为基准，只把 Pe=10 的某一类参数替换进去：

- 在 `phi≈0.60`，只换孔隙率几乎不变，只换曲折度也几乎不变，只换渗透率只把 `R_E` 降到约 0.78 倍。
- 但只换 H+/电化学参数，`R_E` 直接降到 `5.7e-5` 倍。
- 结构参数全换仍只有温和降低；结构 + 电化学全换才接近完整 Pe=10 的塌陷。

这张图给出的因果判断最硬：通道化结构是酸突破的输运原因，但震电响应崩塌的直接物理原因是电化学状态改变。

### 已生成 Figure D: 空间极性翻转图

文件：

- `fig4_spatial_polarity_reversal.png`
- `fig4_spatial_polarity_reversal.pdf`
- `spatial_peak_polarity_summary.csv`

图上讨论：

这张图把“幅值变化”提升为“源项极性重构”。Pe=0.1 和 Pe=1 的 normalized signed peak 在界面两侧保持相似极性结构；Pe=10 的曲线整体翻转，右侧柱状图也显示 R_E / T_E 的 signed peak 符号相对前两组改变。

这与 `zeta` 面板相互呼应：当 Pe=10 的 zeta 由负转正，有限偏移距 interface EM 的 signed response 也随之翻转。因此，interface EM response 不只是强弱指标，还可能记录 zeta-potential zero crossing。

### 已生成 Figure E: 同孔隙率溶蚀模式快照

文件：

- `fig5_dissolution_mode_snapshots_phi080.png`

图上讨论：

这张图用于把“Pe=0.1、1、10 是什么溶蚀模式”讲清楚。三列都取 `phi≈0.8` 附近：

- Pe=0.1: H+ 场仍更像从入口推进的前缘，流速集中仍受右侧残余孔隙网络限制。
- Pe=1: H+ 已沿主通道推进，速度场出现较宽的 channeling 通路。
- Pe=10: 高 H+ 条带和高速通道已经贯穿得更明显，表现为快速 acid breakthrough。

这张图最好与 Figure B 一起讲：Figure E 给出形态证据，Figure B 给出参数和震电证据。

### Figure 1: 机制示意图与时间尺度

画出：

`Pe -> dissolution regime -> phi/k/tau/H+ -> zeta, L(omega), sigma(omega) -> Schakel coefficients -> Liu finite-offset waveform`

突出 dissolution time 和 waveform time 不能混用。

### Figure 2: 三种 Pe 的溶蚀模式对比

建议放三组同一孔隙率附近的 H+ 场、矿物界面、流速场：

- Pe=0.1: phi≈0.79, t=7119.2 s。
- Pe=1: phi≈0.80, t=400.6 s。
- Pe=10: phi≈0.81, t=56.4 s。

旁边配 `k/k0`, `Hout`, `tortuosity`, ICN 的小表。

### Figure 3: 参数桥梁图

横轴用 porosity 或 dissolution time，纵轴画：

- `k/k0`
- `outlet_H_conc` 或 pH
- `zeta`
- `L_abs`
- `sigma_abs`
- `L_abs/sigma_abs`

关键是显示：Pe=10 在低孔隙率阶段就发生 H+ 突破和 `L/sigma` 塌陷。

### Figure 4: 同一孔隙率下的 one-at-a-time 归因

用 `phi≈0.60` 和 `phi≈0.80` 两个代表点，展示只换孔隙率、渗透率、曲折度、H+/电化学时 `RE_abs` 或 `Amax` 的比值。

这张图会非常有杀伤力，因为它直接证明“主控不是渗透率单独变化”。

### Figure 5: 界面转换系数与波形响应

画：

- `RE_abs`, `TTM_abs` 随 porosity / dissolution time。
- `Amax_waveform_spectral` 随 porosity / dissolution time。
- 标记 `phi_max_valid=0.95`，只把有效孔弹区间作为主结果。

### Figure 6: 空间极性与有限偏移距辐射

画三组 `waveform_spatial_peak_diagnostics`：

- signed peak 随 z。
- peak_abs 随 z。
- 标记 R_E/T_E 两侧峰值和 Liu dipole 预期峰值距离。

重点讲 Pe=10 的极性翻转，而不是只讲幅值变小。

## 结果写法草案

### Results subsection 1: Pe controls dissolution pathway rather than final state

Although all three simulations converged toward nearly identical final porosity and permeability, they followed distinct dissolution pathways. Low Pe produced a slow reaction-front-dominated evolution with delayed acid breakthrough, whereas Pe=1 and Pe=10 developed near-connected channelized pathways much earlier. This path dependence is critical because the seismoelectric model responds to the transient combination of hydraulic connectivity and pore-fluid chemistry, not to the final dissolved state alone.

### Results subsection 2: Acid breakthrough collapses electrokinetic conversion efficiency

The dominant divergence in the seismoelectric inputs occurs when channelized acid transport changes the electrochemical state of the pore fluid. At `phi≈0.60`, replacing only the Pe=0.1 permeability by the Pe=10 value reduces `RE_abs` to 0.78 of the baseline, whereas replacing only the H+ controlled electrochemistry reduces it to `5.7e-5`. This shows that permeability enhancement is not the primary cause of the response collapse; it is the acid-driven reduction of zeta potential and increase of dynamic conductivity.

### Results subsection 3: Polarity reversal records zeta-potential zero crossing

The signed finite-offset waveform provides a second, independent diagnostic of the same mechanism. Pe=0.1 and Pe=1 preserve the R_E-negative/T_E-positive polarity pattern, whereas Pe=10 reverses the dominant signed peaks. This reversal occurs when the H+ controlled zeta potential changes sign, indicating that the interface EM response records not only a weakening of coupling but also a reorganization of the electrokinetic source polarity.

### Discussion synthesis

These results reveal a competition between hydraulic amplification and electrochemical screening. Dissolution-induced channeling increases permeability and reduces tortuosity, which could enhance seismoelectric conversion through stronger relative fluid-solid motion. However, the same channels accelerate acid breakthrough, increase pore-fluid conductivity, and drive zeta potential toward zero or sign reversal. In the present simulations, the electrochemical screening term dominates, causing the high-Pe response to become much weaker despite stronger hydraulic connectivity.

## 推荐题目

1. Dissolution Pathways Control Seismoelectric Interface Responses Through Electrokinetic Coupling Collapse
2. Channelized Carbonate Dissolution Suppresses Seismoelectric Interface EM Responses by Acid-Driven Electrochemical Screening
3. Reactive-Transport Memory in Seismoelectric Interface Responses During Carbonate Dissolution
4. From Wormholing to Electrokinetic Short-Circuit: Why Faster Dissolution Weakens Interface EM Responses

## 参考文献支撑

### 溶蚀形态与通道化

- Fredd & Fogler (1999), *SPE Journal*, DOI: https://doi.org/10.2118/56995-PA
- Panga et al. (2005), *AIChE Journal*, DOI: https://doi.org/10.1002/aic.10574
- Szymczak & Ladd (2009), *JGR: Solid Earth*, DOI: https://doi.org/10.1029/2008JB006122
- Soulaine et al. (2017), *Journal of Fluid Mechanics*, DOI: https://doi.org/10.1017/jfm.2017.499
- Menke et al. (2023), *Scientific Reports*, DOI: https://doi.org/10.1038/s41598-023-37725-6

### 水力、电导和电动参数演化

- Noiriel et al. (2004), *Geophysical Research Letters*, DOI: https://doi.org/10.1029/2004GL021572
- Menke et al. (2015), *Environmental Science & Technology*, DOI: https://doi.org/10.1021/es505789f
- Pereira Nunes et al. (2016), *JGR: Solid Earth*, DOI: https://doi.org/10.1002/2015JB012117
- Jouniaux & Pozzi (1995), *Geophysical Research Letters*, DOI: https://doi.org/10.1029/94GL03307
- Soldi et al. (2024), *Geophysical Journal International*, DOI: https://doi.org/10.1093/gji/ggad457
- Cherubini et al. (2019), *JGR: Solid Earth*, DOI: https://doi.org/10.1029/2018JB017057
- Rembert et al. (2022), *Water*, DOI: https://doi.org/10.3390/w14101632

### 震电界面响应理论

- Pride (1994), *Physical Review B*, DOI: https://doi.org/10.1103/PhysRevB.50.15678
- Haartsen & Pride (1997), *JGR: Solid Earth*, DOI: https://doi.org/10.1029/97JB02936
- Garambois & Dietrich (2002), *JGR: Solid Earth*, DOI: https://doi.org/10.1029/2001JB000316
- Schakel & Smeulders (2010), *Journal of the Acoustical Society of America*, DOI: https://doi.org/10.1121/1.3263613
- Liu et al. (2018), *Journal of the Acoustical Society of America*, DOI: https://doi.org/10.1121/1.5020261
- Peng et al. (2019), *Geophysical Journal International*, DOI: https://doi.org/10.1093/gji/ggz249
- Hu et al. (2023), *JGR: Solid Earth*, DOI: https://doi.org/10.1029/2022JB025505

## 需要谨慎写的限制

1. `outlet_H_conc -> C_molL -> zeta/sigma` 是模型假设，不是实测完整离子组成。真实碳酸钙溶蚀会涉及 Ca2+, HCO3-, CO2/HCO3-/CO3 平衡和表面络合。
2. `phi > 0.95` 后主脚本标记 `valid_poroelastic=False`，论文主结论应优先使用有效孔弹区间。
3. 当前界面模型使用均质等效参数，没有直接把局部二维通道几何带入 Schakel 边界条件；通道几何影响是通过 `phi, k0, tau, H+` 等等效参数进入。
4. 波形幅值是当前模型尺度下的合成量，最稳妥的论文指标应包括归一化幅值、比值、极性和空间分布。

## 下一步最值得做的验证

1. 做跨 Pe 的 `same-porosity one-at-a-time` 图，重点展示 H+/电化学替换的主控作用。
2. 增加 `L_abs/sigma_abs` 随 porosity 的图，并与 `Amax` 做 log-log 关系。
3. 对 H+ 到流体电导率/zeta 的映射做敏感性：固定结构，只改变 `C_molL` 或 zeta 模型，检查极性翻转阈值是否稳定。
4. 如果有 Ca2+/HCO3- 或流体电导率输出，优先替换当前 HCl 近似模型。
5. 做 `valid_poroelastic=True` 区间内的主图，后期全溶空阶段放 SI 或只作为趋势参考。
