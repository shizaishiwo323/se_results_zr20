# 结果与讨论中文稿：Pe 控制的溶蚀模式与震电界面响应

本文稿将机制分析整理为中文论文式“结果与讨论”正文。正文嵌入 `paper_figure_sequence_analysis/pe_mode_mechanism/` 中生成的图件，并在关键机制链条处加入文献引用。

## 结果

### Pe 控制溶蚀路径，而不只是控制最终溶蚀状态

三组反应输运模拟最终趋向相近的溶蚀状态：孔隙率接近 1.0，渗透率约增加至初始值的 248 倍，曲折度接近 1.0。然而，三组达到这一终态的路径随 Pe 呈现显著差异。这种路径依赖性对震电响应解释十分关键，因为震电模型响应的是溶蚀过程中瞬时的水力连通性和孔隙流体化学状态，而不只是最终孔隙率或最终渗透率。

![图 1](../paper_figure_sequence_analysis/pe_mode_mechanism/fig5_dissolution_mode_snapshots_phi080.png)

**图 1.** 相近孔隙率下的溶蚀模式对比。三列分别为 Pe=0.1、Pe=1 和 Pe=10，孔隙率均接近 0.8。低 Pe 情况保留以反应前缘推进为主的 H+ 分布，而 Pe=1 和 Pe=10 形成更明显的通道化 H+ 输运和聚焦流速路径。

在孔隙率约为 0.8 的相近阶段，Pe=0.1 仍表现为较宽的反应前缘推进形态，而 Pe=1 和 Pe=10 已形成更加连通的高 H+ 与高流速路径（图 1）。这一现象与已有认识一致：传输与反应速率的相对关系会控制碳酸盐介质中的 compact dissolution、uniform dissolution 以及通道化或 wormhole-like dissolution（Fredd & Fogler, 1999; Panga et al., 2005; Szymczak & Ladd, 2009; Soulaine et al., 2017）。本研究的模拟进一步表明，这种形态差异并不只是水文结构上的结果；它还成为进入震电界面问题之前，控制孔隙流体电化学状态的上游因素。

### 通道化提前酸突破，并重塑电动参数桥梁

参数演化表明，仅用渗透率无法解释不同 Pe 下的震电响应差异。Pe=1 和 Pe=10 比 Pe=0.1 更早进入高渗状态，这与已有研究中“碳酸盐溶蚀可通过孔隙几何和连通性变化导致渗透率非线性增加”的认识一致（Noiriel et al., 2004; Menke et al., 2015; Pereira Nunes et al., 2016）。然而，三组之间最大的差异出现在电化学变量上。在高 Pe 情况下，H+ 浓度较早升高，zeta 电位趋近于零并转为正值，动态电导率增加，而 `|L(omega)|/|sigma(omega)|` 发生塌陷（图 2）。

![图 2](../paper_figure_sequence_analysis/pe_mode_mechanism/fig2_same_porosity_parameter_bridge.png)

**图 2.** 相同终态下的不同演化路径。各面板分别比较三组 Pe 条件下渗透率比值、H+ 浓度、zeta 电位、动态电导率、电动耦合效率 `|L(omega)|/|sigma(omega)|` 和波形峰值幅度随孔隙率的变化。

孔隙率约为 0.75 时的对比可以概括这一差异。Pe=0.1 时，`cH=1.0e-7 mol/L`，`zeta=-0.065 V`，`|sigma|=8.57e-3 S/m`，`|L|/|sigma|=1.78e-7`。在相近孔隙率下，Pe=10 的 `cH=1.74e-2 mol/L`，`zeta=1.60e-3 V`，`|sigma|=1.59e-1 S/m`，`|L|/|sigma|=1.30e-10`。相应地，波形峰值由 `7.81e15` 降至 `9.30e12`。因此，高 Pe 情况在水力上更连通，但在电动转换意义上效率更低。这一结果与电动理论相符：streaming-current coupling 受 zeta 电位、孔隙几何、渗透率和流体电导率共同控制（Jouniaux & Pozzi, 1995; Pride, 1994; Soldi et al., 2024）。它也与碳酸盐实验结果相一致，即反应改造可通过孔隙流体化学变化改变电导率和 streaming-potential coupling（Cherubini et al., 2019）。

### 参数替换表明电化学状态是转换塌陷的主导驱动因素

为区分结构效应和电化学效应，我们在固定孔隙率下进行了 one-at-a-time 参数替换：以 Pe=0.1 状态为基准，分别将部分变量替换为 Pe=10 的对应值。该检验直接评估响应差异究竟主要由孔隙率、渗透率、曲折度，还是由 H+ 控制的电化学状态所驱动。

![图 3](../paper_figure_sequence_analysis/pe_mode_mechanism/fig3_same_phi_parameter_attribution.png)

**图 3.** Pe=10 响应相对于 Pe=0.1 基准的 one-at-a-time 参数归因。在每个孔隙率下，分别替换结构变量和 H+ 控制的电化学变量。电化学替换导致 `R_E` 最显著塌陷。

孔隙率约为 0.60 时，参数归因结果最为清晰。仅替换孔隙率时，`R_E` 变化小于 1%；仅替换曲折度时，`R_E` 只出现小幅增加。将渗透率替换为 Pe=10 的值时，`R_E` 降至 Pe=0.1 基准的 0.78 倍；而仅替换 H+ 控制的电化学状态时，`R_E` 降至基准的 `5.7e-5`。若只替换所有结构变量而不替换电化学状态，`R_E` 仍保持在基准的 0.80 倍；当结构变量和电化学变量同时替换时，`R_E` 降至 `4.5e-5`。这一对比表明，通道化是将酸提前输运到下游的传输机制，但震电转换塌陷的主导模型驱动因素是酸突破产生的电化学状态。

这一差异对解释溶蚀驱动的震电信号十分重要。渗透率增加可以增强流固相对运动并改变动态渗透率，但同一组通道化路径也会提高离子强度，并削弱有效电动源项。在本组模拟中，电化学抑制作用强于水力放大作用。因此，用 `|L(omega)|/|sigma(omega)|` 表征控制机制，比单独使用孔隙率或渗透率更接近震电响应的物理本质。

### 高 Pe 情况下的极性翻转与 zeta 电位变号相关

有限偏移距空间响应为上述机制提供了独立检验。Pe=0.1 和 Pe=1 保持相近的反射侧与透射侧 signed peak 分布，而 Pe=10 的界面两侧主导 signed peak 发生翻转（图 4）。这一翻转出现在模型中的 zeta 电位发生变号之后。该结果表明，interface EM response 不仅记录电动转换强度的减弱，也记录电动源项极性的改变。

![图 4](../paper_figure_sequence_analysis/pe_mode_mechanism/fig4_spatial_polarity_reversal.png)

**图 4.** 有限偏移距极性翻转。Pe=0.1 和 Pe=1 保持相似的 signed spatial pattern，而 Pe=10 表现出主导极性的翻转。该翻转与 H+ 控制的电化学模型所推断的 zeta 电位变号一致。

这一行为与 Pride 和 Schakel 理论框架一致：在该类模型中，电动耦合通过动态耦合系数和电导率进入电磁-孔弹耦合边界值问题（Pride, 1994; Schakel & Smeulders, 2010）。它也与有限偏移距 VSEP 理论相符，即 interface EM 波形及其空间极性受转换系数的符号、幅值和接收几何共同控制（Liu et al., 2018）。因此，图 4 中的符号翻转支持这样一种解释：Pe=10 并非仅仅通过几何变化削弱响应，而是改变了有效电动源项的极性。

## 讨论

### 水力放大与电化学屏蔽之间的竞争

上述结果支持一个双分支机制，用于解释溶蚀控制的震电响应（图 5）。在水力分支中，溶蚀增加渗透率并降低曲折度，从而可能增强流固相对运动并改变动态渗透率。在电化学分支中，通道化酸输运提高 H+ 浓度，使 zeta 电位趋近于零或发生变号，并增加动态电导率。前一分支可能增强震电转换；后一分支则通过减少可用电动耦合并增强电导屏蔽来抑制转换。

![图 5](../paper_figure_sequence_analysis/pe_mode_mechanism/fig1_mechanism_cartoon.png)

**图 5.** 溶蚀模式与震电界面响应之间的概念机制。通道化增强水力连通性，但酸突破通过改变 zeta 电位和电导率导致电动转换塌陷。

在本组模拟中，电化学分支占主导。这解释了一个表面上反直觉的结果：高 Pe 情况使介质在水力上更连通，但在震电响应上反而更弱。该效应并不是 Pe 增大本身的普遍后果，而是因为 Pe=10 的溶蚀路径造成了早期酸突破，使介质由电动活跃状态转向高电导、弱耦合状态。

### 与已有震电和反应输运研究的关系

已有震电研究表明，界面产生的 EM 响应受孔隙率、渗透率、流体盐度或电导率、黏度和电动耦合控制（Pride, 1994; Garambois & Dietrich, 2002; Schakel & Smeulders, 2010; Liu et al., 2018）。反应输运和碳酸盐溶蚀研究则分别表明，溶蚀模式会控制孔隙几何、渗透率、曲折度和溶质突破行为（Noiriel et al., 2004; Szymczak & Ladd, 2009; Menke et al., 2015; Soulaine et al., 2017; Menke et al., 2023）。本研究将这两条研究线索连接起来，展示了溶蚀路径如何在 interface EM 波形生成之前，先改变动态电动参数桥梁。

这一联系将震电监测的解释从静态物性对比问题推进到具有路径依赖的反应输运问题。在酸突破延迟的条件下，即使孔隙率增加，介质仍可保持较强电动活性；而在快速通道化突破条件下，同样的水力连通性增强可能抑制界面响应。这种路径依赖性说明，在反应改造过程中，仅用孔隙率或渗透率预测震电幅值是不充分的。

### 启示与模型限制

本组模拟中最有用的诊断量是 `|L(omega)|/|sigma(omega)|`。该比值同时包含可用于转换的动态电动耦合和促进电导屏蔽的动态电导率。与孔隙率或渗透率相比，它在 Pe 对比中更直接地追踪波形峰值，并提供了一个连接反应输运输出和震电观测量的物理桥梁。

几个限制需要在论文中明确说明。首先，从出口 H+ 浓度到电解质浓度、电导率和 zeta 电位的转换是简化的电化学映射。完整碳酸盐体系还应包含 Ca2+、HCO3-、CO3 形态分布、CO2 平衡和表面络合。其次，主要解释应集中在 `valid_poroelastic=True` 的区间，因为当孔隙率接近完全溶蚀极限时，孔弹框架的适用性会变弱。第三，当前界面模型使用等效介质参数，而没有在 Schakel 边界条件中显式解析局部二维通道几何。这些限制不消除模型结果的内部一致性，但限定了后续需要加强的检验方向。

### 主要认识

核心结果是：溶蚀可以提高水力连通性，同时降低震电活性。在本组模拟中，高 Pe 通道化加速酸突破，使 zeta 电位趋近于零或发生变号，并提高动态电导率。由此造成的 `|L(omega)|/|sigma(omega)|` 塌陷抑制了 interface EM response，并可导致有限偏移距极性翻转。因此，震电响应携带的是反应输运路径信息，而不只是最终孔隙率或渗透率信息。

## 正文引用文献

- Cherubini, A., Garcia, B., Cerepi, A., & Revil, A. (2019). Influence of CO2 on the electrical conductivity and streaming potential of carbonate rocks. *Journal of Geophysical Research: Solid Earth*. https://doi.org/10.1029/2018JB017057
- Fredd, C. N., & Fogler, H. S. (1999). Optimum conditions for wormhole formation in carbonate porous media: Influence of transport and reaction. *SPE Journal*. https://doi.org/10.2118/56995-PA
- Garambois, S., & Dietrich, M. (2002). Full waveform numerical simulations of seismoelectromagnetic wave conversions in fluid-saturated stratified porous media. *Journal of Geophysical Research: Solid Earth*. https://doi.org/10.1029/2001JB000316
- Jouniaux, L., & Pozzi, J. P. (1995). Permeability dependence of streaming potential in rocks for various fluid conductivities. *Geophysical Research Letters*. https://doi.org/10.1029/94GL03307
- Liu, Y., Smeulders, D., Su, Y., & Tang, X. (2018). Seismoelectric interface electromagnetic wave characteristics for the finite offset Vertical Seismoelectric Profiling configuration: Theoretical modeling and experiment verification. *Journal of the Acoustical Society of America*. https://doi.org/10.1121/1.5020261
- Menke, H. P., Bijeljic, B., Andrew, M. G., & Blunt, M. J. (2015). Dynamic three-dimensional pore-scale imaging of reaction in a carbonate at reservoir conditions. *Environmental Science & Technology*. https://doi.org/10.1021/es505789f
- Menke, H. P., Maes, J., & Geiger, S. (2023). Channeling is a distinct class of dissolution in complex porous media. *Scientific Reports*. https://doi.org/10.1038/s41598-023-37725-6
- Noiriel, C., Gouze, P., & Bernard, D. (2004). Investigation of porosity and permeability effects from microstructure changes during limestone dissolution. *Geophysical Research Letters*. https://doi.org/10.1029/2004GL021572
- Panga, M. K. R., Ziauddin, M., & Balakotaiah, V. (2005). Two-scale continuum model for simulation of wormholes in carbonate acidization. *AIChE Journal*. https://doi.org/10.1002/aic.10574
- Pereira Nunes, J. P., Blunt, M. J., & Bijeljic, B. (2016). Pore-scale simulation of carbonate dissolution in micro-CT images. *Journal of Geophysical Research: Solid Earth*. https://doi.org/10.1002/2015JB012117
- Pride, S. R. (1994). Governing equations for the coupled electromagnetics and acoustics of porous media. *Physical Review B*. https://doi.org/10.1103/PhysRevB.50.15678
- Schakel, M. D., & Smeulders, D. M. J. (2010). Seismoelectric reflection and transmission at a fluid/porous-medium interface. *Journal of the Acoustical Society of America*. https://doi.org/10.1121/1.3263613
- Soldi, M., Guarracino, L., & Jougnot, D. (2024). Predicting streaming potential in reactive media: The role of pore geometry during dissolution and precipitation. *Geophysical Journal International*. https://doi.org/10.1093/gji/ggad457
- Soulaine, C., Roman, S., Kovscek, A., & Tchelepi, H. A. (2017). Mineral dissolution and wormholing from pore-scale simulations. *Journal of Fluid Mechanics*. https://doi.org/10.1017/jfm.2017.499
- Szymczak, P., & Ladd, A. J. C. (2009). Wormhole formation in dissolving fractures. *Journal of Geophysical Research: Solid Earth*. https://doi.org/10.1029/2008JB006122
