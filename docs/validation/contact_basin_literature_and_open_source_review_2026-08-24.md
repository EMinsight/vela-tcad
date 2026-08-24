# Contact basin 文献与开源 TCAD 实现调研

日期：2026-08-24

## 结论先行

1. 本次针对 `contact basin`、`contact_basin`、`contact-partitioned quasi-Fermi`和 `local-reference quasi-Fermi` 进行了精确词和组合词检索。在检索到的半导体论文以及 DEVSIM、Charon、Genius 公开代码中，未发现与 Vela 同名的 `contact_basin` 算法或配置项。这是检索范围内的否定结果，不是对所有商业 TCAD 内部实现的绝对证明。
2. Vela `contact_basin` 是一种数值坐标/规范选择，而不是新物理模型。它组合了已有的准费米势变量、对数/Slotboom 变量改善条件数、全局无量纲化，以及多源最短路/Voronoi 分区思想。
3. 开源 TCAD 的可核对代码更多采用“绝对电势+载流子密度”或“全局缩放的准费米势”，并通过正值更新、Kahan 求和、扩展精度或无量纲化降低数值误差；未见按最近接触电压逐节点切换准费米参考的公开实现。

## Vela 实现契约

Vela 对电子和空穴分别建立参考场。以电子为例，节点 `i` 的绝对准费米势表示为：

`Phi_n(i) = R_n(i) + delta_Phi_n(i)`

其中 `R_n(i)` 是节点所属接触 basin 的接触电压，`delta_Phi_n(i)` 是 Newton 求解的小增量。跨 basin 边 `i-j` 在 `i` 参考系中重建：

`delta_Phi_n(j | i) = R_n(j) + delta_Phi_n(j) - R_n(i)`

因而载流子密度、Scharfetter-Gummel 通量和接触电流仍取决于相同的绝对准费米势，物理量对参考分区应保持不变。局部表示的目的是防止小于全局偏压 ULP 的 Newton 更新在浮点加法中丢失。

分区规则：

- 从非金属栅接触中选取种子；按接触节点平均净掺杂的符号，分别作为电子或空穴多数载流子种子。
- 只在输运材料的边上建立邻接图，边权为网格边长。
- 使用确定性的多源 Dijkstra 传播最近接触所有权；等距时用接触顺序破平局。
- 与输运图不连通的节点使用到种子节点的欧氏距离回退。

相关本地实现：`src/solver/NewtonSolver.cpp`、`src/equation/CoupledDDAssembler.cpp`、`tests/test_newton_solver.cpp`。

## 一手证据对照

| 证据 | 公开内容 | 对 Vela 的支持程度 |
|---|---|---|
| 对数/准费米变量论文 | Bonilla et al. 比较传统 DD 和无量纲对数变量，给出准费米势与 `log(n/ni)`、`log(p/ni)` 的关系，并报告对数变量在粗网格下的正值性、无振荡和非线性稳健性。 | 强力支持“改变数值变量可改善稳健性”，但不是 contact-basin 分区。 |
| 欧姆接触边界层分析 | Farrell/Peschka 展示准费米势在欧姆接触附近可出现尖锐边界层，并比较有限差分和 SG 型有限体积离散。 | 支持“接触附近的准费米数值表示需单独关注”，但未提出 Vela 的分区参考。 |
| DEVSIM 开源代码 | `simple_physics.py` 公开示例以 `Potential`、`Electrons`、`Holes` 为未知量，连续性方程使用 `positive` 更新；电荷组合使用 `kahan3`，项目本身支持扩展浮点精度。 | 佐证其它工程路线：正值更新、补偿求和和扩展精度；公开示例不是 contact basin。 |
| Charon 开源代码 | `Charon_Scaling_Parameters.cpp` 以 `V0=kBT/q`、浓度尺度、扩散尺度等构造全局电压、场、电流和时间缩放。 | 强力支持全局无量纲化，但未见逐接触分区局部参考。 |
| Genius-TCAD-Open | 官方仓库明确提供 2-D DD、DC/瞬态/AC/混合模式，源码树包含 `ddm1`、`ddm2`等求解器；公开 MOSCAP 示例使用逐段 DCSWEEP 和 Potential damping。 | 佐证主流做法为偏压延续+全局阻尼；未检索到 contact-basin 同名或等价公开选项。 |
| 多源最短路/Voronoi 文献 | Mitchell-Papadimitriou 的加权区域最短路工作明确说明算法可推广到多源以计算 Voronoi 分区。 | 直接支持 Vela basin 所用的图分区算法骨架，但其应用于准费米参考是 Vela 的工程组合。 |

## 一手来源

- Bonilla et al., *A comparison of formulations and non-linear solvers for computational modelling of semiconductor devices*: https://pmc.ncbi.nlm.nih.gov/articles/PMC12045845/
- Farrell, *Drift-diffusion models for innovative semiconductor devices and their numerical solution*: https://refubium.fu-berlin.de/bitstream/handle/fub188/37941/Habil_Farrell.pdf?isAllowed=y&sequence=3
- DEVSIM 官方仓库: https://github.com/devsim/devsim
- DEVSIM `simple_physics.py`: https://github.com/devsim/devsim/blob/main/python_packages/simple_physics.py
- Charon 官方仓库: https://github.com/tcadsoftware/charon
- Charon 缩放实现: https://github.com/tcadsoftware/charon/blob/main/src/Charon_Scaling_Parameters.cpp
- Genius-TCAD-Open 官方仓库: https://github.com/cogenda/Genius-TCAD-Open
- Genius MOSCAP 公开示例: https://github.com/cogenda/Genius-TCAD-Open/blob/master/examples/MOSCAP/cap.inp
- Mitchell and Papadimitriou, *The Weighted Region Problem: Finding Shortest Paths Through a Weighted Planar Subdivision*: https://doi.org/10.1145/102782.102784

## 对后续开发的建议

1. 文档中将全称固定为 **contact-partitioned quasi-Fermi gauge**，配置键保留 `contact_basin`，避免被误解为 Sentaurus 原生模型名。
2. 回归必须同时验证：物理状态参考不变性、跨 basin SG 通量、接触电流、四端 KCL，以及小于全局偏压 ULP 的 Newton 更新保留。
3. basin 边界应只改变数值坐标，不得在密度、重组或通量中引入不连续物理量；跨分区边的局部参考转换应作为硬性单元测试。
4. 保留“全局缩放+补偿求和+严格连续性/KCL 门槛”；`contact_basin` 只解决准费米坐标的浮点消减，不替代残差收敛和守恒验收。
