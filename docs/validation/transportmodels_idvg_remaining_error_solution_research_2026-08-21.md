# TransportModels DD/DG Id–Vg 剩余误差：跨工具实现与可行解决方案

日期：2026-08-21
范围：仅针对 TransportModels MOS 的 DD/DG 对比验证，不扩展到 AC、瞬态、自热或 MixedMode。

## 1. 结论先行

当前剩余问题应拆成两条线，并按以下顺序处理：

1. **先处理 DD/DG 共享的约 15 mV 阈值偏移。** 现有数据表明，Vela 的 DD 与 DG 在阈值附近相对 Sentaurus 几乎具有相同的负向栅压偏移；这不支持“剩余误差主要由 DG 方程造成”的单一解释。第一嫌疑是 `Fermi + EffectiveIntrinsicDensity(OldSlotboom)` 的统计/BGN 语义，其次是偏压相关的 Enormal、迁移率与高场驱动力。
2. **再补齐 Sentaurus 默认 DG 连续性方程语义。** 原始 Sentaurus 算例没有 `DirectQuantumCorrection`，而 Vela 当前把电子量子势直接并入载流子统计及 SG 漂移势；两者不是同一个模型开关。
3. **DG 界面条件和离散格式后移。** 阶段 7 在 `Vg=1 V, Vd=2 V` 已达到量子势表面 p95 误差 `14.7858 mV`、电子浓度表面 p95 误差 `0.08288 dex`，且 Id–Vd 最大误差为 `1.0809%`。这说明已选 `sentaurus_box + neutral_continuous` 不是当前首要矛盾。

建议先做低成本的 Sentaurus 2×2 语义试验和 Vela BGN 开关试验，再决定是否修改连续性方程。不要先继续拟合 DG 的 `gamma` 或有效质量。

## 2. 现有残差的重新分解

### 2.1 等电流阈值偏移

对 DD 冻结基线与阶段 7 DG 曲线分别在对数电流上插值：

| 电流准则 (A/µm) | DD：Vela−Sentaurus | DG：Vela−Sentaurus | 解释 |
|---:|---:|---:|---|
| `1e-10` | `-14.73 mV` | `-14.48 mV` | 几乎相同 |
| `1e-8` | `-15.42 mV` | `-14.34 mV` | 几乎相同 |
| `1e-6` | `-15.09 mV` | `-13.09 mV` | DG 只减少约 2 mV |
| `1e-5` | `-19.46 mV` | `-14.04 mV` | 开始叠加迁移率/强反型效应 |
| `1e-4` | `-17.84 mV` | `-10.85 mV` | 进入导通区 |

在 `Vg=-0.20 V`，DD 与 DG 的 Vela/Sentaurus 电流比分别为 `1.5695` 和 `1.5068`；在 `Vg=-0.04 V` 分别为 `1.6131` 和 `1.4956`。DG 相对 DD 的额外对数残差只有 `-0.0177 dex` 和 `-0.0328 dex`，明显小于两者共同的约 `0.18–0.21 dex` 残差。

在 `Vg=2.2 V`，DD 与 DG 的 Vela/Sentaurus 电流比分别为 `0.97904` 和 `0.97877`，同样高度一致。这进一步说明高栅压端仍有共享的经典输运/迁移率因素。

上述计算直接使用冻结 DD 曲线和阶段 7 DG 曲线。需要注意：DD 冻结曲线早于阶段 7 的部分 DG 校正；正式定论前应使用阶段 7 的共同材料与求解设置重跑 DD。这是后续任务 0 的组成部分。

### 2.2 与既有 Frozen-Q 结论的关系

早期 `Vg=1 V, Vd=2 V` Frozen-Q 试验把自洽 DG 端点电流误差从 `9.8424%` 降到 `1.3156%`，证明**当时的端点误差**主要来自 DG 场。阶段 2–6 已通过材料驱动、参数、界面和离散格式工作消除了大部分端点误差；阶段 7 的新残差集中在 Id–Vg 阈值附近，不能沿用早期端点分类。

## 3. Sentaurus T-2022.03 手册与原始算例核对

### 3.1 DG 方程与输运耦合

Sentaurus Device User Guide T-2022.03 第 14 章给出以下关键行为：

- 默认采用 potential-like DG 方程；`Formula=0` 为默认形式，`Formula=1` 把质量置于微分算子内部。Fermi 统计下，potential-like 与 density-based 形式即使参数取默认值也不严格等价。
- 若 DG 在内界面两侧都求解，量子势变量连续并满足通量条件；若绝缘体侧不求解，应使用由一维势垒/WKB 推导的非齐次 Neumann 条件。手册明确建议在量子效应改变沟道载流子浓度的界面使用该势垒条件，齐次 Neumann 更适用于远离沟道的界面。
- **默认连续性方程使用指数近似**，其解可理解为经典载流子浓度乘以量子因子。只有在 `Math` 中设置 `DirectQuantumCorrection`，量子势才作为量子力学带边直接加到静电势上。
- 默认模式同时存在 classical density 与 quantum-mechanical density；在 PMI/Tcl 公式中必须明确区分 `eDensity` 与 `eQMDensity`。
- 手册明确警告：DG 改变载流子分布和电场，因此基于经典模型或 Van Dort 校准的迁移率/复合模型可能需要重新校准。

原始 [`IdVgs_des.cmd`](../../build-release/reference_tcad/transportmodels_sentaurus2022/run02/full_raw/IdVgs_des.cmd) 只启用了 `eQuantumPotential`，`Math` 中没有 `DirectQuantumCorrection`。因此 Sentaurus 参考曲线使用默认指数型耦合。

Vela 当前在 [`CoupledDDAssembler.cpp`](../../src/equation/CoupledDDAssembler.cpp) 中使用

\[
\psi_{n,\mathrm{eff}}=\psi-Q_n,
\]

并同时用 `ψ−Qn` 计算 Fermi 载流子浓度和电子 SG 漂移势。这在语义上更接近 Sentaurus 的 `DirectQuantumCorrection`，也是 DEVSIM 示例采用的做法。该差异只影响 DG，不能单独解释 DD/DG 共有的约 15 mV 偏移，但可能解释 DG 在不同栅压下相对 DD 的小幅残差变化。

### 3.2 Fermi + OldSlotboom

原始 Sentaurus 算例同时使用 `Fermi` 和 `EffectiveIntrinsicDensity(OldSlotboom)`。手册指出，OldSlotboom 参数通常按 Maxwell–Boltzmann 统计提取；Sentaurus 在 Fermi 统计下默认加入表观 BGN 修正，`NoFermi` 才会关闭。

Vela 已实现对应修正，但 [`BandgapNarrowing.h`](../../include/vela/physics/BandgapNarrowing.h) 为保持旧算例兼容，将 `fermiStatisticsCorrection` 默认设为 `false`；阶段 7 的字符串配置 `"bandgap_narrowing": "old_slotboom"` 没有显式打开该开关。这是目前最直接、成本最低、且能同时影响 DD/DG 阈值的语义缺口。

### 3.3 Enormal 与高场驱动力

Sentaurus 手册第 15 章规定：

- `Enormal` 默认取电场与最近半导体–绝缘体界面几何法向的投影；界面顶点可用 `NormalFieldCorrection` 修正同顶点电荷屏蔽造成的低估。
- Enhanced Lombardi 把体迁移率、声子散射和表面粗糙散射按 Matthiessen 规则组合，并按距界面的指数因子关闭远场表面项。
- `GradQuasiFermi` 是漂移扩散高场饱和的默认驱动力；接触相邻单元默认改用电场。低密度时还可在准费米梯度与电场之间插值以改善数值稳定性。
- 表面迁移率对界面法向网格非常敏感，手册建议沟道 Si/SiO₂ 下前两层约 `0.1 nm`，Lombardi 扩展项甚至约 `0.05 nm`。

Vela 当前按最近界面段构造单元法向，并以三角形内 `∇ψ·n` 生成单元 Enormal，物理方向正确；但它是**单元值**，没有 Sentaurus 的界面顶点修正，也没有明确复制“接触单元退化为电场”的高场驱动力规则。已有 Frozen-Q 空间审计显示 Vela 沟道中位 Enormal 为 `250642 V/cm`，Sentaurus 为 `159463 V/cm`，而迁移率分别为 `18.94` 和 `44.00 cm²/V·s`。端点电流仍接近，是 Enormal、迁移率和纵向驱动力互相补偿的结果，不能视为模型已对齐。

## 4. 开源 TCAD 实现对照

| 工具 | 与本问题相关的实现 | 可借鉴点 | 不能直接照搬之处 |
|---|---|---|---|
| Sentaurus | potential-like DG；默认指数型连续性方程；可选 direct band-edge；几何 Enormal；Lombardi/IALMob；Fermi-BGN 修正 | 目标语义与验收基准 | 闭源，部分离散细节只能用 A/B 算例反演 |
| Vela | potential-based DG；`ψ−Qn` 同时进入 Fermi 浓度和 SG 漂移势；单元 `∇ψ·n`；Enhanced Lombardi；准费米梯度高场 | 已有 Frozen-Q、固定状态残差和 42 点回归基础 | 默认 DG 输运耦合及 Fermi-BGN 开关与参考 deck 不同 |
| [DEVSIM density-gradient](https://github.com/devsim/devsim_density_gradient) | `EN=EC+Le`，电流 SG 使用 `EN`；代码注明准费米实现假设 Boltzmann；含氧化层/WKB 表面项；用 gamma ramp 求解 | 适合验证 direct-band-edge 分支、WKB 边界和 gamma continuation | README 明确说完整 DD 主要只在 MOSCAP 耦合下测试；不是 Sentaurus 默认指数模式的 oracle |
| [Genius-TCAD-Open](https://github.com/cogenda/Genius-TCAD-Open) | 完整 DG FVM；Si/oxide 未求解侧使用 WKB 非齐次边界；可切换高场迁移率自洽；界面附近 Enormal 使用半导体场和绝缘体位移关系重构 | 可直接借鉴 WKB 边界、偏压内自洽迁移率开关、界面场诊断 | Enormal 公式不等于 Sentaurus 默认几何投影，应作为备选模式而非默认替换 |
| [Charon](https://github.com/tcadsoftware/charon) | 公共源码中可识别到 `QuantumPotentialFlux`/`FieldMag` 对 `∇φ + fit·∇Q` 的组合；MOSFET mobility 把 bulk 与 perpendicular mobility 按 Matthiessen 组合；Shirahata 同时提供 edge/IP 实现并显式使用界面法向与距离 | 模块化“体迁移率/垂直场迁移率/边值离散”结构；同一模型同时输出 edge 与 IP 诊断 | 在公共源码文件清单中未找到完整 density-gradient PDE 求解器，不能把 Charon 当作 DG 方程 oracle |

### 4.1 Genius 的两个具体实现线索

- [`dg_semiconductor.cc`](https://github.com/cogenda/Genius-TCAD-Open/blob/master/src/solver/dg/dg_semiconductor.cc) 在绝缘体界面附近构造
  `E_eff,n = ζ(E·n) + η[(εox/εsi)(Eox·n) − E·n]`，并把它传给表面迁移率。这适合做 Vela 的实验性 `displacement_reconstructed` Enormal 模式。
- [`dg_boundary_is_interface.cc`](https://github.com/cogenda/Genius-TCAD-Open/blob/master/src/solver/dg/dg_boundary_is_interface.cc) 根据势垒高度、氧化层有效质量和穿透深度生成量子势边界通量。这与 Sentaurus 非齐次势垒边界及 DEVSIM/Garcia-Loureiro 实现方向一致。

### 4.2 DEVSIM 的两个具体实现线索

- [`dg_physics.py`](https://github.com/devsim/devsim_density_gradient/blob/main/dg_physics.py) 直接把量子势加入电子带边，并用该带边生成 SG 电流；它可作为 Vela 当前 direct-band-edge 分支的独立对照。
- [`test_2d.py`](https://github.com/devsim/devsim_density_gradient/blob/main/test_2d.py) 先把氧化层 gamma 从 `0.1` ramp 到 `0.25`，再把 Si gamma ramp 到 `3.6`。这种参数 continuation 适合提高新耦合分支的收敛性，但不改变最终物理模型。

### 4.3 Charon 的两个具体实现线索

- [`Charon_Mobility_MOSFET_impl.hpp`](https://github.com/tcadsoftware/charon/blob/main/src/evaluators/Charon_Mobility_MOSFET_impl.hpp) 对 edge 和 integration point 分别计算，并将 bulk/perpendicular mobility 的倒数相加。这提示 Vela 应将同一偏压下的“单元 Enormal、边迁移率、边电流”成套导出，而不只比较节点场。
- [`Charon_Mobility_Shirahata_impl.hpp`](https://github.com/tcadsoftware/charon/blob/main/src/evaluators/Charon_Mobility_Shirahata_impl.hpp) 显式记录氧化层界面几何、法向、距离和边方向投影，可作为 Vela Enormal 几何审计的结构参考。

## 5. 文献证据

- Wettstein 等的 [general-purpose DG implementation](https://doi.org/10.1080/1065514021000012363) 将 DG 与 Schrödinger 结果对照，并专门讨论绝缘体中的模型修改；这支持“先校准电荷/势，再校准电流”的顺序。
- Garcia-Loureiro 等的 [3-D multigate DG implementation](https://doi.org/10.1109/TCAD.2011.2107990) 是 DEVSIM 示例引用的界面/离散依据，支持对非齐次界面项进行独立单元测试。
- [Wettstein 等 2001 年 unstructured-grid DG 工作](https://doi.org/10.1109/16.902727) 指出简化 DG 形式可能匹配端口特性却给出错误的器件内部密度；这正是为什么不能用单个 Id 端点拟合量子势或迁移率。
- [Riddet 等关于 DG 边界条件的研究](https://doi.org/10.1007/s10825-008-0222-6) 表明接触边界条件也会影响量子修正结果，支持把接触邻近的高场/量子边界规则纳入审计。
- [Carapezzi 等的 DG+mobility calibration](https://doi.org/10.1016/j.sse.2020.107902) 使用 Schrödinger–Poisson/微观迁移率结果分别校准 DG 与迁移率，支持避免以迁移率参数补偿量子电荷错误。

## 6. 推荐实施顺序

| 优先级 | 任务 | 具体方法 | 预计工作量 | 通过标准 | 主要风险 |
|---:|---|---|---:|---|---|
| P0 | 重跑共享 DD 控制 | 用阶段 7 完全相同的材料、Fermi/BGN、迁移率和离散设置重跑 21 点 DD Id–Vg | 0.5 天 | 确认约 15 mV 共同偏移不是旧 DD 配置伪差 | 旧/新配置不可直接比较 |
| P0 | Sentaurus 2×2 语义 oracle | 原始 DG deck 分别组合 `{默认, DirectQuantumCorrection}` × `{默认 Fermi-BGN, NoFermi}`，先跑 `Vg=-0.20,-0.04,0.12,0.28,1.0 V` | 0.5–1 天 | 定量分解 DG 连续性语义与 BGN 语义各自贡献 | 需要 VM 中重新跑 SDevice |
| P1 | Vela Fermi-BGN A/B | 在 DD 与 DG 上显式设置 `fermi_statistics_correction=true/false`；不改其他参数 | 0.5–1 天 | 共同阈值偏移明显缩小，且 Id–Vd 不退化 | 修正方向可能暴露材料参数的二次不一致 |
| P1 | 偏压分辨的空间 oracle | 在 5 个关键 Vg 导出 Sentaurus/Vela 的 `Qn,n,φn,Enormal,μn,|∇φn|`，沿源端/沟道中部/漏端三条法向剖面比较 | 1–2 天 | 能把电流误差分为电荷、迁移率和驱动力三项 | 数据映射需保持同一几何位置而非最近点误配 |
| P2 | 实现 DG 输运耦合开关 | 保留 `direct_band_edge`；新增 `sentaurus_exponential`，显式维护 classical/QM density，避免量子势在浓度和 SG 漂移势中重复计入 | 2–4 天 | Frozen-Q 单点及 5 点 Id–Vg 与 Sentaurus 默认 deck 同向改善；Direct 分支不回归 | Fermi 统计下量子因子的精确定义和 Jacobian 较复杂 |
| P2 | Enormal/迁移率离散对齐 | 新增 vertex/edge Enormal 导出；复刻 Sentaurus 几何法向与接触规则；实现可选 `NormalFieldCorrection`；Genius 位移重构仅作实验分支 | 2–3 天 | 5 个 Vg 上 Enormal、μn 和驱动力不再依赖互相补偿；曲线误差改善 | 过度拟合单一 Vg 或单一路径 |
| P3 | WKB 界面条件分支 | 当不在氧化层求解 DG 时，按 Sentaurus/Genius/DEVSIM 增加势垒穿透边界；与当前两侧连续求解做受控 A/B | 2–4 天 | Qn/n 空间误差在多偏压下降，且不靠 gamma 再拟合 | 当前两侧求解已较好，收益可能有限 |
| P4 | DG 方程高级形式 | 再评估 density-based、Formula1、各向异性质量/AutoOrientation | 4–7 天 | Schrödinger–Poisson 或 Sentaurus 多偏压电荷剖面持续改善 | 参数维度增加、可辨识性下降 |
| P4 | 网格收敛 | 对沟道界面法向 0.2/0.1/0.05 nm 做网格序列；分别观察 Qn、Enormal、μn 和 Id | 1–2 天 | 0.1→0.05 nm 的主要量变化低于预设门槛 | 重网格会引入节点映射差异 |

## 7. P0/P1 的最小实验矩阵

建议先跑以下 8 个模型变体、共 40 个关键偏置点（各变体可在一条扫压路径内复用前态），而不是立即重跑全部曲线：

| 平台 | 分支 | 受控变量 | Vg 点 |
|---|---|---|---|
| Vela | DD | Fermi-BGN off/on | `-0.20,-0.04,0.12,0.28,1.0` |
| Vela | DG | Fermi-BGN off/on | 同上 |
| Sentaurus | DG | default/DirectQC × default/NoFermi | 同上 |

每个点至少保存：`Id`、沟道积分电子电荷、表面势、`eQuantumPotential`、`eQMDensity`、`eEnormal`、`eMobility`、`eGradQuasiFermi`。若 P0 结果支持某一根因，再扩展到完整 21 点 Id–Vg 和 21 点 Id–Vd。

## 8. 验收门槛

最终仍沿用阶段 7 的主门槛，并增加防止“端点拟合”的约束：

| 指标 | 门槛 |
|---|---:|
| DG Id–Vg transition 最大绝对对数误差 | `≤0.15 dex` |
| DG Id–Vg on 最大相对误差 | `≤10%` |
| DG Id–Vd 最大相对误差 | 不高于当前 `1.0809%` 的显著回归；建议上限 `2%` |
| 等电流阈值偏移 (`1e-10` 至 `1e-6 A/µm`) | 建议 `|ΔVg|≤5 mV` |
| 表面 Qn p95 绝对误差 | `≤20 mV` |
| 表面 n p95 绝对对数误差 | `≤0.2 dex` |
| 偏压空间一致性 | 5 个 Vg 上电荷、Enormal、迁移率、驱动力均需报告，禁止只以 `Vg=1 V` 端点判定 |

## 9. 最终建议

下一步不应直接进入大规模 DG 方程重构。最优路径是：

1. 重跑同配置 DD，确认共有残差；
2. 完成 Sentaurus `DirectQC × NoFermi` 2×2 oracle；
3. 运行 Vela 已具备的 Fermi-BGN 开关 A/B；
4. 根据结果决定先实现 `sentaurus_exponential`，还是先修 Enormal/迁移率离散；
5. 只有在多偏压空间量仍指向界面时，才进入 WKB 边界和高级 DG 形式。

这一顺序的核心是先消除可验证的语义差异，再做参数校准；否则很容易用 Lombardi 系数或 DG 有效质量抵消错误的统计/连续性方程语义，并在另一段栅压上重新产生误差。

## 10. 主要本地证据

- [阶段 7 完整回归](transportmodels_dg_phase7_regression_2026-08-21.md)
- [DD/DG Id–Vg 残差分解数据](transportmodels_idvg_residual_decomposition_2026-08-21.json)
- [阶段 2–7 汇总](transportmodels_dg_phase2_to_phase7_summary_2026-08-21.md)
- [Frozen-Q oracle](transportmodels_dg_frozen_q_oracle_2026-08-20.md)
- [表面迁移率审计](transportmodels_dg_surface_mobility_audit_2026-08-21.md)
- [原始 Sentaurus Id–Vg deck](../../build-release/reference_tcad/transportmodels_sentaurus2022/run02/full_raw/IdVgs_des.cmd)
- Sentaurus Device User Guide T-2022.03：本机 `D:\工作\学习资料\TCAD软件手册\Sentaurus帮助文档\data\sdevice_ug.pdf`，重点页 362–369、403–407、426–451。
