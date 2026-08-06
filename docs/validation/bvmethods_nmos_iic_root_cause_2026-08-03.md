# BVmethods NMOS IIC 差异根因调查（2026-08-03）

## 结论

当前差异不是由电场或 van Overstraeten 雪崩系数的统一比例误差主导，而是由两条相互独立的载流子路径问题叠加造成：

1. 0--6.40 V 后处理分支接受了局部载流子方程尚未闭合的状态。高偏压下空穴准费米势和空穴密度尾部已经严重偏离 Sentaurus，电子密度尾部也明显偏低；因此 SG 电流、端电流和雪崩源所需的电流支撑均不可信。
2. `DCSweep.cpp` 的后处理 `effectiveNi` 构造仍采用“第一个相邻单元材料”，与核心 `CoupledDDAssembler` 已实现的“优先输运材料”语义不一致。Si/SiO2 共享节点若先遇到氧化层单元，后处理会得到 `ni_eff=0`，使 Sentaurus 雪崩峰所在的界面 SG 边电流和 IIC 源被直接清零。

此外，原分析脚本把 Vela 原生雪崩源积分换算成 SI 时多乘了 `10^6`，造成了“Vela 雪崩源高 28 倍”的错误判断。修正后，6.38 V 的 Vela 生成率峰值实际只约为 Sentaurus 的 `2.83e-5`，方向与载流子/电流支撑不足一致。

## 证据链

### 1. 场与 alpha 不是第一处分叉

6.38 V 同几何边比较中：

- 最大电场：Sentaurus `2.2791e8 V/m`，Vela `2.4460e8 V/m`，比值 `1.073`；
- 电子 alpha 峰值：Sentaurus `4.4008e7 1/m`，Vela `4.6000e7 1/m`，比值 `1.045`；
- 修正单位后的雪崩生成率峰值：Sentaurus `5.5783e36 1/m3/s`，Vela `1.5767e32 1/m3/s`，比值 `2.8265e-5`。

因此源项缺失发生在 `alpha * |J| / q` 的电流支撑侧，而不是电场幅值侧。

### 2. 载流子状态从 1 V 起已经分叉

1909 个完全同坐标半导体节点的比较结果：

| Vd | electron QF p95 | hole QF p95 | electron density p95 | hole density p95 |
|---:|---:|---:|---:|---:|
| 1 V | 0.0945 V | 0.8517 V | 1.596 dex | 14.969 dex |
| 2 V | 0.2035 V | 1.8416 V | 3.474 dex | 31.764 dex |
| 4 V | 0.3866 V | 3.9518 V | 6.806 dex | 67.662 dex |
| 6 V | 0.5352 V | 6.0210 V | 9.201 dex | 102.286 dex |
| 6.38 V | 0.5552 V | 6.4537 V | 9.542 dex | 109.139 dex |

6.38 V 的典型异常：

- 节点 325：Sentaurus `n=2.202e22 m-3`，Vela `n=6.769e8 m-3`，低 `13.51 dex`；
- 节点 428：Sentaurus `p=3.141e19 m-3`，Vela `p=9.147e-135 m-3`，低 `153.54 dex`；
- 空穴准费米势最大绝对误差达到 `12.589 V`。

这解释了 Vela 空穴端电流几乎为零，并使电子 SG 漂移/扩散大项依赖极端抵消得到很小的净电流。

### 3. 原 6.38 V 状态并非严格收敛态

原探针配置使用 `abstol=1e-5`、载流子逐行收敛检查关闭、QF 越界仅告警。6.38 V 状态以 `initial_abstol`、0 次 Newton 迭代被接受，同时仍有 180 个 QF 越界节点。

使用相同初态将 `abstol` 收紧到 `1e-12` 后：

- 初始残差约 `8.86e-6`，不再被接受；
- 第 1 次 Newton 后组合残差约 `4.67e-6`；
- 分块残差由空穴方程主导：`phip=4.665e-6`，而 `phin=1.792e-7`、`psi=3.793e-7`；
- 第 2 次迭代以 `line_search_non_decrease` 失败，13 次回溯均未接受，Newton 步长为 `4.821`。

保持原 `abstol`、只开启载流子逐行强制检查，也会得到 658 个不合格载流子行、最大比值 `1.819`，随后同样在线搜索处失败。说明原分支是由全局绝对阈值放行的非闭合状态，而不是可作为 IIC 基准的物理解。

### 4. 后处理 `ni_eff` 与核心装配器语义不一致

核心 `CoupledDDAssembler` 通过 `buildValidatedEffectiveNodeNi()` 调用 `buildNodeNi()`；后者会在共享节点上优先选择 `ni>0` 或迁移率非零的输运材料。

但 `DCSweep.cpp::buildEffectiveIntrinsicDensityVector()` 只在 `!seen[nodeId]` 时赋值，之后不允许 Si 覆盖先遇到的 SiO2。该向量被传给：

- `sgEdgeCurrentAvalancheSourceRecords()`；
- `sg_avalanche_edges.csv`；
- release BV audit、terminal-current method compare；
- postprocess-only IIC 的雪崩源积分审计。

网格邻接已证实：

- 全网格恰有 51 个同时邻接 Si 与非 Si 材料的节点；按三角形顺序，它们
  的第一个材料全部不是 Si，因此旧实现会系统性地把这 51 个输运界面节点
  赋成零 `ni_eff`；
- 边 208（节点 1381--267）位于 `y=0`，两端节点同时属于先出现的 `R.Gateox/SiO2` 和后出现的 `R.Substrate/Si`；后处理记录两端 `ni=0`；
- 边 616（节点 325--1381）同样位于 Si/SiO2 界面，后处理 SG 电流被清零；
- 这些边处在 Sentaurus 的高生成率表面区域；
- 体内边 2436（节点 428--426）只邻接 Si，Vela 的场和 alpha 与 Sentaurus 接近且能产生非零雪崩源。

因此这不是共享节点物理模型本身缺失，而是后处理存在一份未与核心装配器统一的旧节点材料归属实现。

### 5. 接触峰值边是另一类伪峰

Vela 最大电子电流密度出现在源接触边 6361（节点 2556--2592）。该边电场只有约 `1.25e5 V/m`、alpha 近零，因此不会主导雪崩积分；但共享接触节点 2556 在后处理 `effectiveNi` 中为零，使 SG 分解出现异常大电流峰。它会污染“最大电流密度”统计，应与真正的高场雪崩边分开处理。

## 单位修正

Vela 原生诊断量采用：alpha 为 `cm-1`、粒子通量为 `cm-2 s-1`、二维支撑面积为 `um2`。因此：

- 原生边源积分到线源：乘 `1e-6`，单位为 `m-1 s-1`；
- 每 `um` 深度的雪崩电流：`Iava = q * raw_source_integral * 1e-12`；
- 原生粒子通量到 `m-2 s-1`：乘 `1e4`。

修正后的 6.38 V Vela 数据为：

- `Id = 5.349742e-9 A/um`；
- `Iava = 2.918131e-9 A/um`；
- `Iava/Id = 0.54547`。

原先的 `Iava=2.918e-3 A/um` 和“源项高 28 倍”均来自固定 `10^6` 的换算错误，不应继续作为物理调参依据。

同一量纲假设也存在于 `DCSweep.cpp` 的 release BV audit 和
`avalanche_internal_source_current_audit`：它们将原生
`record.edgeSourceIntegral` 直接按 SI 线源乘 `q*1e-6`。在 unit-scaling
模式下还必须先乘原生源到 SI 线源的 `1e-6`，因此这些 C++ audit 当前
报告的 `qG_A_per_um` 同样高 `10^6`。本轮只修正了离线分析脚本并增加单位
回归测试；在 C++ audit 修复前，不能用其绝对 `qG` 数值作验收。

## 后续修复优先级

1. 统一 `DCSweep` 与 `CoupledDDAssembler` 的 `effectiveNi` 构造，仅保留一个经过共享节点回归测试的实现；先重跑 postprocess-only IIC，确认边 208/616 不再为零。
2. 将载流子逐行收敛和 QF 物理界限纳入分支接受条件，不能再由 `initial_abstol` 单独放行；同时让线搜索失败后的载流子恢复覆盖当前 `4.7e-6` 残差区间。
3. 从 1 V 起重建连续偏压分支，首先闭合空穴方程，再比较电子/空穴 QF、密度和 SG 边通量；不要用 6.38 V 的失真状态直接热启动。
4. 上述两项完成后重新计算 IIC 交点，再进入外接电阻和 voltage-to-current 方法闭合；在此前不调整迁移率或雪崩系数拟合端电流。

## 可复现产物

- 同节点载流子汇总：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_postprocess_20260803/analysis/multibias_sentaurus/same_node_carrier_state_summary.csv`
- 同节点载流子明细：同目录 `same_node_carrier_state.csv`
- 同边场/电流/alpha 明细：同目录 `matched_edge_fields.csv`
- 6.38 V 严格重放：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rootcause_20260803/strict_replay/postprocess_only`
- 载流子逐行强制重放：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rootcause_20260803/carrier_row_enforce/postprocess_only`
- 单位回归测试：`tests/regression/test_bvmethods_nmos_iic_units.py`
