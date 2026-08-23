# TransportModels Id-Vg 深关断数值审计（2026-08-22）

## 结论

既定的五项任务均已实施。reference-aware 准费米坐标错误已经修复；深关断载流子行和全局连续性已经成为硬门槛；DG 热点的解析 Jacobian 与中心有限差分相符；增量准费米表示足以保存低于绝对电势 ULP 的电势差；Newton 参数完成受控标定。

硬门槛仍会拒绝 DD/DG 的若干深关断点。这是预期的保护行为，不能通过提高 stall floor 或减小步长伪装为收敛。剩余问题已经定位为局部连续性闭合与全局接触通量—体源积分闭合之间的不一致，而不是 DG 解析 Jacobian 缺项或端口漂移—扩散相消。

## 1. Reference-aware gauge

- 混合材料非输运节点的电子/空穴 gauge 残差改为物理绝对准费米势，即内部增量加对应 reference。
- 节点源 Jacobian 的载流子密度导数改为使用 reference-relative 准费米变量。
- 新增真实 NMOS/PMOS 混合 Si/SiO2 网格不变性测试，比较残差、完整 Jacobian、连续性分解和氧化层 gauge。
- 完整 `test_mos_mixed_material`：9 个测试、1437 个断言全部通过。

## 2. 深关断硬门槛

严格工作流同时强制：

- `carrier_row_convergence.mode = enforce`
- `eps_row = 1e-3`
- `min_source_scale = 1e-18`
- `global_continuity_closure.mode = enforce`
- `global_continuity_closure.tolerance = 0.1`
- `carrier_row_qualified_stall_acceptance = true`

六个单点中，DD -1 V 通过；DD -0.68/-0.52 V 与 DG -1/-0.68/-0.52 V 被硬门槛拒绝。失败点不再写入可用于 Sentaurus 对比的有效曲线。

## 3. DG 热点 Jacobian 审计

审计状态为 DG -0.4 V 最后收敛状态，热点节点为 706、705、794、342、743、712、714、711、793、716。探针覆盖 `psi`、`phin` 和 `psi_minus_phin` 热点簇方向，并对前四个节点单独审计电子准费米方向。

首次审计发现 JVP 工具错误地把方向向量当成状态打包，错误减去了准费米 reference。修复后扰动范数随步长正确缩放，并加入 reference 坐标不变性回归。

| 扰动 | 热点簇电子方向尺度化误差 | 判断 |
|---:|---:|---|
| 1e-4 V | 2.601e-2 | 有限差分非线性截断误差明显 |
| 1e-6 V | 1.624e-4 | 进入线性区 |
| 1e-8 V | 9.830e-7 | 与解析 Jacobian 一致 |

四个单节点在 1e-8 V 下的误差为 7.62e-7 至 5.28e-6。误差随步长下降，未发现 DG 热点载流子 Jacobian 缺项。

## 4. 增量准费米与局部精度评估

接触参考坐标复算后，DD/DG -1 V 的电子电流边均不再归零，四端 KCL 绝对误差分别改善到 5.34e-19 A/um 和 4.35e-19 A/um。求解状态、连续性通量和端口电流均使用准费米增量；SG 通量使用平衡的 `expm1` 形式，已有 long-double 参考诊断。

新增回归在 1.1 V 绝对准费米电势上施加 5e-18 V 增量。物理输出中的两个电势相等，但增量状态和端口电流仍保持非零，残差法与边通量法一致。

因此当前不引入新的局部多精度未知量。DD -0.68 V 等失败点仍存在数百个局部连续性违规，KCL 误差不是端口电流累加舍入造成的。

## 5. Newton floor 与 line search 标定

在相同硬门槛下运行 18 个单点：DD/DG、-1/-0.68/-0.52 V，以及三组参数组合。

| 组合 | 结果 |
|---|---|
| qf cap 0.01 V，全步 | 仅 DD -1 V 通过；其余为全局连续性失败 |
| qf cap 0.005 V，全步 | 未提高通过率，多个点局部违规更多 |
| qf cap 0.01 V，初始阻尼 0.5 | DD -0.68 V 局部违规降到 0，但 74 次迭代后仍全局连续性失败；DG -0.52 V 从 17 次恶化到 91 次 |

严格工作流最终采用：

- `stall_residual_floor = 2e-11`
- `line_search = true`
- `damping_factor = 1.0`
- `quasi_fermi_update_limit_V = 0.025`

floor 位于已通过点观测到的 1.54e-11 至 1.89e-11 数值平台之上，并且只有局部载流子行与全局连续性同时通过时才能接受 stall。

## 验证与数据

- `test_newton_solver`：88 个测试、1196 个断言全部通过。
- `test_mos_mixed_material`：9 个测试、1437 个断言全部通过。
- 严格配置 Python 回归：2 个测试全部通过。
- JVP 数据：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/dg_hotspot_jvp_audit/`
- Newton 标定：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/newton_calibration_summary.csv`
- reference A/B：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/qf_reference_ab_summary.csv`

## 后续问题

下一步不应继续放宽 Newton 参数。优先审计“局部载流子行已满足但全局连续性仍失败”的 DD -0.68 V 状态，逐项核对接触节点连续性残差、接触边 SG 通量、接触边界行替换以及 SRH 体积分的符号和控制体权重。
