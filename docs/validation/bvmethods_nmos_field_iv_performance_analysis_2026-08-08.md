# BVmethods NMOS 场量、I-V 与性能对比（2026-08-08）

## 结论

在漏极总电流 `1e-4 A/um` 的相同边界上，Sentaurus 新运行得到
`6.384111661907364 V`，Vela 得到 `6.395904200106065 V`，相对误差为
`0.184717%`。两者的电势拓扑、高场热点位置和 I-V 击穿拐点一致。

当前 Vela 的主要性能瓶颈是高场非线性收敛和外层边界算法重复调用完整
漂移扩散求解；单次 Newton 更新本身也明显慢于本次 Sentaurus 运行的有效均值。

## 场量与 I-V

- 电势逐节点绝对误差中位数：`0.000763928 V`。
- 电势逐节点绝对误差 P95：`0.015810852 V`。
- 节点电场图中，Sentaurus 节点向量幅值峰值为 `2.842592059e8 V/m`；
  Vela 三角形电势梯度面积加权节点重构峰值为 `2.769127564e8 V/m`。
- 上述节点电场峰值与此前边投影电场峰值使用不同离散定义，不可直接混合比较。
- 在 `1e-4 A/um` 判据电流处，两套结果电压差为约 `11.79 mV`。

图与对齐数据由 `scripts/plot_bvmethods_nmos_boundary_comparison.py` 生成，位于：

`build-release/reference_tcad/bvmethods_sentaurus2018/run01/sentaurus_boundary_state_20260808/report_20260808`

## 性能观测

| 运行 | 墙钟时间 | Newton 更新 | 每更新有效墙钟时间 |
|---|---:|---:|---:|
| Sentaurus 完整预偏置和电流边界 | 103.24 s | 370 | 0.2790 s |
| Vela 电压转电流边界恢复段 | 1032.37 s | 279 | 3.7003 s |
| Vela 外接电阻恢复段 | 3885.10 s | 1356 | 2.8651 s |

原始墙钟时间不能视为严格同口径基准：Sentaurus 包含从 0 V 开始的完整预偏置，
而两项 Vela 时间只覆盖已有检查点之后的高场段；两者运行环境也不同。即使在这个
对 Vela 更有利的范围下，Vela 两项恢复段仍分别为 Sentaurus 完整运行的约 `10.0x`
和 `37.6x`，每次 Newton 更新为 Sentaurus 有效均值的约 `13.3x` 和 `10.3x`。

Sentaurus 日志可拆分的 `77.75 s` 核心步骤中：

- RHS、残差与装配：`45.13 s`，`58.0%`；
- 线性求解：`25.67 s`，`33.0%`；
- Jacobian：`5.58 s`，`7.2%`；
- 其余：约 `1.37 s`，`1.8%`。

Vela 电流切换最终需要 2 次新的完整 DD 边界评估、共 279 次 Newton 更新；
外接电阻最终括区需要 7 次完整 DD 评估、共 1356 次更新。检查点恢复完成后，
约 99.9% 的记录时间位于这些新高场 Newton 求解中，故检查点 I/O 不是瓶颈。

## 代码路径判断

日志直接证明的瓶颈优先级为：

1. 高场 Newton 更新数过多；
2. 外接电阻的嵌套标量根导致多次完整 DD 求解；
3. 单次 Newton 更新成本偏高。

Vela 每次 Newton 会执行连续性闭合/行缩放、Jacobian 装配、稀疏数值分解和
线性求解，并在线搜索候选点上重新计算残差。`LinearSolver` 只缓存符号分析，仍在
每次调用执行 `SparseLU::factorize()`。不过 Vela 尚未记录这些阶段的独立耗时，
因此“雪崩源/Jacobian 装配”和“SparseLU 数值分解”的内部排序仍是待插桩验证的推断。
最终接受状态中，绝大多数 Newton 更新只需一次线搜索尝试，线搜索不是当前第一嫌疑。

## 优化顺序

1. 给残差、雪崩源/Jacobian 装配、连续性缩放、SparseLU 分解/求解、线搜索、
   诊断和检查点 I/O 加低开销分阶段计时。
2. 用括区恢复和割线/切线预测减少外层 DD 评估，并将外接电阻/电流接触改为
   分块或 Schur 联立约束，目标是去除嵌套标量根。
3. 在不调整迁移率或雪崩参数的前提下，改进高场初值预测、连续性缩放和
   continuation，先降低 Newton 更新数。
4. 再依据插桩结果缓存不变装配项、减少稀疏矩阵复制，并评估合适的稀疏求解器。
5. 最后在同机、同线程、同预偏置范围下重新测量端到端速度。

## 可复现产物

- `scripts/compare_bvmethods_nmos_boundary_fields.py`
- `scripts/plot_bvmethods_nmos_boundary_comparison.py`
- `scripts/build_bvmethods_nmos_performance_report.py`
- `report_20260808/field_node_comparison.csv`
- `report_20260808/iv_curve_comparison.csv`
- `report_20260808/runtime_summary.csv`
- `report_20260808/performance_summary.json`
- `report_20260808/artifact.json`
- `report_20260808/report.html`
