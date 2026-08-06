# BVmethods NMOS postprocess-only IIC 闭合报告（2026-08-03）

> **更正（2026-08-03）：** 本文件的早期分析存在雪崩源积分单位换算错误，
> 且未识别 `DCSweep` 后处理 `effectiveNi` 的共享节点材料归属问题。有关
> “Vela 雪崩源高 28 倍”的结论已经作废。请以
> `docs/validation/bvmethods_nmos_iic_root_cause_2026-08-03.md` 为准；修正后
> 6.38 V 的 Vela 生成率峰值约为 Sentaurus 的 `2.8265e-5`，主要差异方向是
> 载流子/SG 电流支撑不足。

## 结论

Sentaurus 虚拟机已恢复连接。基于官方 `BVmethods/pp4_des.cmd` 的 IIC 物理设置，已经重新生成、下载并导入精确多偏压 TDR，覆盖：

- `1/2/4/5/6 V`；
- `6.32/6.34/6.36/6.37/6.38/6.39/6.40 V`；
- 用于包围精细电流交点的 `6.42/6.45/6.50/6.60/6.70/6.80/6.90/7.00/7.10 V`。

因此，原先“虚拟机无响应、缺少多状态 TDR”的限制已经解除。21 个状态均包含电势、准费米势、载流子、迁移率、电场、电流密度、电子/空穴雪崩系数、离化积分和 ImpactIonization。

当前 Vela 结果仍未闭合 IIC 基准。Sentaurus 精确状态表明电场峰值已经比较接近，但 Vela 的电子/空穴电流分支和雪崩源空间积分仍存在决定性误差，不能通过统一缩放雪崩系数解决。

## Sentaurus 精确检查点

| Vd (V) | Id (A/um) | Iava (A/um) | Iava/Id | Phi electron | Phi hole |
|---:|---:|---:|---:|---:|---:|
| 1.00 | 2.993217e-9 | 1.864703e-10 | 0.062298 | 0.212956 | 0.093862 |
| 2.00 | 3.735698e-9 | 9.707659e-10 | 0.259862 | 0.350896 | 0.203823 |
| 4.00 | 2.673564e-7 | 1.523612e-7 | 0.569881 | 0.675587 | 0.510957 |
| 5.00 | 1.934757e-6 | 1.458475e-6 | 0.753829 | 0.813925 | 0.718452 |
| 6.00 | 6.326760e-6 | 5.530908e-6 | 0.874209 | 0.991509 | 0.940903 |
| 6.38 | 9.301763e-6 | 8.691665e-6 | 0.934411 | 1.051093 | 1.019594 |
| 6.70 | 1.262131e-5 | 1.253246e-5 | 0.992960 | 1.105374 | 1.086027 |
| 6.80 | 1.383571e-5 | 1.400495e-5 | 1.012233 | 1.122999 | 1.107055 |

在 `6.70–6.80 V` 两个精确检查点之间对 `Iava-Id` 线性插值，电流交点为 `6.734425890 V`。

这与官方 Workbench 稀疏曲线提取的 `6.377494278 V` 不是同一个采样结果：官方值必须继续保留为“官方稀疏提取基准”，`6.734425890 V` 则记录为“当前精确检查点电流交点”。在完全复现官方 CurrentPlot 采样点、BreakAtIonIntegral 停止条件和提取表达式以前，不应用其中一个覆盖另一个。这个差异已经从求解器问题收敛为可复现的采样/提取约定问题。

## 同网格空间比较

Sentaurus TDR 与 Vela 网格的 1909 个半导体节点坐标完全一致，最大坐标误差为 `0 um`。电势可直接做同节点比较；电场、电流密度和 alpha 使用同一几何边比较：Sentaurus 节点向量投影到边方向，标量取边两端平均，Vela 使用 `sg_avalanche_edges.csv` 中已换算为 SI 的量。

`6.38 V` 的同边结果如下：

| 物理量 | Sentaurus 峰值 | Vela 峰值 | Vela/Sentaurus | 中位数绝对 log10 比误差 | 相关系数 |
|---|---:|---:|---:|---:|---:|
| 电场 (V/m) | 2.279055e8 | 2.445959e8 | 1.07323 | 0.0463 dex | 0.9970 |
| electron alpha (1/m) | 4.400812e7 | 4.599973e7 | 1.04526 | 1.3988 dex | 0.5855 |
| hole alpha (1/m) | 3.498695e7 | 6.663211e7 | 1.90448 | 1.9442 dex | 0.6614 |
| electron current density (A/m2) | 1.369125e9 | 4.498811e14 | 3.28590e5 | 1.1797 dex | 0.0256 |
| hole current density (A/m2) | 4.972675e10 | 1.212508e4 | 2.43834e-7 | 13.5428 dex | -0.0038 |
| avalanche generation (1/m3/s) | 5.578268e36 | 1.576702e38 | 28.2651 | 5.3497 dex | -0.0372 |

这组结果说明：

1. 电场峰值和整体空间形状已经较好对齐；
2. electron alpha 的峰值接近，但空间分布仍不一致；
3. hole alpha 峰值约高 1.9 倍；
4. 最大误差来自电子/空穴 SG 电流支撑，继而放大雪崩源积分；
5. 下一阶段仍应优先修复异常接触边和电子/空穴输运分支，不应先拟合 van Overstraeten 系数。

同节点电势的中位绝对误差在全部偏压下约 `11.4–11.6 mV`，但 95 分位约为 `0.899–0.901 V`，最大值固定为 `0.901156 V`。这表明大部分节点已接近，同时仍有一组接触/界面节点存在系统性偏移，需要按区域和边界类型继续分解。

## Vela VTK 缩放限制

当前 `DCSweep.cpp` 在写多偏压 VTK 时没有把 `sweep.scaling` 传给 `writeDDSolutionVTK`，因此启用 `unit_scaling` 的运行会退回默认 legacy scaling。由此产生的 VTK 派生节点电场、alpha 和电流密度不能作为有效比较数据；电势和准费米势状态本身仍有效。

本报告已排除这些无效派生量，并删除旧的 `field_summary.csv`、`same_node_fields.csv` 和 `top_same_node_errors.csv`。修正该调用并重跑 Vela 探针后，才能补上严格“同节点电场/alpha”比较；在此之前，以同边 `sg_avalanche_edges.csv` 比较为可信依据。

## 后续闭合顺序

1. 将 `sweep.scaling` 传入多偏压 VTK 写出路径，增加 unit-scaling 回归测试，并重跑 1/2/4/5/6/6.32–6.40 V Vela 探针。
2. 对电势尾部 `0.899–0.901 V` 的节点按材料、接触、共享界面节点分类，确认是否来自 metal-gate/Barrier/flatband 或共享节点边界语义。
3. 固定异常接触边，逐边核对电子和空穴准费米势差、Bernoulli 参数、迁移率、广义 SG 通量及 conventional-current 符号。
4. 在端电流和同边电流密度闭合后再检查 alpha、ImpactIonization 与积分雪崩电流。
5. 最后分别复现官方稀疏 IIC 提取和密集检查点电流交点，明确 BV 验收采用的采样与插值规范。

## 可复现产物

- Sentaurus 精确多偏压执行脚本：`scripts/run_bvmethods_nmos_multibias_sentaurus_vm.py`
- 同节点/同边比较脚本：`scripts/compare_bvmethods_nmos_iic_multibias_fields.py`
- 21 个导入状态：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/sentaurus_iic_multibias_exact_extended_20260803/imported`
- 精确扩展曲线：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_postprocess_20260803/analysis/multibias_sentaurus/sentaurus_exact_extended_curve.csv`
- 同边汇总：同目录 `matched_edge_field_summary.csv`
- 同边逐条数据：同目录 `matched_edge_fields.csv`
- 同节点电势汇总：同目录 `same_node_potential_summary.csv`
- 机器可读执行摘要：同目录 `result.json`
