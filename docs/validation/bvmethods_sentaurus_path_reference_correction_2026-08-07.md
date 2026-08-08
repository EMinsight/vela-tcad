# BVmethods Sentaurus 路径参考修正

## 技术结论

当前图中 Vela 连续路径与 Sentaurus Visual 流线的巨大几何差异，主要不能归因于
Vela 路径算法。此前构造的 Sentaurus 参考线不是 `ComputeIonizationIntegrals`
内部路径：官方 `pp4_des.cmd` 使用 `Avalanche(Eparallel)`，而参考线使用
Sentaurus Visual 的 `ElectricField-V`；同时，TDR 积分平台上的最大电场节点并不
等于对应 WriteAll 路径的内部种子。

## 数值证据仍支持前三条 IIC 路径接近

比较终态 10.4482667308 V 下按积分值排序的前三条路径：

| Rank | Sentaurus Imean | Vela Imean | 相对误差 | Sentaurus log 最大场 (V/cm) | Vela 分组种子场 (V/cm) | 相对误差 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.808320 | 1.788127 | -1.12% | 1.8871e6 | 1.881244e6 | -0.31% |
| 2 | 1.578535 | 1.542601 | -2.28% | 3.4187e5 | 3.385238e5 | -0.98% |
| 3 | 1.468030 | 1.477063 | +0.62% | 3.6249e6 | 3.568072e6 | -1.57% |

因此，积分值和峰值场量级已经接近；尚未闭合的是内部路径坐标、载流子专属方向和
终止语义。

## ElectricField-V 代理为何失效

从 Sentaurus Visual 导出的三条代理流线重新采样原始积分场后得到：

| 代理线 | 目标 Imean | 流线上实际最大 MeanIonIntegral | 结论 |
|---:|---:|---:|---|
| R1 | 1.808 | 0.580 | 未穿过目标积分平台 |
| R2 | 1.579 | 1.544 | 落到 rank-2 平台 |
| R3 | 1.468 | 1.544 | 与 R2 落到同一平台 |

另外，Sentaurus WriteAll 路径 3 的最大场为 `3.4187e5 V/cm`，但此前从 rank-2
积分平台挑出的“最大场节点”为 `2.1677e6 V/cm`，相差 6.34 倍。这直接否定了
“平台内最大场节点就是路径种子”的假设。

## 当前路径定义

- Sentaurus 控制算例：`Avalanche(Eparallel)`、`ComputeIonizationIntegrals`、
  `BreakAtIonIntegral(3 1.)`。
- Vela 当前终态：`continuous_cell`、`sentaurus_eparallel_adaptive`、
  `numbered_peak_groups`，路径停止场为 0。
- Sentaurus Visual 代理：TDR 插值后的 `ElectricField-V` 双向 RK4 流线，不能替代
  上述内部 Eparallel IIC 路径。

## 后续验证顺序

1. 以 WriteAll 每条路径的 `Maximum Field` 为约束，在 TDR 的局部峰候选中重新识别
   路径种子，不再从积分平台直接取最大场节点。
2. 分别用 `eCurrentDensity-V`、`hCurrentDensity-V` 及载流子专属 Eparallel 方向生成
   候选流线，并沿线回采 `eIonIntegral/hIonIntegral/MeanIonIntegral`；只有恢复目标平台
   和日志积分的候选才保留。
3. 对保留下来的路径比较种子坐标、转向点、长度、停止原因和逐段离化积分，再决定
   是否需要修改 Vela 的峰分组、方向插值或停止条件。
4. 在完成上述闭合前，不使用 ElectricField-V 代理图评价 Vela 几何误差。

## 限制

Sentaurus TDR 和 WriteAll 日志均不直接公开 `ComputeIonizationIntegrals` 的内部有序
折线。Sentaurus Visual `extract_streamlines` 能导出可视化流线坐标，但其输入向量场、
插值和积分器与 sdevice 内部 IIC 路径追踪并非同一接口，因此只能作为候选路径验证工具。
