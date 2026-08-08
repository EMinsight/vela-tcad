# BVmethods Sentaurus 载流子电流方向 IIC 候选路径恢复

## 结论

已按 Sentaurus `WriteAll` 输出的三档 `Maximum Field` 约束重新搜索局部场峰种子，并分别沿
`eCurrentDensity-V` 与 `hCurrentDensity-V` 双向追踪候选流线。候选只有在同时恢复对应
`eIonIntegral`、`hIonIntegral` 和 TDR `MeanIonIntegral` 平台后才保留。

本轮保留的参考候选为：rank-1 使用空穴电流方向，rank-2 与 rank-3 使用电子电流方向。
此前仅用 `ElectricField-V` 方向得到的三条代理线不再作为 IIC 几何参考。

## 搜索范围与验收条件

- 数据源：`iic_v10p448267_0000_des.tdr`，终态漏极电压 10.4482667308 V。
- 从 188 个节点局部电场极大值中，按三档 WriteAll 最大场各取 8 个最近局部峰。
- 对局部峰节点以及其相邻单元内的小偏移点分别生成电子、空穴电流方向流线。
- 实际追踪 164 条流线：48 条顶点种子、90 条 rank-3 单元内偏移种子、26 条 rank-1/rank-2 定向补扫种子。
- 验收阈值：种子节点场相对误差不超过 5%；三个积分平台中任一相对误差不超过 1%。
- `MeanIonIntegral` 使用 TDR 中独立存储的平台值；它不是 WriteAll 中 e/h 两个积分的简单算术平均。

## 最终保留路径

| Rank | 候选 | 种子节点/单元 | 种子坐标 (um) | 方向场 | 种子场误差 | eIon 目标/恢复 | hIon 目标/恢复 | Mean 目标/恢复 | 最大积分误差 |
|---:|---|---|---|---|---:|---:|---:|---:|---:|
| 1 | C0001 | node 924 | (0.17656250, 0.00224609) | hCurrentDensity-V | 0.001% | 1.71545 / 1.714549 | 1.90119 / 1.900193 | 1.794564 / 1.793625 | 0.053% |
| 2 | C0058 | node 98 | (0.03531250, 0.00077881) | eCurrentDensity-V | 0.779% | 1.39435 / 1.394349 | 1.76272 / 1.762718 | 1.544187 / 1.544187 | 0.0001% |
| 3 | C0125 | node 1375, cell 2917 nudge | (0.06500081, 0.00048545) | eCurrentDensity-V | 1.325% | 1.28569 / 1.283154 | 1.65037 / 1.647116 | 1.428772 / 1.425956 | 0.197% |

## 排除的歧义候选

- Rank-2：最大场最接近的局部峰 node 2184 沿电子方向只得到
  `0.90917 / 1.14936 / 1.00687`，不能恢复 rank-2 平台，因此排除；node 98 虽然是该档
  第 2 近的局部峰，但三项积分均精确进入目标平台，因此保留。
- Rank-3：node 285 的种子场仅高于目标 1.029%，但电子/空穴方向都进入 rank-2 平台
  (`1.39435 / 1.76272 / 1.54419`)，因此排除；node 1375 的单元内偏移种子恢复 rank-3
  平台，因而保留。
- Rank-1：node 924 的空穴方向同时满足近乎精确的最大场约束和三项平台约束；同一点的
  电子方向平台不足，因此只保留空穴方向。

## 产物

- 全部种子：`daily_report_2026-08-07/sentaurus_iic_candidate_seeds_20260807.csv`
- 全部已评分扫描：`daily_report_2026-08-07/sentaurus_iic_candidate_scan_scored_20260807.csv`
- 三条最佳候选：`daily_report_2026-08-07/sentaurus_iic_candidate_best_selected_20260807.csv`
- 三条有序路径：`daily_report_2026-08-07/sentaurus_iic_selected_path_rank{1,2,3}_20260807.csv`
- 网格与积分剖面图：`daily_report_2026-08-07/bvmethods_sentaurus_current_directed_iic_paths_20260807.png`

## 限制与后续使用方式

WriteAll 最大场约束用于筛选种子节点处的 TDR 电场，不应再用流线所穿过的全域最大总电场
替代。Sentaurus Visual 电流密度流线仍不是 sdevice 内部 `ComputeIonizationIntegrals` 有序折线的
直接导出接口；但本轮候选同时满足载流子方向、种子场和三项积分平台，已经可以作为 Vela
路径转向点、空间支撑和路径排序的几何对比参考。后续不应把未恢复目标平台的流线用于评价
Vela 路径误差。

## 后续几何对比补充

将 Visual 流线步数上限从 5000 提高到 50000 后，rank-2 和 rank-3 的完整电流流线会继续进入
更高的离化积分平台，证明无限延伸的 Visual 流线不能直接代表独立 IIC rank。几何对比现采用
围绕种子、且 e/h/Mean 均不超过目标平台 0.5% 的连通段；完整扩展流线保留为诊断数据。
详见 `bvmethods_sentaurus_vela_iic_path_comparison_2026-08-07.md`。
