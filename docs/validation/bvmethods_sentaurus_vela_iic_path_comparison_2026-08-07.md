# BVmethods Sentaurus 与 Vela IIC 路径对比

## 结论

三条路径的积分终值接近，并不代表几何路径已经闭合：

- Rank-2 是当前几何最接近的一条，Vela 路径的 72.4% 位于 Sentaurus 候选路径 10 nm 内，
  对称平均距离为 21.7 nm。
- Rank-3 的种子位置和局部高场支撑已经基本闭合：种子相差 1.25 nm，Sentaurus 候选路径
  全部位于 Vela 路径 10 nm 内；但 Vela 路径继续向衬底延伸，只有 3.1% 的 Vela 路径位于
  Sentaurus 候选 10 nm 内。这是明显的路径停止/物理通道保留差异。
- Rank-1 不只是低场尾部不同：两个种子相差 87.3 nm，10 nm 双向覆盖率均低于 10%，说明
  Vela 当前选中了场强相近、但空间位置不同的局部峰或物理通道。

## 数据和比较口径

- 偏压：10.4482667308 V。
- Vela：`adaptive_minority_qf_branch_20260806/bias_10p448267/postprocess_only`，前三个数值 rank。
- Sentaurus：按 WriteAll 最大场与 e/h/Mean 平台筛选的 current-density Visual 候选。
- Rank 映射：按终态离化积分平台的数值次序匹配，而不是按路径坐标人工配对。
- Sentaurus Visual 流线步数上限从 5000 提高到 50000；三条流线分别在积分器索引
  8158、46441、14866 处自然终止，不再受原始 5000 步上限截断。
- 扩展后的 rank-2/rank-3 电流流线会进入更高积分平台，因此不能把整条流线当成独立 rank。
  本报告围绕种子保留 e/h/Mean 均不超过目标平台 0.5% 的连通段。完整扩展流线单独保留，
  用于证明 Sentaurus 内部 IIC 具有物理通道分离、停止或分叉语义。
- 几何距离使用两条折线各 1001 个等弧长采样点。`S near V` 和 `V near S` 分别表示
  Sentaurus/Vela 采样点到另一条路径距离不超过 10 nm 的比例。

## 定量结果

| Rank | 种子距离 | Sentaurus/Vela 长度 | Vela/Sentaurus 长度比 | 对称平均距离 | Hausdorff 距离 | S near V (10 nm) | V near S (10 nm) | Mean 积分差异 |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 87.3 nm | 0.075 / 1.117 um | 14.86x | 228.1 nm | 951.5 nm | 9.6% | 1.7% | -0.31% |
| 2 | 53.0 nm | 0.345 / 0.215 um | 0.62x | 21.7 nm | 145.5 nm | 48.6% | 72.4% | -0.10% |
| 3 | 1.25 nm | 0.029 / 1.348 um | 46.32x | 231.5 nm | 1055.8 nm | 100.0% | 3.1% | +3.58% |

路径长度和 Hausdorff 距离会放大低场尾部差异；种子距离和双向 10 nm 覆盖率更适合判断
高场通道是否相同。因此 rank-3 不能简单归类为“整条路径完全错误”：其高场起始段正确，主要
错误发生在离开热点后的延续与终止。Rank-1 则不同，起始通道本身尚未闭合。

## 分 rank 诊断

### Rank-1：峰分组/路径排序优先于停止条件

- Sentaurus 候选种子：`(0.1765625, 0.0022461) um`，空穴电流方向。
- Vela 种子 node 990：`(0.1677344, 0.0891094) um`。
- 两者峰值场量级接近，但空间位置相差 87.3 nm；说明仅按峰值大小排序不足以恢复同一通道。
- 后续应先修正局部峰显著度、鞍点连通和物理通道排序，再讨论低场停止长度。

### Rank-2：主体通道最接近，但种子与端点仍不同

- 两条路径在表面横向通道上有明显重合，Vela 路径 72.4% 位于 Sentaurus 候选 10 nm 内。
- 种子沿同一表面通道相差约 53 nm，说明 rank 标识所对应的“代表峰”仍不一致。
- Vela 路径长度约为平台裁剪后 Sentaurus 候选的 62%，应核对电子方向的转向点、分支选择和
  路径终止位置。

### Rank-3：种子与局部走向闭合，低场长尾过长

- 种子仅相差 1.25 nm，Sentaurus 候选全部位于 Vela 路径 10 nm 内。
- Vela 继续沿衬底方向延伸至约 1.35 um，而目标平台连通段约 0.029 um。
- 这表明峰搜索和局部电子方向插值已基本正确，剩余问题主要是路径离开强离化区后的停止场、
  通道合并或 rank 独立保留语义。

## 产物

- 全路径叠加图：`daily_report_2026-08-07/bvmethods_sentaurus_vela_iic_path_compare_20260807.png`
- 高场区局部图：`daily_report_2026-08-07/bvmethods_sentaurus_vela_iic_path_hotspot_compare_20260807.png`
- 几何指标：`daily_report_2026-08-07/sentaurus_vela_iic_path_geometry_compare_20260807.csv`
- 平台裁剪记录：`daily_report_2026-08-07/sentaurus_iic_platform_clip_summary_20260807.csv`
- 平台连通段：`daily_report_2026-08-07/sentaurus_iic_selected_path_rank{1,2,3}_20260807.csv`
- 未裁剪扩展流线：`daily_report_2026-08-07/sentaurus_iic_selected_path_rank{1,2,3}_extended_20260807.csv`

## 限制

Sentaurus 路径仍是 Visual 中由 `eCurrentDensity-V`/`hCurrentDensity-V` 生成、并用目标积分平台
约束的候选，而不是 sdevice 内部 `ComputeIonizationIntegrals` 折线的直接导出。因此上述几何
误差适合定位 Vela 的峰选择、方向延续和停止语义，不应表述为对 sdevice 内部坐标的精确误差。
