# BVmethods NMOS 节点星形 Eparallel 高压分支重算

## 执行范围

- 模式：`impact_ionization.eparallel_field_recovery=nodal_vertex_star`。
- 偏压：6.4、6.5、6.6、6.7、6.8、6.9、7.0、7.1 V。
- 固定项：SG 输运、Fermi-Dirac、OldSlotboom BGN、Masetti 高场迁移率、BTBT E2、van Overstraeten 参数和端电流提取方法。
- 雪崩耦合：`postprocess_only`，因此场恢复模式不反馈改变 DD 状态。

连续运行因工具 240 秒时限分为四段，但每段均从上一段接受态继续；八个目标偏压均已接受并生成状态、端电流和 SG 雪崩边记录。

## 电流源交点

| 口径 | 交点 / V |
|---|---:|
| Vela `nodal_vertex_star` | 7.053274758 |
| Vela 旧 `edge_adjacent_cells`（同一 DD 状态） | 6.920073790 |
| Sentaurus 密集同偏压曲线 | 6.734425890 |
| Sentaurus 官方稀疏 ABA 线性插值 | 6.377494278 |

- 新模式相对 Sentaurus 密集交点高 `+0.318848868 V`。
- 新模式相对旧 Vela 模式高 `+0.133200968 V`。
- 6.4 V 时新模式 `Iava/Id=0.89004`，7.0 V 时为 `0.99051`，7.1 V 时为 `1.00764`。
- Vela/Sentaurus 雪崩电流比从 6.4 V 的 `0.96449` 降至 7.0 V 的 `0.95069`。

官方 6.377494278 V 是稀疏 ABA-coupled 点的线性插值结果，不能与密集同偏压交点混用；本轮主要使用 6.734425890 V 作为同口径参考。

## 热点量闭合

在 6.4--7.0 V 的共同偏压点：

- 热点位置距离保持 `0.000147 um`；
- 峰值电子 `alpha` 比 Vela/Sentaurus 为 `0.98842--0.98888`；
- 峰值电子电流密度比为 `0.93574--0.93321`；
- 峰值电子产生率比为 `0.89856--0.89710`；
- 电子积分源比由 `0.94648` 降至 `0.93134`；
- 有效热点面积比由 `1.05333` 改善至 `1.03817`。

因此节点星形模式已经闭合局部 `Eparallel/alpha` 的空间语义，但电子电流支撑仍低约 6.5%--6.8%，积分电子源随偏压增长仍比 Sentaurus 慢。

## 产生率支撑

### P1 节点投影面积比（Vela/Sentaurus）

| 峰值阈值 | 6.4 V | 7.0 V |
|---:|---:|---:|
| 10% | 1.0833 | 1.0012 |
| 30% | 1.0825 | 1.1011 |
| 50% | 1.1361 | 1.0935 |
| 80% | 1.5155 | 1.1580 |

80% 阈值对离散投影和峰值归一化高度敏感；10%--50% 支撑更适合判断肩部趋势。

### 6.4 到 7.0 V 支撑面积增长

| 阈值 | Sentaurus | Vela P1 投影 | Vela 原生边支撑 |
|---:|---:|---:|---:|
| 10% | 1.1752 | 1.0861 | 1.1023 |
| 30% | 1.0882 | 1.1069 | 1.1247 |
| 50% | 1.0379 | 0.9989 | 1.1211 |
| 80% | 1.1507 | 0.8792 | 1.1132 |

径向累计源在所有共同偏压和检查半径上的最大差异约为 `2.50` 个百分点。整体空间分布已较接近，但 0.01--0.04 um 中尺度区域的 Vela 累计源略偏高，外侧肩部随偏压扩展仍不完全一致。

## 结论与下一步

1. 新节点星形恢复不应回退：它已复现 Sentaurus 节点场插值，并消除了原有 `Eparallel` 增长偏慢。
2. 旧模式交点较低，部分来自边两单元窄模板对肩部场和 `alpha` 的高估；因此旧交点更接近 Sentaurus 并不代表离散语义更正确。
3. 当前主要剩余项为：
   - 节点 `alpha(Eparallel)` 与边平均后求 `alpha` 的非线性交换误差；
   - 电子 SG 电流空间支撑仍低约 6.5%--6.8%；
   - 积分源在 0.01--0.04 um 区域的偏压增长和各向异性扩展。
4. 下一步应在节点星形模式上实现/验证“节点场 + 节点电流 + 节点 alpha + P1 积分”的完整同位置求值，再重算交点；不应通过调整 van Overstraeten 参数或迁移率补偿。

## 结果文件

- `analysis/branch_closure/branch_current_source_compare.csv`
- `analysis/branch_closure/summary.json`
- `analysis/hotspot_slope_criteria/hotspot_same_bias_compare.csv`
- `analysis/hotspot_slope_criteria/hotspot_slope_compare.csv`
- `analysis/generation_support/threshold_support_compare.csv`
- `analysis/generation_support/radial_cumulative_source_compare.csv`
