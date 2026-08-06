# BVmethods NMOS 载流子专属 Eparallel 与相邻单元插值核对

## 核对目标

固定既有 SG 电子电流、van Overstraeten 参数和雪崩源映射，逐边核对 2663、2666、2672、2675、2681 在 6.4 V 与 7.0 V 的：

1. Vela 载流子专属 `Eparallel` 恢复；
2. 相邻三角形电势梯度到边电场的插值；
3. Sentaurus 节点 `ElectricField`、`eCurrentDensity` 与 `eAlphaAvalanche`；
4. 冻结 6.4 V 电子电流方向后的反事实场增长；
5. 仅改变场恢复算子时对 `alpha_n` 增长的影响。

## 方法

- 完整复现 Vela 当前实现：相邻于目标边的三角形 `grad(psi)` 按面积加权，电子电流由端点所有有效入射 SG 边通量做长度加权最小二乘恢复，再计算 `max(E dot Jn / |Jn|, 0)`。
- 从 Sentaurus 节点 `eAlphaAvalanche` 按当前 300 K van Overstraeten 参数反演电子 `Eparallel`。
- 分别对 Vela 和 Sentaurus 电势执行两种 P1 场恢复：
  - 当前边局部模板：只使用目标边的两个相邻三角形；
  - 节点星形模板：先对每个端点全部 Si 邻接三角形做面积加权，再平均两个端点场。
- 固定 6.4 V 电子电流方向，用 7.0 V 电场重算 `Eparallel`，分离电场与电流方向贡献。

## 自一致性验证

- Python 重建的 Vela `Eparallel` 与原始 `sg_avalanche_edges.csv` 最大相对误差：`2.22e-16`。
- 由重建场代入 van Overstraeten 公式得到的 `alpha_n` 与 Vela CSV 最大相对误差：`4.44e-16`。
- 因此以下差异不是审计脚本的单位、符号或边方向错误。

## 逐边增长分解（6.4 V 到 7.0 V）

| 边 | 当前 Vela Eparallel | Vela 节点星形 Eparallel | Sentaurus 导出 Eparallel | 电流方向单独变化 | 当前 Vela alpha | 节点星形 Vela alpha | Sentaurus alpha |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 2663 | 1.0921 | 1.1258 | 1.1259 | -0.016% | 1.0967 | 1.1639 | 1.1580 |
| 2666 | 1.0938 | 1.1290 | 1.1294 | +0.005% | 1.1008 | 1.1734 | 1.1631 |
| 2672 | 1.0954 | 1.1369 | 1.1386 | +0.009% | 1.1053 | 1.1959 | 1.1672 |
| 2675 | 1.0964 | 1.1307 | 1.1332 | +0.030% | 1.1092 | 1.1834 | 1.1654 |
| 2681 | 1.0974 | 1.1294 | 1.1350 | +0.114% | 1.1136 | 1.1841 | 1.1703 |

## 结论

1. **SG 电子电流方向不是主因。** 固定电场时，电子电流方向变化对五条边的影响仅 `-0.016%` 到 `+0.114%`；固定 6.4 V 电流方向得到的增长与当前 Vela `Eparallel` 增长几乎一致。
2. **van Overstraeten 参数和非线性实现不是主因。** 当前 Vela 场到 `alpha_n` 的映射逐位闭合。
3. **第一处不同的离散算子是电势梯度的空间支撑。** 当前 Vela 只平均目标边两侧单元，因此 6.4 V 时五条边的场比 Sentaurus 节点场高约 22% 到 31%；到 7.0 V 高约 18% 到 25%。这个偏差随偏压缩小，导致当前 Vela 的场增长斜率偏慢。
4. **Sentaurus 导出场可由节点星形 P1 恢复复现。** 对 Sentaurus 电势使用“端点全部 Si 邻接单元面积加权，再平均端点”的恢复，所得 `Eparallel` 与 Sentaurus 导出节点向量的差异仅 `-0.49%` 到 `+0.61%`。
5. **同一节点星形算子也使 Vela 与 Sentaurus 闭合。** Vela 节点星形 `Eparallel` 相对 Sentaurus 导出场的差异为约 `+0.72%` 到 `+1.88%`，且其增长因子从当前的 `1.092-1.097` 修正到 `1.126-1.137`，与 Sentaurus 的 `1.126-1.135` 基本一致。
6. 使用节点星形场预测的 Vela `alpha_n` 增长为 `1.164-1.196`，已消除当前 `1.097-1.114` 的系统性偏慢；与 Sentaurus `1.158-1.170` 尚有 `0.6%-2.9%` 的剩余增长误差，后续需按节点而不是边平均继续核对 `alpha(Eparallel)` 的求值位置。

## 已实现功能

1. 新增 `impact_ionization.eparallel_field_recovery`：
   - `edge_adjacent_cells`：默认值，保留原有边两单元面积加权行为；
   - `nodal_vertex_star`：对端点完整 Si 单元星形模板做面积加权，再平均端点向量。
2. Newton 和 Gummel 配置解析均支持新选项；非法值以及非 `eparallel` 组合会被拒绝。
3. 新模式复用现有 P1 单元梯度、材料筛选和节点 SG 电流恢复，没有引入经验缩放参数。
4. 新增离散回归测试，显式验证默认窄模板为 `(0.5, 0)`，节点星形模板为 `(1/3, 1/6)`。
5. `test_impact_ionization` 全部 51 个用例、679 个断言通过。

## 实际 C++ 重放验证

使用相同 6.4 V、7.0 V 收敛态，只切换 `eparallel_field_recovery=nodal_vertex_star`。实际 C++ 输出与审计预测的最大相对误差为：

- `Eparallel`：`4.44e-16`；
- `alpha_n`：`6.66e-16`；
- `electron_flux_proxy` 相对旧模式变化：严格为 `0`。

| 边 | 新模式 Eparallel 增长 | Sentaurus Eparallel 增长 | 新模式 alpha 增长 | Sentaurus alpha 增长 |
|---:|---:|---:|---:|---:|
| 2663 | 1.12577 | 1.12586 | 1.16393 | 1.15803 |
| 2666 | 1.12905 | 1.12945 | 1.17338 | 1.16313 |
| 2672 | 1.13692 | 1.13863 | 1.19591 | 1.16721 |
| 2675 | 1.13065 | 1.13322 | 1.18345 | 1.16539 |
| 2681 | 1.12937 | 1.13505 | 1.18412 | 1.17035 |

新模式已消除 `Eparallel` 增长系统性偏慢；剩余 `alpha_n` 增长高约 `0.5%-2.5%`，主要来自 Sentaurus 节点 alpha 求值位置与当前边向量平均后求 alpha 的非线性交换误差。下一步应重算完整 6.4--7.0 V 分支，并复查 10%/30% 产生率肩部、累计雪崩源和电流源交点。

## 产物

- `scripts/audit_bvmethods_nmos_eparallel_edge_interpolation.py`
- `edge_eparallel_audit.csv`：逐偏压逐边恢复、Sentaurus 反演与算子比较；
- `adjacent_cell_projection_audit.csv`：相邻单元电场及投影；
- `incident_sg_stencil_audit.csv`：目标端点入射 SG 通量模板；
- `edge_growth_decomposition.csv`：场幅值、电流方向、节点恢复和 alpha 增长分解；
- `implemented_mode_validation.csv`：新 C++ 模式逐边预测闭合；
- `implemented_mode_growth.csv`：新模式与 Sentaurus 增长对比；
- `summary.json`：机器可读摘要。
