# TransportModels 分区准费米参考实现与验证（2026-08-22）

## 结论

已实现可选的 `newton.quasi_fermi_reference = "contact_basin"` 模式。电子准费米势以 n 型源/漏欧姆接触为种子，空穴准费米势以 p 型欧姆接触为种子；各输运节点按材料连通图上离最近种子的距离选择参考。当前 MOS 算例因此形成源区 `0 V`、漏区 `1.1 V` 的电子准费米参考，而栅极金属不参与输运分区。

该改动解决了本轮针对的问题：源区约 `6.81e-17 V` 的 Newton 电子准费米更新不再因写回 `1.1 V + delta` 而全部舍入为零。深关断冻结态的源端电流由数值零恢复为约 `-1.98006e-16 A/um`，四端 KCL 残差仅为漏电流的约 `0.0221%`。

完整严格 Newton 求解尚未达到最终收敛判据：在 7 次已接受更新后，第 8 次更新因线搜索不下降停止。此时电子、空穴连续性块已在 `1e-16` 量级，剩余范数由 Poisson 块 `1.47685e-11` 主导，并仍有 169 个载流子行违反局部门槛。故本次工作确认“局部参考＋小增量”问题已经解决，但不能据此宣称整个深关断工作点已经严格收敛。

## 实现范围

- `CoupledDDAssembler` 支持逐节点电子/空穴准费米参考场；标量参考仍保留为兼容回退。
- 解向量打包、解包、载流子密度、SRH、边 SG 通量、迁移率驱动力、解析 Jacobian、规范条件和边界条件均改为读取节点局部参考。
- 跨分区边在 node-0 锚点下用 `long double` 计算参考变换和准费米差，避免先重构两个大绝对数再相减。
- Newton 重启文件保存每个节点的 `reference` 和 `increment`；读入时用 `old_reference + old_increment - new_reference` 重分区。
- 接触电流和接触 SG 诊断直接消费逐节点参考表示。
- `contact_basin` 通过输运材料邻接图上的多源最短路分配参考；不可达节点使用最近接触的几何回退。

## 亚 ULP 更新 A/B

输入为同一 `Vg = -0.68 V` 失败状态，关闭准费米更新上限，比较全局/分区坐标。

| 模式 | 源区非零原始步行数 | 写回后丢失行数 | 最大原始步 | 最大实际写回步 | 线搜索 |
|---|---:|---:|---:|---:|---|
| `contact_majority`（全局 `1.1 V`） | 25 | 25 | `6.80876e-17 V` | `0 V` | 未启用 |
| `contact_basin` | 25 | 0 | `6.80876e-17 V` | `6.80876e-17 V` | 未启用 |
| `contact_basin` | 25 | 0 | `6.80876e-17 V` | `2.12774e-18 V` | 接受，6 次尝试，阻尼 `0.03125` |

全局坐标下源区更新只有局部 ULP 的约 `0.14%` 到 `30.7%`，写回物理绝对势后全部消失。分区坐标以源区 `0 V` 为参考，增量本身就是存储值，因此 25 行全部保留。

原始表：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/source_newton_linear_audit_20260822/summary.csv`。

## 分区与重启检查

冻结状态共 3315 个节点，电子参考分布为：

| 电子参考 | 节点数 | 典型区域 |
|---:|---:|---|
| `0 V` | 1524 | 源区盆地及其连通邻域 |
| `1.1 V` | 1576 | 漏区盆地及其连通邻域 |
| `-0.68 V` | 215 | 栅极/非输运边界状态 |

节点 2506 的电子参考为 `0 V`，电子准费米增量为约 `6.80854e-17 V`；漏区节点 825 的参考为 `1.1 V`。这证明重启文件没有再次把源区小增量合并进漏端绝对参考。

状态文件：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/contact_basin_full_noglobal_20260822/rejected_states/attempt_1_bias_m0p680000_final.csv`。

## 深关断端口与 SRH 守恒

对最后一个已接受状态作冻结重放，单位均为 `A/um`：

| 量 | 数值 |
|---|---:|
| 源端总电流 | `-1.98005910e-16` |
| 漏端总电流 | `2.03319036e-16` |
| 栅端总电流 | `3.89440591e-22` |
| 衬底总电流 | `-5.35852076e-18` |
| 四端 KCL 残差 | `-4.50047386e-20` |
| `|Id| / |KCL|` | `4517.73` |

KCL 残差约为漏电流的 `2.21e-4`，远小于“至少低一个数量级”的深关断解析门槛，冻结态可判定为端口电流已解析。

硅区 SRH 净产生等效电流为 `5.34851895e-18 A/um`，衬底空穴电流为 `5.35014271e-18 A/um`，幅值相对差约 `0.0304%`。这说明源端、漏端及衬底的极低电流已经按连续性方程闭合，而非由单个端口的消减噪声构成。

原始数据：

- `build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/contact_basin_full_noglobal_20260822/frozen_replay/terminal_balance.csv`
- `build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/contact_basin_full_noglobal_20260822/frozen_replay/srh_balance.csv`

## 严格 Newton 的剩余问题

完整求解最后记录：

| 指标 | 数值 |
|---|---:|
| 已接受 Newton 迭代 | 7 |
| 最终组合残差 | `1.47685332e-11` |
| Poisson 块 | `1.47685332e-11` |
| 电子连续性块 | `4.60104118e-17` |
| 空穴连续性块 | `9.81681214e-17` |
| 载流子违规行 | 169 |
| 最大违规比 | `1.64030` |
| 停止原因 | `line_search_non_decrease` |

因此下一步不应继续修改准费米存储坐标，而应审计 Poisson 主导方向的线搜索目标、违规载流子行与 Poisson 更新之间的耦合，并决定是否采用分块 merit function 或对近收敛 Poisson 步使用更合适的尺度。

## 测试

最终版本重新编译后通过：

- `test_newton_solver.exe`：1218 assertions，91 cases。
- `test_dc_sweep.exe`：3471 assertions，98 cases。
- `test_mos_mixed_material.exe`：1437 assertions，9 cases。

新增回归覆盖：

- 全局参考和节点参考下密度、残差、Jacobian 的不变性。
- 源区 `5e-18 V` 增量在全局 `1.1 V` 坐标下丢失、在源区 `0 V` 局部坐标下保留。
- `contact_basin` 配置解析。

## 复现入口

- 亚 ULP A/B：`scripts/run_transportmodels_source_newton_linear_audit.py`
- 完整/冻结态验证：`scripts/run_transportmodels_contact_basin_validation.py`
