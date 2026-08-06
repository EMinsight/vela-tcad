# BVmethods NMOS Band2Band(E2) 实现与闭合记录（2026-08-03）

## 结论

Vela 已实现 Sentaurus `Band2Band(E2)` 局部带间隧穿产生率，并将同一个电子-空穴对产生源装配到 Newton 与 Gummel 连续性方程。使用 Sentaurus O-2018.06 硅参数

\[
G_{BTBT}=A|F|^2\exp(-B/|F|),\quad
A=3.4\times10^{21}\ \mathrm{cm^{-1}s^{-1}V^{-2}},\quad
B=22.6\times10^6\ \mathrm{V/cm}
\]

重建 BVmethods NMOS 分支后，6.4 V 的 Vela 漏电流为 `9.751490e-6 A/um`，Sentaurus 为 `9.485716e-6 A/um`，相对偏差 `+2.80%`。4--6.4 V 的误差随偏压升高由 `+10.77%` 收敛到 `+2.80%`；2 V 为 `-18.59%`，该点绝对电流仅为数 nA/um，仍受低场输运/接触基线差异主导。

## 实现内容

1. 新增 `BandToBandTunnelingModel`，支持 `none`、`e2`/`sentaurus_e2`，内部统一使用 SI；JSON 同时接受 SI 键和 Sentaurus `cm` 参数别名。
2. `RecombinationModelConfig`、`NewtonConfig`、`GummelConfig` 和 DC sweep 配置链均携带 `band_to_band`。
3. Newton 电子、空穴连续性残差减去完全相同的 BTBT 对产生源；Gummel 两个密度方程右端加入相同源。
4. 默认 `frozen_field` 只在线性化中冻结 `dG/dpsi`，每一次非线性残差仍由当前电势重算场和产生率。`potential_finite_difference` 保留完整电势导数，用于小网格核验。
5. BTBT 场和源按半导体三角形计算：单元内求 `|grad(psi)|`，积分 `G*area` 后等权质量集总到三个顶点。Si/SiO2 共享节点不会混入氧化层电场或氧化层控制体积。
6. VTK 输出新增 `Band2BandGeneration`，并修复 DC sweep 调用遗漏 `UnitScalingConfig` 的问题。

## 离散根因定位

最初的节点 WLS 版本在 6.4 V、完整 E2 参数下得到 `1.01623e-4 A/um`，约为 Sentaurus 的 `10.71` 倍；而将参数缩放到 `0.1 A` 时恰好得到约 `1.00718e-5 A/um`。这不是 E2 系数本身需要经验缩放，而是 Si/SiO2 共享节点的 WLS 邻域跨入氧化层：E2 对电场呈指数依赖，界面场污染再乘以混合材料节点体积会显著放大总源。

改成半导体单元积分后，仍使用完整的官方参数，6.4 V 电流变为 `9.75149e-6 A/um`。因此已由离散支撑修复闭合，不需要调整迁移率、接触电压或 E2 参数。

## 多偏压结果

| Drain bias (V) | Sentaurus (A/um) | Vela E2 (A/um) | 相对偏差 | Newton 迭代 | 最大场 (V/m) |
|---:|---:|---:|---:|---:|---:|
| 2.0 | 3.735696e-9 | 3.041395e-9 | -18.59% | 2 | 4.779221e8 |
| 4.0 | 2.673564e-7 | 2.961594e-7 | +10.77% | 6 | 7.606256e8 |
| 5.0 | 1.934757e-6 | 2.043179e-6 | +5.60% | 10 | 8.818396e8 |
| 6.0 | 6.326760e-6 | 6.542639e-6 | +3.41% | 11 | 9.734086e8 |
| 6.4 | 9.485716e-6 | 9.751490e-6 | +2.80% | 15 | 1.008977e9 |

各点均使用完整 `A`、`B` 参数和自洽 BTBT 源。2/4/5 V 从对应无-BTBT 收敛态按 `0.01 A -> 0.1 A -> 1.0 A` 分级激活；6.4 V 使用 `0.001 A -> 0.01 A -> 0.1 A -> 1.0 A`；6.0 V 复用下降分支已接受状态进行同偏压复核。

从 6.0 V 直接跳到 5.0 V 的单次大步下降扫描会触发大量线搜索重试，因此多偏压闭合采用同偏压分级激活。这说明当前 E2 物理闭合已完成，但要把它直接用于连续的大步长 BV 扫描，仍应增加自适应偏压步长或沿 0.05--0.1 V 中间点延续。

## 验证

- 固定 Sentaurus 状态上，以导出节点场代入 E2 公式，1/2/4/5/6/6.4 V 的预测/导出 `Band2BandGeneration` 中位比为 `0.791/0.835/0.874/0.885/0.890/0.892`；剩余差异符合导出节点场与求解器内部单元场恢复差异。
- 4/5/6/6.4 V 的 Sentaurus `q*integral(G_BTBT)` 与漏电流此前已闭合至约 `1.1%/0.4%/0.5%/0.5%`，证明高压漏电流主要由 E2 解释。
- E2 定向测试覆盖公式、SI/cm 单位换算、场导数、电子/空穴成对源、精确 Jacobian 与混合材料排除。
- 完整相关回归：
  - `test_recombination`: 48 assertions / 17 cases
  - `test_impact_ionization`: 644 assertions / 46 cases
  - `test_dd_gummel`: 206 assertions / 12 cases
  - `test_newton_solver`: 1072 assertions / 76 cases
  - `test_dc_sweep`: 2921 assertions / 76 cases

以上测试全部通过，`git diff --check` 通过（仅报告现有 LF/CRLF 提示）。

## 输出位置

- 修正后的分级激活和状态：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_semiconductor_cell_20260803/`
- 初始节点场公式对比：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rebuild_fd_gummel_20260803/btbt_e2_fixed_state_compare.csv`
- 可复现运行脚本：`scripts/run_bvmethods_nmos_btbt_e2_validation.py`
- 固定状态公式检查脚本：`scripts/compare_sentaurus_e2_btbt.py`

## 剩余工作

1. 将自适应偏压延续接入完整 0--7 V BTBT 分支，避免 6 -> 5 V 大跳步线搜索退化。
2. 对 2 V 的低电流差异继续分离接触/SG 输运基线与 BTBT 贡献；不应通过调整 E2 参数补偿。
3. 在 E2 分支基础上重新计算 IIC 积分及雪崩后处理，确认 BTBT 改变后的载流子电流空间支撑能否关闭剩余 IIC 差异。
