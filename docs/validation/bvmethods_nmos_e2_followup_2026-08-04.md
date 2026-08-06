# BVmethods NMOS E2 全分支与 IIC 复核（2026-08-04）

## 结论

本轮完成了自洽 `Band2Band(E2)` 模型的 `0--7 V` 连续反偏分支，共 71 个
目标偏压点，全部收敛且没有 QF bounds 告警。E2 已将 4--7 V 的端电流闭合到
Sentaurus 的约 `2.2%--10.8%`，但 IIC 尚未闭合：在 Sentaurus 官方稀疏提取
电压 `6.377494278 V`，Vela 端电流只高 `5.67%`，积分雪崩电流却只有
Sentaurus 的 `49.88%`。

因此，当前主要障碍已经从“高压分支无法推进/缺少 E2 电流”收敛为“雪崩源的
空间支撑与积分约低一半”。不能通过继续调整 E2 系数或端电流提取器解决该问题。

## 1. 0--7 V 自适应 E2 分支

运行采用 0.1 V 目标点、内部自适应缩步和 0.5 V 可恢复分段。为避免偏压发生变化
后复用旧状态，绝对残差阈值设为 `5e-10`；实际 0.1 V 偏压变化产生的初始残差约
为 `1e4`，有效 Newton 更新后才进入 `1e-10` 数值噪声区。载流子局部行判据在
主分支使用 `eps_row=2e-3`，IIC 密集探针使用 `eps_row=1e-3`。

| Vd (V) | Vela Id (A/um) | Sentaurus Id (A/um) | 相对误差 | max(E) (V/m) |
|---:|---:|---:|---:|---:|
| 1.0 | 2.378272e-9 | 2.993217e-9 | -20.54% | 2.838119e8 |
| 2.0 | 3.041395e-9 | 3.735698e-9 | -18.59% | 4.779221e8 |
| 4.0 | 2.961594e-7 | 2.673564e-7 | +10.77% | 7.606256e8 |
| 5.0 | 2.043179e-6 | 1.934757e-6 | +5.60% | 8.818396e8 |
| 6.0 | 6.542639e-6 | 6.326760e-6 | +3.41% | 9.734086e8 |
| 6.4 | 9.751490e-6 | 9.485716e-6 | +2.80% | 1.008977e9 |
| 7.0 | 1.690718e-5 | 1.654544e-5 | +2.19% | 1.062330e9 |

数值质量汇总：

- 目标点数：`71`；未收敛点：`0`；
- QF bounds 告警总数：`0`；
- 主分支最大载流子行比值：`1.975513e-3`，低于 `2e-3` 门槛；
- 早期区段在执行效率修正前最多运行 200 次 Newton；修正后关键高压点通常为
  3--10 次目标点迭代，困难区间由内部自适应子步穿越；
- 分段首点不再重复求解上一分段的终点，已完成分段可直接跳过并从断点恢复。

## 2. 2 V 低电流基线拆分

| 量 | 数值 (A/um) |
|---|---:|
| Vela 无 E2 电流 | 2.909320e-9 |
| Vela 自洽 E2 电流 | 3.041395e-9 |
| Sentaurus 电流 | 3.735698e-9 |
| E2 带来的电流增量 | 1.320749e-10 |
| 加入 E2 后仍剩余的差值 | 6.943026e-10 |

E2 增量只占完整 Vela 2 V 电流的 `4.34%`，而加入 E2 后相对 Sentaurus 仍低
`18.59%`。同时，最大电场从 `4.779220323e8 V/m` 变为
`4.779220555e8 V/m`，相对变化仅 `4.86e-8`。这说明 2 V 误差主要属于欧姆接触/
广义 SG 输运基线，而不是 E2 产生率或静电场。

## 3. E2 状态上的 IIC 密集复核

从连续分支的 6.3 V 状态出发，重新运行 17 个带逐边雪崩诊断的检查点，覆盖
`6.30--7.00 V`，并显式包含 `6.377494278 V`。

| Vd (V) | Vela Id (A/um) | Vela Iava (A/um) | Vela Iava/Id | Sentaurus Iava/Id |
|---:|---:|---:|---:|---:|
| 6.377494 | 9.541683e-6 | 4.226083e-6 | 0.442908 | 0.938353 |
| 6.50 | 1.073032e-5 | 4.856270e-6 | 0.452575 | 0.955697 |
| 6.70 | 1.293065e-5 | 6.077146e-6 | 0.469980 | 0.992960 |
| 6.80 | 1.416131e-5 | 6.788044e-6 | 0.479337 | 1.012233 |
| 7.00 | 1.690718e-5 | 8.437945e-6 | 0.499075 | 1.051626 |

在 `6.377494278 V`：

- Vela/Sentaurus 端电流比为 `1.056724`；
- Vela/Sentaurus 积分雪崩电流比为 `0.498779`；
- Vela 最大电子/空穴雪崩系数分别为 `4.208491e7` 和 `3.209271e7 1/m`；
- 截至 7 V，Vela 最大 `Iava/|Id|` 仍只有 `0.499075`，没有包围电流交点。

Sentaurus 精确密集检查点的 `Iava-Id=0` 交点位于
`6.734425890 V`（6.7--6.8 V 线性插值）。这与官方 Workbench 稀疏曲线/
`BreakAtIonIntegral` 提取的 `6.377494278 V` 是两种不同提取口径，必须继续分别保留；
Vela 当前两种口径都尚未闭合。

## 4. 根因判断与后续顺序

1. **先闭合雪崩源积分。** 在 6.377 V 端电流已接近的前提下，逐边拆分电子/
   空穴 SG 通量、alpha、边长/对偶面积、单元到节点/边的源映射和半导体区域筛选，
   定位约 `0.50x` 积分因子的来源；不得先统一缩放 van Overstraeten 系数。
2. **实现或复核路径离化积分。** 对齐 Sentaurus 的 `eIonIntegral`、`hIonIntegral`、
   `MeanIonIntegral` 和 `BreakAtIonIntegral` 停止/提取语义，避免把体积分雪崩电流交点
   与官方 IIC 路径积分阈值混为同一指标。
3. **继续 2 V 输运审计。** 固定同一异常接触/有源边，核对接触到内部的准费米势
   降、Fermi--Dirac 广义 Einstein 因子、Bernoulli 参数、迁移率与几何权重。
4. **IIC 闭合后再进入外电路方法。** 外接电阻和 voltage-to-current 均依赖同一高压
   DD/雪崩状态；在源积分仍低约一半时实现只会把物理误差带入外电路分支。

## 5. 可复现产物

- 主分支脚本：`scripts/run_bvmethods_nmos_e2_followup.py`
- 汇总脚本：`scripts/analyze_bvmethods_nmos_e2_followup.py`
- 71 点主分支：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_adaptive_0_7_20260804/branch_0_7.csv`
- 17 点 IIC 逐边诊断：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_reclosure_20260804/postprocess_only`
- 电流分支对比：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_followup_20260804/analysis/e2_branch_compare.csv`
- 2 V 拆分：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_followup_20260804/analysis/e2_2v_baseline.csv`
- IIC 对比：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_followup_20260804/analysis/e2_iic_compare.csv`
- 机器可读摘要：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_followup_20260804/analysis/e2_followup_summary.json`
