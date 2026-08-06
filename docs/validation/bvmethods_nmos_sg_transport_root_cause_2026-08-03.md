# BVmethods NMOS 6.4 V SG 输运逐边审计与根因结论

日期：2026-08-03

## 结论

固定 6.4 V 后，对 Sentaurus 最大电子电流边 2226（节点 344--1335）逐项核对。Vela 的电子电流密度仅为 Sentaurus 的 `3.674714328e-11`，而同一边上的准费米势差、迁移率、广义 Einstein 因子和几何权重均接近一致。唯一与电流缺口同阶的量是电子密度，其 Vela/Sentaurus 比为 `3.754547112e-11`。

因此，高压分支的第一处主导分叉不是 SG Bernoulli 公式、迁移率或几何权重，而是绝对准费米势位置偏移造成的载流子人口不足。进一步积分 Sentaurus 的 `Band2BandGeneration` 后确认：4 V 以上的漏电流几乎全部由当前 Vela 配置缺失的 `Band2Band(E2)` 产生项解释。

## 6.4 V 热点边分量审计

| 分量 | Vela/Sentaurus |
|---|---:|
| 电子电流密度 | 3.674714328e-11 |
| 电子准费米势边差 | 9.703256516e-01 |
| 电子密度对数均值 | 3.754547112e-11 |
| 电子迁移率 | 1.099862166 |
| 广义 Einstein 因子 | 9.739154933e-01 |
| 几何权重 | 1.000000000 |

Sentaurus 与 Vela 的广义 Bernoulli 参数分别为 `15.3862` 和 `15.5195`。这一差异不足以解释约 10 个数量级的电子电流密度差；电流比与载流子密度比的一致性直接指向准费米势绝对位置/产生项。

## 随偏压出现分叉的位置

| Vd (V) | Sentaurus Id (A/um) | Vela Id (A/um) | Vela/Sentaurus | Sentaurus 积分 BTBT (A/um) | Vela-Sentaurus 电子 QF 均值偏移 (V) |
|---:|---:|---:|---:|---:|---:|
| 1.0 | 2.993217e-9 | 2.378269e-9 | 7.945529e-1 | 6.535077e-17 | 3.654420e-7 |
| 2.0 | 3.735696e-9 | 2.909320e-9 | 7.787893e-1 | 7.383209e-11 | 3.514087e-2 |
| 4.0 | 2.673564e-7 | 3.985787e-9 | 1.490814e-2 | 2.642992e-7 | 4.326829e-1 |
| 5.0 | 1.934757e-6 | 4.590063e-9 | 2.372424e-3 | 1.942542e-6 | 5.309422e-1 |
| 6.0 | 6.326760e-6 | 5.399425e-9 | 8.534266e-4 | 6.360026e-6 | 5.886145e-1 |
| 6.4 | 9.485716e-6 | 5.399419e-9 | 5.692158e-4 | 9.535500e-6 | 6.082959e-1 |

`Band2BandGeneration` 使用原始半导体三角形面积、1 um 器件深度和电子电荷积分。4/5/6/6.4 V 的积分结果分别与 Sentaurus 漏电流吻合至约 1.1%、0.4%、0.5% 和 0.5%。这构成了缺失 BTBT 是主因的闭合证据。

## 官方算例物理模型核对

本地保存的官方 `pp4_des.cmd` 明确启用：

```text
Mobility(DopingDep HighFieldSaturation(GradQuasiFermi) Enormal)
Recombination(SRH(DopingDep) Band2Band(E2) Avalanche(Eparallel))
Fermi
AvalPostProcessing
```

原始文件：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/raw/pp4_des.cmd`。Sentaurus 2018 用户手册位于虚拟机：`/usr/synopsys/sentaurus/O_2018.06-SP2/tcad/O-2018.06-SP2/manuals/PDFManual/data/sdevice_ug.pdf`。

## 本轮代码修复

发现并修复了一个独立的一致性缺陷：连续性方程和接触电流已使用 Fermi--Dirac 广义 SG，但 `sg_avalanche_edges` 与 `currentDensityAvalancheSourceIntegrals` 仍使用旧的 Boltzmann/可变 `ni` 重构。现已：

1. 将载流子统计、Nc/Nv 传入雪崩后处理电流重构；
2. 在 Fermi--Dirac 模式下使用广义 Einstein 因子与广义 Bernoulli 参数；
3. 保留 Boltzmann 路径兼容性；
4. 在逐边 CSV 中输出是否使用 Fermi--Dirac、广义 Einstein 因子和广义 Bernoulli 参数；
5. 新增直接与广义 SG 通量对照的退化态单元测试。

修复后 6.4 V 漏端提取电流与连续性残差仍一致：`5.3994188360868002e-9` 与 `5.3994188360867969e-9 A/um`。雪崩源积分只发生微小变化，说明雪崩主区域接近非简并；该修复消除了诊断公式不一致，但不会也不应人为填补 BTBT 物理缺口。

## 产物

- 分析脚本：`scripts/audit_bvmethods_nmos_sg_transport.py`
- 完整逐边表：`transport_audit_6p4/edge_transport_audit.csv`
- Sentaurus 顶部电流边：`transport_audit_6p4/top_sentaurus_electron_current_edges.csv`
- 分量统计：`transport_audit_6p4/component_ratio_summary.csv`
- 多偏压演化：`transport_audit_6p4/hotspot_bias_evolution.csv`

产物根目录：`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rebuild_fd_gummel_20260803/transport_audit_6p4`。

## 后续实施顺序

1. 实现可配置的 Sentaurus 兼容 `Band2Band(E2)` 局部产生率，包括硅默认参数、单位换算和电场驱动定义；
2. 将同一对产生项以一致符号装配到电子、空穴连续性方程，并实现 Newton/Gummel 所需线性化或受控冻结策略；
3. 输出节点/单元 BTBT 产生率及面积积分电流，先在固定 Sentaurus 状态下闭合 2/4/5/6/6.4 V；
4. 再重建 Vela 自洽分支，检查准费米势绝对位置、载流子密度和漏电流；
5. BTBT 闭合后再量化 `SRH(DopingDep)` 与 `Enormal` 对低压和表面通道的次级影响；
6. 最后恢复 IIC 雪崩积分判据，验证击穿电压，而不是用迁移率或经验接触电压拟合端电流。
