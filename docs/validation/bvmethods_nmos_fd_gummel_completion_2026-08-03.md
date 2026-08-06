# BVmethods NMOS：Fermi–Dirac Gummel 密度块与 6.4 V 分支推进

日期：2026-08-03

## 结论

本轮完成了 Gummel 密度求解路径的 Fermi–Dirac 广义 Einstein/Scharfetter–Gummel
离散，并把相同统计模型接入欧姆接触、初始态、准费米势反演、Poisson 载流子导数和
Newton 载流子行恢复路径。此前在 `5.725 V` 后因 Fermi–Dirac Gummel 恢复不可用而
中断的 `postprocess_only` 基本 DD 分支已经越过：`5.73125 V` 使用
Gummel-to-Newton 交接收敛，之后纯 Newton 连续推进到 `6.4 V`。

这意味着本轮目标已经完成，但不等于 IIC 基准已经闭合。在 `6.4 V`，电场峰值和
electron alpha 峰值已经接近 Sentaurus；端电流、电子/空穴 SG 电流空间支撑和雪崩源
积分仍有数量级差异，是下一阶段的主要闭合对象。

## 实现内容

### 1. 密度形式的广义 SG 离散

`DDAssembler` 现在持有每节点 `Nc`、`Nv` 和载流子统计配置。对每条半导体输运边：

- 由旧密度反演两端约化费米能级 `eta`；
- 用两端密度和 `eta` 冻结广义 Einstein 因子；
- 电子采用与耦合 DD 相同的漂移势：
  `dpsi + Vt*log((ni1/Nc1)/(ni0/Nc0))`；
- 空穴采用：
  `dpsi + Vt*log((ni0/Nv0)/(ni1/Nv1))`；
- Bernoulli 权重使用 `Vt*EinsteinFactor`，矩阵通量系数同步乘该因子；
- 极低密度下若 Fermi–Dirac 导数下溢，使用非简并极限 `dn/deta -> n`，避免零导数。

Boltzmann 路径保持原有离散，不改变其数值行为。

### 2. Fermi–Dirac 状态语义统一

- Gummel 不再拒绝 `carrier_statistics=fermi_dirac`；
- `effectiveNi` 统一使用共享构造与材料归属规则；
- 欧姆接触通过统一的平衡载流子状态计算 `psi/n/p`；
- 初始密度和每轮 `phin/phip` 反演均使用 `Nc/Nv` 与相同统计模型；
- Poisson 密度导数、SRH/Auger 平衡载流子乘积与准费米势驱动迁移率使用同一统计定义；
- Newton 的 Gummel 密度恢复显式传递 Fermi–Dirac 配置。

### 3. 测试覆盖

新增两条关键回归：

- 完整 Fermi–Dirac Gummel 零偏平衡态收敛，并验证密度/准费米势正反变换一致；
- Fermi–Dirac Gummel 密度恢复能够修复 Newton 的退化载流子行。

全量工程编译通过。直接执行的相关测试结果：

| 测试程序 | 结果 |
|---|---:|
| `test_carrier_statistics` | 22 assertions / 5 cases，通过 |
| `test_sg_flux` | 208 assertions / 25 cases，通过 |
| `test_dd_gummel` | 202 assertions / 11 cases，通过 |
| `test_newton_solver` | 1068 assertions / 75 cases，通过 |
| `test_mos_mixed_material` | 659 assertions / 5 cases，通过 |
| `test_impact_ionization` | 613 assertions / 44 cases，通过 |

## 从 5.725 V 恢复到 6.4 V

### 断点跨越

以原 `5.725 V` 接受态为初值，在 `5.73125 V` 启用 Gummel-to-Newton：

- Gummel 密度阶段运行 80 次后交给 Newton；
- Newton 4 次迭代后以 `carrier_row_qualified_residual_floor` 接受；
- 载流子行违规数为 0，最大合格残差比为 `9.032816e-4`；
- QF 边界违规数为 0；
- `Id = 4.999969e-9 A/um`，最大电场 `9.488176e6 V/cm`。

随后 `5.75 V` 由纯 Newton 在 7 次迭代内收敛。再以 `0.05 V` 步长连续运行
`5.80, 5.85, ..., 6.40 V`，13 个检查点全部接受，所有点的载流子行违规和 QF
违规均为 0。

### 6.4 V 接受态

| 量 | Vela 结果 |
|---|---:|
| Newton 迭代数 | 20 |
| 载流子行最大合格残差比 | 5.555744e-4 |
| 载流子行违规数 | 0 |
| QF 违规数 | 0 |
| drain current | 5.399425e-9 A/um |
| 最大电场 | 1.007593e7 V/cm |
| 最大 electron alpha | 4.602750e7 1/m |
| 最大 hole alpha | 6.436261e7 1/m |
| 雪崩边源原始积分 | 1.856332e22 |
| 换算的积分雪崩电流 | 约 2.974e-9 A/um |

SG 端电流与连续性残差端电流一致到浮点舍入误差：drain 分别为
`5.3994188360868002e-9` 和 `5.3994188360867969e-9 A/um`。这说明 Vela 内部端电流
提取路径已经自洽；它与 Sentaurus 的差异来自物理解/离散支撑，而不是两个 Vela
端电流提取器互相矛盾。

## 6.4 V 与 Sentaurus 的严格同网格比较

| 物理量 | Sentaurus 峰值 | Vela 峰值 | Vela/Sentaurus | 空间相关系数 |
|---|---:|---:|---:|---:|
| electric field (V/m) | 2.283188e8 | 2.450294e8 | 1.07319 | 0.99700 |
| electron alpha (1/m) | 4.404434e7 | 4.602750e7 | 1.04503 | 0.58549 |
| hole alpha (1/m) | 3.502760e7 | 6.436261e7 | 1.83748 | 0.60287 |
| electron current density (A/m2) | 1.395069e9 | 1.400273e7 | 1.00373e-2 | 0.00293 |
| hole current density (A/m2) | 5.073870e10 | 6.014482e-4 | 1.18538e-14 | 0.92045 |
| avalanche generation (1/m3/s) | 5.700771e36 | 1.605402e32 | 2.81611e-5 | -0.03696 |

载流子状态比较表明实现后分支已经明显改善，但高分位误差仍大：

- electron QF：中位绝对误差 `9.30e-5 V`，95 分位 `0.556 V`；
- hole QF：中位绝对误差 `0.0984 V`，95 分位 `0.819 V`；
- electron density：中位 `0.217 dex`，95 分位 `9.552 dex`；
- hole density：中位 `2.395 dex`，95 分位 `13.902 dex`。

Sentaurus 在 `6.4 V` 的 drain current 为 `9.485716e-6 A/um`，积分雪崩电流约
`8.896588e-6 A/um`。本轮运行是无自洽雪崩反馈的 `postprocess_only` 基本 DD 分支，
其端电流不能直接作为最终 IIC 曲线；但由于 electron current density 峰值也只有
Sentaurus 的约 1%，且空间相关性很低，差异不能仅归因于是否开启雪崩反馈。

## 后续闭合顺序

1. 固定 6.4 V 同网格状态，定位 Sentaurus 与 Vela 电子电流峰值边、主导电流路径和
   接触邻边，逐边核对 QF 差、广义 Einstein 因子、Bernoulli 参数、迁移率和面积权重。
2. 对空穴连续性单独审计。当前空穴端电流和雪崩区电流密度过小，是 hole-alpha
   虽接近但空穴雪崩贡献无法闭合的直接原因。
3. 在基本 SG 电流空间支撑闭合后，再核对 `alpha*|J|/q` 的单元/边映射、体积积分和
   IIC 电流提取，避免用经验缩放 alpha 掩盖输运差异。
4. 最后开启自洽雪崩，重建高场分支并验证外接电阻和 voltage-to-current 两种 BV
   路径；在基本 DD/IIC 未闭合前不调整迁移率或接触电压拟合端电流。

## 可复现产物

- 5.73125 V Gummel-to-Newton 状态：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rebuild_fd_gummel_20260803/step_5p73125`
- 5.8–6.4 V 连续分支：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rebuild_fd_gummel_20260803/trunk_5p8_6p4_newton`
- 6.4 V 完整 SG/雪崩诊断：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rebuild_fd_gummel_20260803/probe_6p4_full`
- 6.4 V Sentaurus/Vela 同网格比较：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/iic_rebuild_fd_gummel_20260803/analysis_6p4`
