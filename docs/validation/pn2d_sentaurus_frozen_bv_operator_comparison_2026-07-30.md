# PN2D Sentaurus 冻结状态 BV 算子对比

日期：2026-07-30

## 目的

在完全不推进状态、不向连续性方程写回碰撞电离源项的条件下，将
Sentaurus 的电势、电子/空穴准费米势和电子/空穴浓度导入 Vela，
使用当前生产 BV 算子依次重建驱动力、迁移率、电流、电离系数和源项，
定位两套实现的第一个显著差异。

## 执行合同

- 网格：27 节点、32 个三角形的同拓扑 coarse_0p05 网格。
- 偏压：`-19.7 V`、`-19.8 V`。
- 导入状态：`psi`、`phin`、`phip`、`n`、`p`。
- Vela 碰撞电离配置：生产 `van_overstraeten`、
  `quasi_fermi_gradient/cell_gradient`、
  `cell_reconstructed`、`triangle_gss_gradqf_truncated`。
- 强制设置：`coupling_mode=postprocess_only`。
- Sentaurus 节点量到 Vela 边/单元支撑的比较采用端点平均；总
  `ImpactIonization` 采用线性三角形节点积分。

两组状态的节点 ID 和坐标完全一致，最大坐标误差为 `0 um`。每个偏压
生成 192 条生产支撑过程记录；`solver_coupled` 条数为 0，非零连续性
残差反馈条数为 0。

## 结果

| 指标 | -19.7 V | -19.8 V |
|---|---:|---:|
| 五个导入状态场回读误差 | 0 | 0 |
| QFP 梯度独立重建相对 L2（e/h） | 1.36e-16 / 1.66e-16 | 1.35e-16 / 1.51e-16 |
| 迁移率中位相对误差（e/h） | 8.05% / 6.64% | 8.63% / 7.25% |
| 迁移率峰值比（e/h） | 1.95 / 1.47 | 1.95 / 1.47 |
| 电流幅值相对 L2（e/h） | 2.20e7 / 1.64e7 | 7.02e6 / 5.21e6 |
| 电离系数相对 L2（e/h） | 0.505 / 0.512 | 0.504 / 0.510 |
| 电离系数峰值比（e/h） | 1.071 / 1.116 | 1.066 / 1.108 |
| Vela 总源积分 / Sentaurus 总源积分 | 2.4855e7 | 7.9271e6 |

Vela 总源积分从 `-19.7 V` 的 `2.4970e16 1/(m*s)` 降至 `-19.8 V`
的 `9.6922e15 1/(m*s)`；Sentaurus 则从 `1.0046e9 1/(m*s)` 升至
`1.2227e9 1/(m*s)`。两者不仅绝对量不同，偏压趋势也相反。

`-19.7 V` 的代表性热点是 cell 12、edge 15→12：

- Vela GSS midpoint electron density：`1.1447e18 m^-3`；
- QFP 梯度：`3.7949e7 V/m`；
- Vela 重建电子电流幅值：约 `1.9276e4 A/m^2`；
- Sentaurus 端点平均电子电流幅值：约 `2.491e-4 A/m^2`。

该支撑上的电流相差约 `7.7e7` 倍，并直接进入
`G=(alpha_n*|Jn|+alpha_p*|Jp|)/q`。

## 结论

1. Sentaurus 状态导入和 Vela 状态到 QFP 梯度的重建闭合到机器精度，
   因而此次巨大源项差异不是 CSV 导入、单位换算或状态被求解器修改造成的。
2. 按 5% 阈值，第一个统计偏离出现在迁移率，但约 7%--9% 的典型迁移率
   偏差不足以解释百万至千万倍的源项偏差。
3. 主导差异出现在 Vela `cell_reconstructed` 电流代理。它在跨结区支撑上
   将 GSS midpoint density 与单元 QFP 梯度组合，得到远大于 Sentaurus
   导出电流的电流幅值；该误差随后被 alpha 乘入局部生成率和积分源项。
4. `alpha_n/alpha_p` 的峰值仅相差约 7%--12%，因此现有证据不支持把
   Van Overstraeten 参数本身列为首要根因；其空间分布仍有约 50% 的
   L2 差异，需要在修正电流输入后继续检查。

## 外部 Jn/Jp 只读替换结果

保持 Vela alpha、生产 `triangle_gss_gradqf_truncated` 几何权重和
`postprocess_only` 不变，仅将生成率中的电流幅值替换为 Sentaurus
`Jn/Jp` 端点平均矢量幅值：

| 反事实源积分 / Sentaurus P1 面积分 | -19.7 V | -19.8 V |
|---|---:|---:|
| Vela 原始电流代理 + Vela alpha | 2.4855e7 | 7.9271e6 |
| Sentaurus Jn/Jp + Vela alpha | **1.03982** | **1.03762** |
| Sentaurus Jn/Jp + Sentaurus alpha + Vela 几何 | 1.04925 | 1.04675 |
| Sentaurus ImpactIonization + Vela 几何 | **1.00000** | **1.00000** |

替换电流后，积分误差从百万至千万倍降至 `3.76%--3.98%`。直接把
Sentaurus `ImpactIonization` 投影到 Vela 源项支撑则与线性三角形面积
积分完全闭合，说明当前生产几何权重不是积分总量差异的来源。

Sentaurus 节点导出的 `alpha*|J|/q` 与 `ImpactIonization` 的相对 L2
闭合约为 `8.6%`；投影到 Vela 局部支撑后，Sentaurus alpha/current
组合的局部 L2 约为 `11.7%`。这给出了节点场导出和端点投影本身的误差
底限。使用 Vela alpha 时局部 L2 仍约 `59%`，但其正负空间偏差在总
积分上大幅抵消，最终积分只高约 4%。

该反事实确认：原始百万至千万倍差异由 Vela 电流代理主导；alpha 的
空间离散差异是后续次级问题，源项几何总量闭合。

## 电流代理因子分解

生产 triangle-GSS 电流代理的精确闭合式是：

```text
J_proxy = q * mobility * triangle_gss_midpoint_density * edge_QFP_drive
```

独立重算与过程记录闭合到 `1.5e-16` 以下。注意 alpha 使用单元 QFP
梯度，而电流乘积使用边 QFP 压降；两者在当前正交热点边上数值相同，
但语义支撑不同。

| 控制量 | -19.7 V | -19.8 V |
|---|---:|---:|
| 原始 SG 边电流 vs Sentaurus 边投影 L2（e/h） | 8.42% / 8.18% | 8.09% / 7.86% |
| SG/Laux 重建矢量电流中位误差（e/h） | 1.50% / 0.64% | 1.54% / 1.08% |
| SG/Laux 矢量源积分 / Sentaurus | **1.00954** | **1.00948** |
| triangle-GSS midpoint / SG 边 midpoint 最大比 | 7.253e7 | 2.335e7 |
| 生产代理电流 / 原始 SG 电流最大比 | 6.923e7 | 2.236e7 |

最大热点在两个偏压下均为 electron、cell 17、edge 35：

| 物理量 | -19.7 V | -19.8 V |
|---|---:|---:|
| triangle-GSS midpoint density (m^-3) | 1.1447e18 | 4.3523e17 |
| SG 边 midpoint density (m^-3) | 1.5783e10 | 1.8639e10 |
| Sentaurus 电流反推 density (m^-3) | 1.4790e10 | 1.7576e10 |
| production proxy current (A/m2) | 1.9276e4 | 7.3297e3 |
| Vela raw SG current (A/m2) | 2.7843e-4 | 3.2781e-4 |
| Sentaurus current (A/m2) | 2.4904e-4 | 2.9600e-4 |

在热点上，SG 边 midpoint 和 Sentaurus 电流反推密度均接近低浓度端，
而 triangle-GSS midpoint 选择了高浓度端。保持同一 mobility 和边 QFP
drive，只改用 SG 边 midpoint，`-19.7 V` 热点电流为
`2.6577e-4 A/m2`，已接近 Sentaurus 的 `2.4904e-4 A/m2`。

将 mobility 单独替换成较大的 Sentaurus 端点平均值会进一步放大代理
电流；算术、几何或对数 midpoint 也不能恢复 SG 的漂移/扩散抵消。
因此首要差异不是 mobility 或 QFP drive，而是
`gss_logistic` midpoint 在 triangle source current 中的载流子方向/
语义选择。原始 SG 边电流和由其重建的二维矢量均与 Sentaurus 接近。

八个最大 triangle-GSS 过程记录贡献了两个偏压总源项的
`99.99999%` 以上，所以该局部 midpoint 选择足以解释总曲线异常。

## 三角形单元总杂质重构对电流的影响

`cell_reconstructed_total_impurity` 不是对 Poisson 掺杂或载流子状态做
平滑。它只把每个三角形三个节点的 donor+acceptor 总杂质浓度取平均，
作为该单元 Masetti 迁移率模型的掺杂输入。迁移率随后进入 SG 输运电流，
也进入碰撞电离电流代理，所以它会间接但真实地影响电流。

既有正向 IV 对比给出的定量影响是：

- 20 V 端电流误差由 `net_doping` 的 `+3.785%` 降至
  `cell_reconstructed_total_impurity` 的 `-0.262%`；
- 六个正向锚点的中位绝对电流误差由 `3.911%` 降至 `0.271%`；
- 固定 Sentaurus 状态下，结区电子/空穴边电流比由
  `1.1564/1.0926` 改善为 `0.9914/0.9959`。

因此它是正向 IV 结区电流一致性的关键因素。反向 BV 中影响明显较小：
既有资格验证中总 avalanche source 只下降约 `0.97%--1.29%`，gain-2
电压仅移动约 `0.00048 V`。这不足以解释本次电流代理造成的
`1e6--1e7` 量级误差。

## 后续建议

外部电流替换和电流代理因子分解均已完成。下一项不应直接翻转生产公式，
而应回到 GSS 0.47 方程定义，核对 `aux2` 项究竟是可直接与 QFP drive
相乘的 midpoint density，还是载流子有向通量分解系数；同时核对电子/
空穴势能符号和无向网格边的方向约定。只有完成该公式所有权审查后，
才能授权一个最小、可选启用的 midpoint/actual-SG-current 候选。

本报告比较的是局部 alpha、生成率和二维面积积分，不是 Sentaurus IIC
路径积分；IIC 仍需单独匹配路径定义和积分规则。

## 产物

- 机器可读结论：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/result.json`
- 分阶段汇总：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/stage_summary.csv`
- 支撑级对比：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/support_comparison.csv`
- 外部电流替换：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/external_current_substitution.csv`
- 电流代理因子分解：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/current_proxy_factorization.csv`
- SG 矢量控制：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/sg_vector_current_control.csv`
- 完整报告：`build-release/pn2d-sentaurus-frozen-bv-operator-20260730/report.md`
