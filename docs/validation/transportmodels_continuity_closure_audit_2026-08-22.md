# TransportModels 接触通量—SRH 连续性闭合审计（2026-08-22）

## 审计范围

审计对象为 TransportModels Id-Vg 深关断区，重点复核：

1. 接触节点连续性项与 Dirichlet 行替换；
2. 接触边 Scharfetter-Gummel（SG）电子/空穴通量；
3. 硅区 SRH 控制体积分及其符号；
4. 局部载流子行收敛与全局端口闭合的关系。

主审计点为 DD -0.68 V，并以 DG -0.52 V 和已通过的 DD -1 V 交叉验证。

## 离散闭合关系

Newton 全局闭合诊断使用未进行接触边界行替换的物理连续性项。对自由节点求和后，内部边通量成对消去，得到：

```text
接触节点 SG 通量总和 = 自由节点复合/产生源积分总和
```

实际 Newton 残差中的接触准费米行仍是 Dirichlet identity。因此局部载流子行门槛只检查自由节点，必须另设全局闭合门槛检查被替换接触行背后的物理通量。

现有合成单元测试对该符号和接触/自由节点划分已有覆盖。本次真实 MOS 数据没有发现闭合公式本身的符号错误。

## 求解器内部闭合结果

下表数值为统一缩放后的连续性量；闭合比不受缩放影响。

| 模型 | Vg/V | Newton 变体 | 局部违规 | 电子闭合比 | 空穴闭合比 | 电子接触通量 | 电子积分源 |
|---|---:|---|---:|---:|---:|---:|---:|
| DD | -1.00 | 全步 | 0 | 4.196e-2 | 7.001e-6 | -1.164e-13 | -1.115e-13 |
| DD | -0.68 | 全步 | 87 | 9.737e-1 | 7.157e-11 | -4.336e-12 | -1.142e-13 |
| DD | -0.68 | 初始阻尼 0.5 | 0 | 1.029e0 | 3.317e-15 | +3.949e-12 | -1.142e-13 |
| DG | -0.52 | 全步 | 0 | 9.358e-1 | 2.886e-15 | -1.772e-12 | -1.137e-13 |

关键证据：

- 每个审计点的电子与空穴积分源完全相同，符合 SRH 成对产生/复合。
- DD -0.68 V 和 DG -0.52 V 的空穴闭合达到约 1e-11 至 1e-15，证明 SRH 公式、符号和控制体权重一致。
- DD -0.68 V 的 0.5 阻尼状态已经没有局部载流子行违规，但电子闭合仍失败。
- DD -0.68 V 的电子接触通量在全步与 0.5 阻尼之间由 -4.336e-12 变为 +3.949e-12，积分源仅由 -1.141604219e-13 变为 -1.141604141e-13。电子通量的符号/数量级敏感性是数值分辨率特征，不是物理源项变化。

## 物理端口与 SRH 积分

| 模型 | Vg/V | 电子 reference | ΣIe/(A/um) | ΣIh/(A/um) | SRH 产生/(A/um) | 电子误差 | 空穴误差 | |KCL|/SRH |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| DD | -0.68 | drain, 1.1 V | +2.033e-16 | +5.349e-18 | 5.349e-18 | 3.702e1 | 9.806e-5 | 3.702e1 |
| DG | -0.52 | drain, 1.1 V | -3.906e-16 | +5.330e-18 | 5.328e-18 | 7.231e1 | 3.350e-4 | 7.431e1 |

空穴端口电流能够复现 SRH 产生电流，电子端口和四端 KCL 未达到深关断绝对精度。

## 接触边 SG 通量

| 模型 | Vg/V | 接触 | 活跃电子边 | 零电子电流边 | 物理输出 phin 差为零的边 | Ie/(A/um) |
|---|---:|---|---:|---:|---:|---:|
| DD | -0.68 | source | 30 | 30 | 30 | 0 |
| DD | -0.68 | drain | 30 | 0 | 30 | +2.033e-16 |
| DG | -0.52 | source | 30 | 1 | 1 | -3.956e-14 |
| DG | -0.52 | drain | 30 | 0 | 2 | +3.917e-14 |

接触边求和与 `terminal_balance.csv` 的电子端口电流在约 1e-16 相对精度内一致，因此端口后处理没有漏边或重复计数。

DD -0.68 V 的源漏平均 n 型掺杂均为 7.985e20 cm^-3。`configureQuasiFermiReferences` 只在掺杂严格更大时替换参考接触；由于网格接触顺序中 drain 先于 source，电子全局 reference 被固定为 drain=1.1 V。该 reference 只能让漏端附近的未知量靠近零，另一端仍位于约 -1.1 V 的增量基值。

## 根因判断

| 候选原因 | 结论 | 证据 |
|---|---|---|
| SRH 广义 Fermi 公式 | 排除为主因 | 电子/空穴积分源完全相同，空穴闭合到 1e-15 |
| SRH 符号或控制体积 | 排除为主因 | 空穴端口与 SRH 产生误差仅 1e-4 至 3e-4 |
| 接触边枚举/方向 | 排除为主因 | 边求和与端口后处理一致 |
| 接触 Dirichlet 行替换 | 设计正确，但必须配合全局门槛 | 局部门槛通过时电子全局闭合仍可失败 |
| 解析 Jacobian 缺项 | 已由前序 JVP 审计排除 | 1e-8 V 热点 JVP 误差约 1e-6 |
| 电子准费米/SG 极低通量分辨率 | 当前主因 | 电子闭合单独失败，且接触通量随阻尼变号 |
| 单一全局准费米 reference | 高优先级诱因 | 对称源漏只能有一端以零附近增量表示 |

## 审计中发现的可复现性缺口

当前 restart/failed-state CSV 保存物理 `phin`、`phip`，但不保存：

- `phinIncrement`、`phipIncrement`；
- `electronQfReference_V`、`holeQfReference_V`。

因此重新读取状态时，低于绝对电势 ULP 的内部准费米增量会丢失。现有同一次运行生成的 `contact_edges.csv` 和全局闭合字段可信，但失败状态不能被逐边精确复放。

## 建议的实施顺序

1. 扩展 restart/failed-state 格式，持久化准费米 increment 和 reference，并保持旧 CSV 向后兼容。
2. 用同一失败状态进行 double、long-double/补偿求和的接触通量 A/B，区分“状态量化”与“求和舍入”。
3. 评估按接触盆地或连通分区设置电子准费米 reference，避免源漏两端只能保护一端。
4. 若分区 reference 仍不足，再考虑局部高精度电子 SG 通量或增量式线性求解变量。
5. 保留现有全局连续性硬门槛，不通过提高 Newton floor 接受未解析点。

## 产物

- 可重复脚本：`scripts/audit_transportmodels_continuity_closure.py`
- 机器可读汇总：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/continuity_closure_audit/audit.json`
- 求解器闭合表：同目录 `solver_closure.csv`
- 物理端口表：同目录 `physical_balance.csv`
- 接触边表：同目录 `contact_edge_balance.csv`

## 准费米重启与精度 A/B（追加审计）

已扩展 `DDSolution` 重启格式，在原有物理 `phin/phip` 之外保存：

- `electron_qf_increment_V`、`hole_qf_increment_V`；
- 逐节点 `electron_qf_reference_V`、`hole_qf_reference_V`。

读取器继续兼容旧格式，并校验 `physical_qf = reference + increment`。Newton
热启动在 reference 与当前求解器一致时直接恢复自由节点 increment，不再先由物理
准费米势反推，从而避免在 1.1 V 基值附近丢失子 ULP 增量。

在 DD `Vg=-0.68 V` 状态上执行了以下冻结状态 A/B：

1. 单一全局电子 reference；
2. source/drain/n+ poly gate 最近接触盆地 reference；
3. 普通 double 顺序求和；
4. Neumaier 补偿求和；
5. long-double SG 边通量与 long-double 接触求和。

| 接触 | 全局 double 电子电流/(A/um) | 补偿-double 差值 | long-double 差值 | 分区电子电流/(A/um) |
|---|---:|---:|---:|---:|
| drain | 2.033391076734552e-16 | -2.465e-32 | -1.726e-31 | 2.033391076734552e-16 |
| source | 0 | 0 | 0 | 0 |
| gate | 0 | 0 | 0 | 0 |
| substrate | -8.378052409861056e-21 | 0 | -7.523e-36 | -6.140117563354096e-21 |

分区中 source、drain、gate 分别包含 1600、1661、54 个最近接触节点。重构后的
物理准费米势与旧物理输出列最大差 `2.206e-16 V`，约为 1.1 V 附近一个 ULP；
取消抵消所需的 increment/reference 分量则被独立保存。

### 结论

- drain 的补偿求和相对改变量约 `1.21e-16`，long-double 相对改变量约
  `8.49e-16`，可排除普通 double 多边求和是 `1e-16 A/um` 端口误差主因。
- source 在全局/分区 reference 与 double/补偿/long-double 下均为严格零，说明
  该保存状态的 source 接触边电子准费米 increment 差已经为零；后处理提高精度
  无法恢复不存在的状态差。
- drain 和 gate 在分区后完全不变，验证同一 reference 盆地内的局部 SG 通量保持；
  substrate 的微小电子分量变化表明最近几何接触分区会把抵消转移到盆地交界边，
  尚不能直接作为生产 Newton 坐标模式。
- 根因进一步收缩为：Newton 更新/线性解在源端局部增量形成之前已经受到状态尺度
  或收敛门槛限制，而不是重启量化、接触边累加或 double SG `expm1` 本身。

追加产物：

- 执行脚本：`scripts/run_transportmodels_qf_reference_precision_ab.py`
- A/B 汇总：`build-release/reference_tcad/transportmodels_sentaurus2022/reports/idvg_deep_off_precision_20260822/qf_reference_precision_ab_20260822/summary.csv`
- 完整执行记录：同目录 `execution.json`
- 无损全局状态：同目录 `global_solve/final_state.csv`
- 分区坐标状态：同目录 `partition_state.csv`
