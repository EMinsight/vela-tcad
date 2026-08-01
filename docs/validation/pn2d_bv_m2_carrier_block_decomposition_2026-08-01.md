# PN2D BV M2 carrier-block 线性解分解

日期：2026-08-01

## 技术摘要

本轮只读诊断通过，typed outcome 为：

`carrier_block_linear_solve_decomposed`

在共享 M2 网格、逐节点掺杂和 SG/Laux 配置不变的条件下，对 `-18`、
`-19.5`、`-19.7`、`-20 V` 的 Vela baseline 与 joint Sentaurus-QFP
冻结状态进行了两套独立运行。完整 carrier 线性闭合误差最大为
`7.195e-16`，两轮五类原始输出逐字节一致。未修改物理参数、生产默认值、
continuation 或验收门限。

核心结论是：joint-QFP 冻结状态把 carrier Newton 步集中到结区两侧的两个
超软奇异方向。雪崩 Jacobian 提供了几乎全部电子—空穴交叉耦合，并显著
放大这些方向上的步长；没有发现 row scaling、复合项、接触行或线性求解闭合
错误，也没有发现更新方向翻转。

## 冻结合同

- 网格和掺杂：共享 M2 `mesh.json` 与 `doping.csv`，结区为 `x=1.0 um`。
- 偏压：`-18`、`-19.5`、`-19.7`、`-20 V`。
- 状态：Vela baseline 与只替换电子/空穴 QFP 的 joint Sentaurus-QFP。
- SVD/条件数区域：仅保留非接触约束的 216 个 carrier QFP 行和列。
- 接触 identity 行从 SVD 和条件数中排除。
- 反事实矩阵：`full`、`row_scaled`、`no_cross_carrier`、
  `no_recombination`、`no_avalanche`、`transport_only`。
- 所有反事实矩阵使用同一个冻结 full residual，不重算状态。
- 独立重复：两次完整运行；每个 case 的 summary、columns、singular modes、
  solve variants、solve nodes 五类文件必须联合哈希一致。

## 数值正确性

| 检查 | 最大值 | 结论 |
|---|---:|---|
| full carrier 线性相对闭合误差 | `7.195e-16` | 通过 |
| row-scaled 与 full 步相对差 | `2.098e-13` | 通过 |
| SVD RHS/step 能量闭合偏差 | `4.885e-15` | 通过 |
| 纯 transport 交叉块 / full Frobenius 范数 | `0` | 通过 |
| 两套独立运行 | 16 个 case 联合哈希一致 | 通过 |
| C++ 诊断测试 | 14 assertions | 通过 |
| 完整 `test_newton_solver` | 1037 assertions / 70 cases | 通过 |
| runner 集成测试 | 1 case | 通过 |

左乘 continuation row weights 不改变精确线性解。raw 与 production
row-scaled 矩阵的尺度跨度使相对 SVD 阈值只能解析部分秩，因此不能把它们的
`resolved_condition_number` 当成完整条件数。满秩的 L2 行列平衡矩阵用于跨状态
比较。

## 条件数、残差和步长

| Bias | State | L2 平衡条件数 | electron residual | hole residual | carrier 步范数 (V) |
|---:|---|---:|---:|---:|---:|
| -19.5 | baseline | `65.80` | `1.299e-10` | `1.591e-9` | `1.989e-3` |
| -19.5 | joint QFP | `381.53` | `4.787e-8` | `5.553e-8` | `1.429` |
| -19.7 | baseline | `77.16` | `1.404e-10` | `1.823e-9` | `1.997e-3` |
| -19.7 | joint QFP | `426.91` | `5.836e-8` | `6.672e-8` | `1.385` |
| -20.0 | baseline | `105.06` | `1.758e-10` | `1.799e-9` | `2.009e-3` |
| -20.0 | joint QFP | `543.26` | `8.367e-8` | `9.435e-8` | `1.296` |

在 knee 区间，joint-QFP 的 L2 平衡条件数是 baseline 的约 `5.17x` 至
`5.80x`。electron residual 放大约 `368x` 至 `476x`，hole residual 放大
约 `35x` 至 `52x`，carrier 步范数放大约 `645x` 至 `718x`。因此大步长既
包含冻结 Sentaurus QFP 对 Vela 方程的不平衡，也包含该不平衡沿软模态的放大；
不能只用 raw condition number 解释。

## 主导奇异方向

实际 production sparse-solver 步直接投影到 free carrier-block 的右奇异向量。
在 knee 区间，joint-QFP 前两个模态承载：

| Bias | 第一模态能量 | 第二模态能量 | 合计 | 相对奇异值范围 |
|---:|---:|---:|---:|---:|
| -19.5 | `60.24%` | `36.10%` | `96.34%` | `6.86e-16`–`8.92e-16` |
| -19.7 | `56.57%` | `39.06%` | `95.63%` | `7.63e-16`–`9.12e-16` |
| -20.0 | `58.70%` | `36.18%` | `94.88%` | `8.72e-16`–`1.14e-15` |

主导模态分别以 hole 与 electron 分量为主，形成成对的结区两侧软方向，而不是
单个错误行。最大列范数稳定出现在远 n 区 `x=1.5 um`；最大实际更新和超软
模态却在结区附近，说明“最大局部 Jacobian 列”不是大步长根因，关键是小奇异
方向的全局组合。

## 电子—空穴耦合来源

| Bias | State | 去 cross-carrier 步差 | 方向余弦 | 去 avalanche 步差 | 方向余弦 | 去 recombination 步差 |
|---:|---|---:|---:|---:|---:|---:|
| -19.5 | baseline | `2.19%` | `0.99976` | `2.22%` | `0.99975` | `3.20e-6` |
| -19.5 | joint QFP | `38.01%` | `0.99493` | `25.11%` | `0.99869` | `9.99e-7` |
| -19.7 | baseline | `2.87%` | `0.99959` | `2.93%` | `0.99957` | `3.19e-6` |
| -19.7 | joint QFP | `40.26%` | `0.99427` | `26.64%` | `0.99876` | `1.17e-6` |
| -20.0 | baseline | `4.96%` | `0.99877` | `5.06%` | `0.99872` | `3.17e-6` |
| -20.0 | joint QFP | `46.42%` | `0.99293` | `31.53%` | `0.99857` | `1.52e-6` |

纯 transport 的电子—空穴交叉块严格为零。`full - no_avalanche` 贡献了约
100% 的交叉块范数，SRH 交叉项相对可忽略。交叉块的 Frobenius 范数只占 full
矩阵约 `1e-13`，但它作用于相对奇异值约 `1e-15` 的软方向，因此能改变
joint-QFP 步长 `38%`–`46%`。方向余弦仍高于 `0.992`，现有证据支持“幅值
放大”，不支持符号翻转。

`no_cross_carrier` 与 `no_avalanche` 不是可相加的分解：前者只删除交叉块并
保留 avalanche 对角块，后者同时删除 avalanche 的对角和交叉贡献。

## 与结区节点净掺杂的关系

本轮结果支持“相关，但不是直接局部根因”的判断：

- joint-QFP 的 4 个偏压、每个偏压 top-10 实际更新节点共 `40/40` 个均位于
  `x=1.0 um` 结区的 `0.25 um` 范围内；baseline 对应为 `0/40`。
- knee 区最强两个 joint-QFP 模态通常位于 `x=0.75 um` 的 p 型肩部和
  `x=1.25 um` 的 n 型肩部，节点净掺杂仍分别是 `-1e17` 与 `+1e17 cm^-3`。
- 这些节点距补偿结区列通常为 4 个三角网格图距离，且不直接属于包含
  `x=1.0 um` 补偿节点的三角形。
- 少量次级模态会落到直接接触补偿结区三角形的 `x=0.9375/1.0625 um`
  节点，但不是 95% 步能量的主要来源。

配置和实现复核确认：当前 BV 使用节点净掺杂。Poisson 固定电荷直接使用节点
`ND-NA`；SG 边迁移率使用两个端点节点净掺杂的算术平均。只有显式选择
`cell_reconstructed_total_impurity` 才会启用三角形三节点总杂质平均，而当前
配置没有选择该分支。因此，不存在可归因的“补偿结区三角形掺杂平滑公式”。
更准确的空间定位是结区肩部/耗尽区全局模态，以及节点净掺杂、边端点平均和
控制体积之间的离散关系。

## 限制

- joint-QFP 是冻结替换状态，不是 Vela 自洽分支；其较大 residual 是判别实验
  的预期组成，不能直接视为生产求解失败。
- 本实验只分解 Vela 方程和 Jacobian，尚未取得 Sentaurus 的对应离散
  carrier-block Jacobian，因此不能从条件数本身判定哪个仿真器的离散更正确。
- Frobenius 范数会掩盖近零奇异方向上的小耦合放大，必须与模态投影共同解读。

## 下一步

保持 SG/Laux 与生产默认值不变，执行结区肩部 `x=0.75–1.25 um` 的只读审计：

1. 将两个超软模态分别投影到 transport、avalanche diagonal 和 avalanche
   cross-carrier Jacobian 分项，确定哪个分项改变对应奇异值和 RHS 投影。
2. 对共享 M2 网格逐节点/逐边核对净掺杂、端点平均、控制体积剂量和
   `p -> compensated -> n` 结区所有权；不修改输入，只比较 Vela 与 Sentaurus
   导出语义。
3. 若需做掺杂反事实，只允许 frozen-state、opt-in 的诊断变体，并保持离散总剂量
   不变；不得直接改生产节点掺杂、控制体积策略或 SG/Laux 默认值。

## 证据

- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/result.json`
- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/case_summary.csv`
- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/dominant_columns.csv`
- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/dominant_singular_modes.csv`
- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/determinism.csv`
- `build-release/pn2d-bv-m2-carrier-block-decomposition-20260801/report.html`
