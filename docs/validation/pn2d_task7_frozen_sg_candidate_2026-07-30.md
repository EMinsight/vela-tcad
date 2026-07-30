# PN2D Task 7 冻结状态完整 SG 电流候选预审

日期：2026-07-30

## 目标和执行合同

在完全相同的 Sentaurus `-19.7/-19.8 V` 冻结状态和 Van Overstraeten
模型下，对比三条单轴支路：

1. 当前生产 triangle-GSS scalar proxy；
2. sign-correct midpoint-only 负对照；
3. opt-in `element_edge_sg_gss_laux` 完整 SG 电流向量及
   `element_vertex_box_measure` 匹配几何。

所有支路均为 `postprocess_only`。状态没有推进，碰撞电离源项没有进入
连续性方程、Jacobian 或 continuation，生产默认值未修改。

两次独立重放分别写入：

- `build-release/pn2d-task7-frozen-sg-candidate-a-20260730`
- `build-release/pn2d-task7-frozen-sg-candidate-b-20260730`

## 支撑对齐

Sentaurus 本次导出的是节点 `Jn/Jp/ImpactIonization`，不是原生单元量。
因此不能把同一个节点值复制到每个相邻单元，并把它当作 element-local
reference。

本预审执行以下守恒投影：

- 将 Vela element-vertex 源项对相邻单元求和，得到物理节点源项；
- 将 Vela GSS/Laux 单元电流向量按 element-vertex measure 加权回物理
  节点；
- 电子粒子通量乘以 `-q`、空穴粒子通量乘以 `+q` 后，再与 Sentaurus
  conventional-current vector 比较；
- active source 定义为参考节点源项不低于本偏压峰值的 `1e-6`；
- `1e-5/1e-7/1e-8` 三档敏感性检查得到完全相同的 9 个 active 节点和
  相同判定。

未经节点投影的 element-to-node 行保留为 unsupported diagnostic，不参与
候选门限。

## 三支路积分结果

| 支路 | -19.7 V | -19.8 V |
|---|---:|---:|
| 当前 triangle baseline / Sentaurus | 2.4855322e7 | 7.9270586e6 |
| sign-correct midpoint-only / Sentaurus | 0.486029 | 0.487386 |
| 完整 SG/Laux vector / Sentaurus | **1.009537** | **1.009483** |

sign-correct midpoint-only 仍缺少约一半源项，继续作为负对照。完整 SG
向量候选的积分误差为 `0.00410-0.00412 dex`。

## 完整 SG 向量候选的 matching-support 指标

| 指标 | -19.7 V | -19.8 V | 门限 |
|---|---:|---:|---:|
| electron current median / P95 (dex) | 0.00366 / 0.00487 | 0.00367 / 0.00484 | 0.05 / 0.15 |
| hole current median / P95 (dex) | 0.00340 / 0.00436 | 0.00339 / 0.00433 | 0.05 / 0.15 |
| active node source median / maximum (dex) | 0.00404 / 0.00518 | 0.00395 / 0.00514 | 0.10 / 0.30 |
| electron vector median / P95 angle | 0.0117° / 0.0241° | 0.0106° / 0.0228° | diagnostic |
| hole vector median / P95 angle | 0.00247° / 0.00713° | 0.00241° / 0.00680° | diagnostic |
| nonzero vector direction agreement | 100% | 100% | 100% |

全部固定状态门限通过。

## 确定性和闭合

- 两次运行的 27 个数值 CSV/过程产物逐字节一致；
- 两组 imported state、baseline process、SG-vector process、因子分解和
  sign-correct 负对照均一致；
- 两次运行均为 observation-only，solver-coupled 和 residual feedback
  记录为零；
- 完整 SG 向量源项总和分别只比 Sentaurus 高 `0.9537%` 和 `0.9483%`。

## 判定

本阶段 typed outcome：

`complete_sg_vector_fixed_state_prequalified`

含义：

1. 完整 `element_edge_sg_gss_laux` 电流向量和匹配几何通过 Task 7
   固定状态预审；
2. sign-correct midpoint-only 被正式封存为负对照；
3. 此结果只授权下一阶段 Task 7 的单轴、自洽、双运行 exact-lattice
   候选验证；
4. 尚未满足 knee-window、`V_slope/V_break`、非单调区间、全局电流和
   Task 6 自洽反馈门限，因此不授权 Task 8，也不授权生产默认值修改。

## 产物

- `build-release/pn2d-task7-frozen-sg-candidate-score-20260730/result.json`
- `build-release/pn2d-task7-frozen-sg-candidate-score-20260730/task7_frozen_candidate_scorecard.csv`
- `build-release/pn2d-task7-frozen-sg-candidate-score-20260730/determinism.csv`
- `scripts/score_pn2d_task7_frozen_sg_candidate.py`
