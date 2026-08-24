# TransportModels DD/DG 固定契约统一验收（2026-08-24）

## 结论

本次固定了 Sentaurus 2022 材料参数、Fermi/OldSlotboom、`sentaurus_default` SRH 密度耦合、偏压点阵及最终 DG 量子势契约。DD 使用同一材料与偏压契约重新计算 21 点 Id–Vg 和 21 点 Id–Vd；DG 使用已完成的 42 点证据并通过契约逐字段审计。

统一验收：**未通过**。

## 主曲线指标

| 分支 | Id–Vg 过渡区最大对数误差 | Id–Vg 导通区最大相对误差 | Id–Vd 最大相对误差 | Id–Vd 2 V 误差 | 结果 |
|---|---:|---:|---:|---:|---|
| DD | 0.031852 dex | 2.930% | 1.361% | 1.343% | 主曲线通过；深关断未解析 |
| DG | 0.032874 dex | 2.710% | 2.432% | 0.976% | 通过 |

## 深关断前三点

| 分支 | Vg (V) | 相对误差 | 对数误差 | Id/abs(KCL) | 状态 |
|---|---:|---:|---:|---:|---|
| DD | -1.00 | 9.016% | 0.041036 dex | 8.460 | numerically_unresolved |
| DD | -0.84 | 40.317% | 0.224151 dex | 1.428 | numerically_unresolved |
| DD | -0.68 | 14.440% | 0.067731 dex | — | numerically_unresolved |
| DG | -1.00 | — | 0.004265 dex | 378.273 | pass |
| DG | -0.84 | — | 0.004763 dex | 377.241 | pass |
| DG | -0.68 | — | 0.011539 dex | 537.642 | pass |

深关断采用独立门槛：对数误差不超过 0.15 dex，且 `Id/abs(KCL) >= 10`；不与导通区相对误差混用。

## 固定基线

- 契约：`D:\code-repo\vela-tcad\configs\regression\transportmodels_dd_dg_sentaurus2022_v1.json`
- 契约 SHA-256：`0513b486661fb7ae788f108782d79a5ded75e3f3c7f55216c51eaf0eec14cbec`
- 材料 SHA-256：`e00fa2d0585aa75d624d8275b91e76b7e412ca0d60c33bf591ba9ffa44bac12f`
- DD/DG 受控差分：通过，仅 electron_quantum_potential 不同
- DG 配置逐字段审计：通过
- DD 运行器 SHA-256：`2cc86fa83fc0ec70e74b45314ca7e8bd201972a38adda37255e452e2fdba745e`
- DG 证据运行器 SHA-256：`2cc86fa83fc0ec70e74b45314ca7e8bd201972a38adda37255e452e2fdba745e`

## 产物

- JSON 报告：`D:\code-repo\vela-tcad\docs\validation\transportmodels_dd_dg_fixed_contract_v1_2026-08-24.json`
- 对比图：`D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_dg_fixed_contract_v1_2026-08-24\transportmodels_dd_dg_fixed_contract_comparison.png`
- DD 工作流清单：`D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dd_dg_fixed_contract_v1_2026-08-24\runs\dd\workflow_manifest.json`
