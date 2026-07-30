# PN2D GSS `aux2` 方程所有权与符号审查

日期：2026-07-30

## 目的与边界

本审查回答两个问题：

1. GSS `aux2` 表达式究竟定义物理 midpoint density，还是有向通量分解系数；
2. Vela `triangle_gss_gradqf_truncated` 是否正确继承了电子、空穴和边方向约定。

审查只读取归档方程、开源参考代码和 `-19.7/-19.8 V` 冻结状态。没有推进
求解状态，没有向连续性方程、Jacobian 或 continuation 写回源项，也没有修改
生产默认值。

## 参考定义

GSS 用户指南 9.100、9.103、9.107、9.108 在等温条件下给出

```text
alpha = (psi_i - psi_j) / (2 Vt)
n_mid = n_i aux2(alpha)  + n_j aux2(-alpha)
p_mid = p_i aux2(-alpha) + p_j aux2(alpha)
aux2(x) = 1 / (1 + exp(x))
```

本机 `Genius-TCAD-Open` 修订
`543da8452d5dfd33e6f8c457f962f6f670f0fce7` 的
`include/math/jflux1q.h` 标注来源为 `GSS 0.4x`，其中 `nmid/pmid` 与上述
方程逐项一致。文件 SHA-256 为
`fce16bffb62644054dd4e3934dacb043f21861cfcebaf9ef3c418ee831f9e625`。

`aux2` 在这里确实构成一个 edge midpoint carrier density，但它属于完整
SG 漂移—扩散分解。以等温电子式为例：

```text
J_n = q mu_n [
    n_mid * (psi_i - psi_j) / h
  + Vt * aux1(alpha) * (n_j - n_i) / h
]
```

因此它不是可以脱离 `aux1` 扩散项和有向 SG 边通量、单独与
`abs(delta QFP)/h` 相乘的碰撞电离电流定义。

独立代码证据同样明确：Genius DDM1 的碰撞电离生成率使用
`GIIn = alpha_n * abs(Jn) / q` 和 `GIIp = alpha_p * abs(Jp) / q`，其中
`Jn/Jp` 是前面已经计算的完整有向 SG 电流，而不是 midpoint-only proxy。

## Vela 当前符号映射

令 `alpha=(psi_i-psi_j)/(2Vt)`。当前 triangle-GSS 路径实际计算：

```text
n_mid_current = n_i aux2(-alpha) + n_j aux2(alpha)
p_mid_current = p_i aux2(alpha)  + p_j aux2(-alpha)
```

这两个式子都与 GSS 参考式相反。旧实现测试能证明代码与当时抄录的表达式
一致，但测试名称中的 “published carrier orientation” 并不能证明抄录时
已经完成静电势、有效导带/价带势和载流子符号的映射。

同时需要注意：参考式和当前式在同时交换端点、密度和电势后都保持不变。
因此单纯的 endpoint-swap 测试不能发现本次符号错误；必须使用方程定义或
有效带边到静电势的映射作为外部基准。

## 冻结状态数值审查

可重复脚本：
`scripts/audit_pn2d_gss_aux2_ownership.py`。

| 指标 | -19.7 V | -19.8 V |
|---|---:|---:|
| 当前 triangle 公式重建闭合 L2 | 1.11e-16 | 1.43e-16 |
| GSS 参考 midpoint 对 Vela SG-edge midpoint L2 | 6.64e-17 | 9.40e-17 |
| GSS 参考式端点交换最大相对误差 | 0 | 0 |
| 热点当前 midpoint / GSS 参考 midpoint | 7.2531e7 | 2.3350e7 |
| GSS 参考 midpoint-only 源项 / Sentaurus | 0.48603 | 0.48739 |
| 实际 SG/Laux 向量源项 / Sentaurus | 1.00954 | 1.00948 |

两点的最大热点均为 electron、cell 17、edge 35、node 11→14。以
`-19.7 V` 为例：

- 当前 triangle midpoint：`1.1447218451e18 m^-3`；
- GSS 参考 midpoint：`1.5782596317e10 m^-3`；
- Vela SG-edge midpoint：`1.5782596317e10 m^-3`；
- 当前 proxy：`1.9276121943e4 A/m^2`；
- GSS 参考 midpoint proxy：`2.6576521842e-4 A/m^2`；
- 完整 Vela SG 电流：`2.7842713665e-4 A/m^2`；
- Sentaurus 电流：`2.4904353319e-4 A/m^2`。

符号修正会消除该热点的七个数量级放大，但把参考 midpoint 继续放入
triangle scalar proxy 后，总源项只有 Sentaurus 的约 `48.6%-48.7%`。
这说明“只翻转 midpoint 符号”不是完整修复：它仍丢失完整 SG 通量、
二维电流向量重建和相匹配的源项支撑。相反，已经存在的
element-edge SG/GSS-Laux 向量对照在同一冻结状态下只高约 `0.95%`。

## 判定

审查结论为：

`gss_aux2_sign_transcription_and_operator_ownership_mismatch_confirmed`

具体含义：

1. `aux2` 组合可解释为 GSS SG 分解中的 midpoint density；
2. 当前 triangle 路径对电子和空穴的静电势符号均与参考式相反；
3. GSS 碰撞电离消费完整 SG 电流，不支持把 `mu*n_mid*abs(grad QFP)`
   当作独立、完整的离散电流算子；
4. endpoint orientation 不是根因，因为两个表达式自身都满足端点交换
   不变性；根因是势变量/载流子符号映射和算子所有权；
5. sign-only 修正候选被否决，不能进入完整 BV 曲线验证；
6. 若后续授权最小 opt-in 候选，应优先复用真实
   `element_edge_sg_gss_laux` 电流向量及其匹配源项几何，同时保留
   sign-correct midpoint-only 作为负对照。生产默认值保持不变。

## 产物

- `build-release/pn2d-gss-aux2-ownership-audit-20260730/result.json`
- `build-release/pn2d-gss-aux2-ownership-audit-20260730/gss_aux2_ownership_details.csv`
- `build-release/pn2d-gss-aux2-ownership-audit-20260730/gss_aux2_ownership_summary.csv`
