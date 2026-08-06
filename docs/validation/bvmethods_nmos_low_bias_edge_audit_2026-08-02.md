# BVmethods NMOS metal-gate 与低偏压逐边审计（2026-08-02）

> **更正（同日后续审计）**：本文最初把 `edge 3457` 的场差解释成求解器
> 平衡态差异。随后通过势梯度/电场交叉校验确认，首因是 Sentaurus TDR
> 区域节点场值被导入器按错误顺序绑定，同时 Vela 掺杂也来自同一错误映射。
> 修正后 1--10 mV 端电流误差降至约 `0.07 dex`。请以后续报告
> `bvmethods_nmos_contact_state_audit_2026-08-02.md` 为准；本文保留用于记录
> 定位过程和为何必须审计导入映射。

## 结论

本轮完成了 metal-gate DD 边界修正、Sentaurus `Barrier=-0.55` 映射核对、
0--10 mV 双求解器重算，以及准费米势/SG 通量逐边定位。

10 mV 时 Vela 漏极电流为 `3.89550e-5 A/um`，Sentaurus 为
`8.564e-10 A/um`，Vela 高 `4.5487e4` 倍（`4.658 dex`）。第一条有物理
传输权重的异常边是硅区 `edge 3457 (node 914 -> 915)`，它从漏极接触节点
直接进入第一层内部节点。

该边不是 SG 算子本身首先失配：把 Sentaurus 的 `psi/phin/phip` 状态送进
完全相同的 Vela SG 算子后，电子粒子线通量为 `-1.7880e9 m^-1 s^-1`；
由 Sentaurus 节点电流密度投影得到的参考值为 `3.1593e9 m^-1 s^-1`，幅值
只差约 `0.247 dex`。Vela 自身状态在同一边给出
`-1.4601e21 m^-1 s^-1`，相对 Sentaurus 状态上的 Vela SG 值高
`11.912 dex`。

因此，五个数量级端电流差的第一处可观测根因位于漏极接触之后第一层
硅节点的静电势/载流子状态，而不是 metal gate 上错误钉扎准费米势，也不
是第一条异常边上的 SG 通量公式。后续应优先审计 ohmic 接触内建势、掺杂/
BGN 参考能级以及接触邻接节点的 Poisson 平衡态。

## 1. Metal-gate DD 边界语义

修正后的语义为：

- metal gate 只对 Poisson 未知量施加 Dirichlet 边界；
- 不在绝缘层 gate 节点钉扎 `phin/phip`，也不恢复或校验这些节点的载流子行；
- Gummel、Newton、DCSweep、恢复路径和延拓路径均携带同一 contact spec；
- DD 解校验和准费米势边界检查排除 metal gate。

回归测试覆盖了“氧化层 metal gate 仅约束 `psi`，准费米势保持自由”的行为。

## 2. Barrier 到 flatband 的映射

官方 deck 使用：

```text
{ Name="gate" Voltage=0.0 Barrier=-0.55 }
```

Sentaurus 日志的 gate potential 为 `+0.55 V`。因此本算例的映射为：

```text
psi_gate = Voltage - Barrier
         = 0 - (-0.55)
         = +0.55 V
```

Vela 的 metal-gate 约定是 `psi_gate = bias - flatband_voltage`，故应配置
`flatband_voltage = -0.55 V`。边界测试已固定该符号约定。

## 3. 低偏压端电流

| Vd (V) | Vela (A/um) | Sentaurus (A/um) | `|Vela/S|` | 差值 (dex) |
|---:|---:|---:|---:|---:|
| 0.001 | 3.92949e-6 | 1.029e-10 | 3.8187e4 | 4.582 |
| 0.002 | 7.85212e-6 | 2.015e-10 | 3.8968e4 | 4.591 |
| 0.005 | 1.95753e-5 | 4.737e-10 | 4.1324e4 | 4.616 |
| 0.010 | 3.89550e-5 | 8.564e-10 | 4.5487e4 | 4.658 |

两者均近似线性响应，但斜率已在 1 mV 处相差约 `3.8e4`，所以该误差不是
雪崩开启后才产生，而是低场平衡态/电导问题。

## 4. 第一条异常边

10 mV、`edge 3457 (914 -> 915)`：

| 量 | Vela 状态 | Sentaurus 状态 |
|---|---:|---:|
| `psi(node 914)` (V) | 0.476188 | 0.438345 |
| `psi(node 915)` (V) | 0.476762 | -0.222467 |
| `Delta phin` (V) | -7.9621e-6 | -2.5469e-7 |
| Vela SG 粒子线通量 (`m^-1 s^-1`) | -1.4601e21 | -1.7880e9 |
| Sentaurus 原生投影参考 (`m^-1 s^-1`) | -- | 3.1593e9 |

`Delta phin` 只相差 `1.495 dex`，但 SG 通量相差 `11.912 dex`。决定性差异是
内部节点 915 的 `psi`：Vela 比 Sentaurus 高约 `0.699 V`。在 300 K 下，
该势差通过 Boltzmann 因子产生约 12 个数量级的电子浓度/导电能力差，和逐边
通量差一致。

网格中还存在从漏极接触出发、`couple_m=0` 的几何边；逐边脚本将它们单独
标识并排除在有效 SG 边排序之外，避免把无对偶传输权重的边误判为算子异常。

## 5. 可复现产物

- `scripts/run_bvmethods_nmos_low_bias_sentaurus_vm.py`：复制官方 IIC deck，
  独立运行 0/1/2/5/10 mV Sentaurus 工况并回收结果。
- `scripts/analyze_bvmethods_nmos_low_bias_edges.py`：合并区域场、调用 Vela SG
  边探针、限定硅区有效边、计算漏极 hop 距离并输出第一异常边。
- `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/low_bias_edge_compare_20260802/first_abnormal_edge_summary.csv`
  汇总五个偏压点的定位结果。
- 同目录的 `edge_compare_*.csv` 保存每条边的端点势、准费米势差、物理粒子
  线通量和对数差。

## 6. 下一验证顺序

1. 固定 `edge 3457` 和相邻三角形，逐项核对 Vela/Sentaurus 的净掺杂、
   OldSlotboom BGN、`ni_eff`、接触平衡载流子和 ohmic `psi` 公式。
2. 先要求 0 V 下 node 914/915 的 `psi/n/p` 与 Sentaurus 对齐，再重复
   1--10 mV 扫描；不要先调迁移率拟合端电流。
3. 平衡态对齐后，用 Sentaurus 状态上的 Vela SG 结果作为算子门槛：首层
   有效边误差目标 `<0.3 dex`，端电流目标先收敛到 `<1 dex`。
4. 低场 IV 对齐后再恢复 1/2/4/5/6 V 的电场、离化系数、离化积分和
   雪崩电流比较，最后执行 0--7 V BV 扫描。
