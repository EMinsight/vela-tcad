# BVmethods NMOS 原始物理模型等价 BV 执行计划（2026-08-12）

## 范围

本任务只验证原始全耦合物理组合：Fermi、Old Slotboom、Masetti
`DopingDep`、`HighFieldSaturation(GradQuasiFermi)`、`Enormal`、
`SRH(DopingDep)`、`Band2Band(E2)` 和 `Avalanche(Eparallel)`。
主 BV 判据为漏极电流 `1e-4 A/um`。不包含 ABA 路径内部支撑、Transient、
Sentaurus Continuation 专项复现或性能优化。

## 冻结的四组消融

| 组 | SRH 寿命 | Enormal | 目的 |
|---|---|---|---|
| A | 显式常数 `1e-7 s` | 关闭 | 共同基线 |
| B | `SRH(DopingDep)` | 关闭 | SRH 单因素 |
| C | 显式常数 `1e-7 s` | 开启 | Enormal 单因素 |
| D | `SRH(DopingDep)` | 开启 | 原始全模型 |

所有组使用相同 SDE 网格、接触、`Barrier=-0.55 V`、300 K、E2 和
Van Overstraeten 参数。禁止迁移率、雪崩或端电流经验缩放。

Sentaurus 模板：
`reference_tcad/bvmethods_sentaurus2018/source/full_physics_ablation_sdevice.cmd`。
常数寿命控制：
`reference_tcad/bvmethods_sentaurus2018/source/full_physics_constant_srh.par`。

## 当前实现进度

- 已增加 Scharfetter 掺杂相关寿命公式和电子/空穴独立参数。
- 已增加 `total_impurity` 与 `net_doping` 两种掺杂基准；最终基准须由
  Sentaurus A/B 同状态场量决定，不能依靠最终 BV 拟合。
- 已将局部寿命接入 Newton 残差、解析 Jacobian、Gummel 线性化和诊断。
- 常数寿命配置保持兼容；复合单元测试为 25 个算例、72 个断言，全通过。
- Enormal 仍处于趋势原型状态，且部分自洽雪崩源路径目前显式拒绝表面迁移率；
  这是下一阶段的首要实现缺口。

## 验收顺序

1. Sentaurus A/B/C/D 四组生成并导出 0/1/2/4/5/6 V 与目标电流状态。
2. 对 A/B 做寿命和 SRH 的冻结状态逐节点比较。
3. 完成 Enormal 法向场、参数公式及雪崩源/Jacobian 组合支持。
4. 对 A/C 做界面法向场、迁移率和电流增量比较。
5. 对 D 做全模型同状态比较和自洽 BV 扫描。
6. voltage-to-current 为主验收，external-resistor 为独立交叉验证。

最终门槛：BV 误差不超过 2%，高场区场量中位误差不超过 10%，热点位置
不超过一个局部单元，击穿前电流误差不超过 0.3 decade，全局电子/空穴
连续性不守恒不超过 1%。
