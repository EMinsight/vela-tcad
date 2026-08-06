# BVmethods NMOS Fermi-Dirac、Ohmic 与广义 SG 实现验证

日期：2026-08-02

## 1. 结论

Vela 的 coupled Newton 漂移扩散路径已增加一致的 Fermi-Dirac 统计：

1. 节点电子/空穴密度使用归一化完全 Fermi-Dirac 积分 `F_{1/2}`；
2. ideal Ohmic 接触用同一密度模型数值求解 `n-p=Nnet`，不再用经验电压修正；
3. 连续性方程使用广义 Einstein/Scharfetter-Gummel 边通量；
4. SRH/Auger 的平衡乘积、初值/恢复、端电流和连续性诊断使用同一统计配置；
5. `boltzmann` 保持默认兼容行为，纯 density-form Gummel 对
   `fermi_dirac` 明确拒绝，避免静默混用经典 Einstein 关系。

BVmethods 实网格的五阶段 0 V 求解全部收敛；0、1、2、5、10 mV 严格低偏压
也全部收敛。固定异常边 3457 的同算子 SG 误差由约 11.9 dex 降至
0.166--0.170 dex，说明修复命中了接触附近第一处错误的统计/输运语义。

## 2. 依据

Sentaurus BVmethods 原始输入
`build-release/reference_tcad/bvmethods_sentaurus2018/run01/raw/pp4_des.cmd`
同时启用 `Fermi`、`EffectiveIntrinsicDensity(OldSlotboom)` 和
`HighFieldsaturation(GradQuasiFermi)`；对应 `n4_des.log` 明确记录
`Fermi Statistic` 与 `OldSlotboom with bandgap narrowing (Fermi)`。

外部资料交叉核对：

- [Sentaurus Device 入门：Electrode/Ohmic 接触](https://ghzphy.github.io/Sentaurus_Training/sd/sd_1.html)
  说明默认理想 Ohmic 接触采用热平衡和电中性条件；
- [Sentaurus Device 物理模型设置](https://ghzphy.github.io/Sentaurus_Training/sd/sd_2.html)
  给出 `Fermi` 与 `EffectiveIntrinsicDensity(OldSlotboom)` 的 deck 语义；
- [QTCAD Poisson 理论文档](https://docs.nanoacademic.com/qtcad/theory_spin_fem/poisson/)
  给出 `n=Nc F_{1/2}((EF-Ec)/kT)`、归一化定义和 Bednarczyk 近似，并将
  Ohmic 边界表述为非线性电中性问题；
- [Bessemoulin-Chatard, A finite volume scheme for convection-diffusion equations](https://arxiv.org/abs/1011.2299)
  给出非线性扩散/退化统计下的广义 Scharfetter-Gummel 构造。

本地开源实现对照：

- `D:/code-repo/Genius-TCAD-Open/include/math/mathfunc.h`：Bednarczyk
  `F_{1/2}`、导数和反函数；
- `D:/code-repo/Genius-TCAD-Open/src/solver/ddm1/ddm1_boundary_ohmic.cc`：
  Fermi 电中性 Ohmic 方程；
- `D:/code-repo/tcad-charon/src/evaluators/Charon_FermiDirac_Integral_*`：
  Fermi 积分与反函数；
- `D:/code-repo/tcad-charon/src/bc_strategies/Charon_BC_OhmicContact_impl.hpp`：
  Ohmic 平衡状态；
- `D:/code-repo/tcad-charon/src/evaluators/Charon_Degeneracy_Factor_impl.hpp`
  与 `Charon_DiffCoeff_Default_impl.hpp`：退化因子和广义 Einstein 关系。

本轮尝试通过 SSH 读取
`/usr/synopsys/sentaurus/O_2018.06-SP2/tcad/O-2018.06-SP2/Applications_Library/GettingStarted/sdevice/BVmethods`，
但虚拟机 `192.168.119.128:22` 超时；因此执行对照使用此前从同一 O-2018.06
项目复制到仓库的原始 cmd/par/log/tdr 产物。

## 3. 实现公式

Vela 的 intrinsic-potential 约定下：

```text
eta_n = (psi - phin)/Vt + ln(ni_eff/Nc)
eta_p = (phip - psi)/Vt + ln(ni_eff/Nv)
n = Nc F_1/2(eta_n)
p = Nv F_1/2(eta_p)
```

非退化极限 `F_{1/2}(eta) -> exp(eta)`，因而严格退化回原有 Boltzmann
公式。Ohmic 接触求解：

```text
g(psi) = n(psi, phin=0) - p(psi, phip=0) - Nnet = 0
```

使用括区间的 Newton/bisection 混合迭代。广义 SG 的边有效热电压为：

```text
Vt* = Vt (eta_1-eta_0)/(ln(c_1)-ln(c_0))
```

相等端点采用 `F_{1/2}/F'_{1/2}` 极限；平坦准费米势有显式零通量短路，
避免大数相消。

## 4. 配置

```json
{
  "solver": {
    "method": "newton",
    "carrier_statistics": "fermi_dirac"
  }
}
```

BVmethods 辅助脚本支持：

```text
--carrier-statistics fermi_dirac
```

低偏压脚本还新增 `--node-doping-file`，保证平衡态与偏压扫描使用同一份
sorted-node-order 修正掺杂。

## 5. BVmethods 数值结果

产物：

- 0 V 分层：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/fermi_dirac_20260802/equilibrium_stages`
- 严格低偏压：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/fermi_dirac_20260802/low_bias_strict`
- 逐边对照：
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/fermi_dirac_20260802/low_bias_strict_edge_compare`

五个 0 V 阶段 `constant -> SRH -> OldSlotboom -> Masetti -> Masetti high-field`
均收敛。低偏压端电流为：

| Vd (V) | Vela F-D (A/um) | Sentaurus (A/um) | Vela/S | error (dex) |
|---:|---:|---:|---:|---:|
| 0.001 | 8.716406e-11 | 1.028661e-10 | 0.847355 | -0.071935 |
| 0.002 | 1.706412e-10 | 2.015401e-10 | 0.846686 | -0.072278 |
| 0.005 | 4.000987e-10 | 4.737424e-10 | 0.844549 | -0.073375 |
| 0.010 | 7.199892e-10 | 8.563496e-10 | 0.840765 | -0.075325 |

整体低场电导与 Boltzmann 基线只变化约 0.14%，没有用统计模型掩盖迁移率
差异。局部 edge 3457 则发生决定性改善：

| Vd (V) | QF drop error old (dex) | QF drop error F-D (dex) | SG error old (dex) | SG error F-D (dex) |
|---:|---:|---:|---:|---:|
| 0.001 | 1.5628 | -0.1086 | 11.8864 | 0.1698 |
| 0.002 | 1.5555 | -0.1090 | 11.8896 | 0.1694 |
| 0.005 | 1.5332 | -0.1101 | 11.8985 | 0.1683 |
| 0.010 | 1.4950 | -0.1120 | 11.9120 | 0.1664 |

edge 3457 已不再属于异常边；1 mV 时其 Vela/Sentaurus 准费米势降分别为
`-2.98675e-12 V` 与 `-3.83562e-12 V`。

## 6. 尚未闭合的差异

0 V 接触节点 914：

| quantity | Boltzmann Vela | F-D Vela | Sentaurus |
|---|---:|---:|---:|
| psi (V) | 0.5527992 | 0.6458991 | 0.5884232 |
| n (cm^-3) | 3.256917e20 | 3.256917e20 | 3.256917e20 |
| p (cm^-3) | 87.01 | 2.374 | 190.44 |

F-D 已保证多数载流子、电中性和边输运内部一致，但绝对接触势/少数载流子
仍未复现 Sentaurus。这把下一问题收窄到 Sentaurus 的 DOS、能带参考和
OldSlotboom(Fermi) 带边分配映射，而不是 Ohmic 经验电压或迁移率。下一步应
从 O-2018.06 `models.par` 的 `eDOSMass/hDOSMass`、`Bandgap.Bgn2Chi`、
`dEg0(OldSlotboom)` 联合重放节点 914 的 `n,p,psi`，再决定 Vela 是否需要显式
的 `DeltaEc/DeltaEv`，而不能继续只用单一 `ni_eff` 表示 BGN。

## 7. 测试

- `test_carrier_statistics`: Fermi 积分、导数、反函数、密度/QF 互逆、
  高掺杂 Ohmic 中性、平坦 QF 零通量；
- `test_newton_solver`: coupled residual 与 F-D 端电流广义 SG 一致；
- `test_sg_flux`, `test_dd_gummel`: Boltzmann/SG/Gummel 兼容回归。
