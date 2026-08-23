# SLOT-LDMOS 直接增广负载线 Newton 验证（2026-08-22）

## 结论

Vela 已具备直接装配并求解设备—串联电阻增广系统的能力：

\[
\begin{bmatrix}
J & F_{V_\mathrm{in}} \\
R\,I_x & 1+R\,I_{V_\mathrm{in}}
\end{bmatrix}
\begin{bmatrix}\Delta x\\\Delta V_\mathrm{in}\end{bmatrix}
=-
\begin{bmatrix}
F(x,V_\mathrm{in})\\
V_\mathrm{in}+RI(x,V_\mathrm{in})-V_\mathrm{out}
\end{bmatrix}.
\]

该路径在解析折叠测试中能越过电压控制 Jacobian 奇异点，并在真实
SLOT-LDMOS 检查点上从外压 1187.82348632813 V 收敛到
1188.82348632813 V。验证点使用 `R=1e12 ohm*um`、自洽雪崩、
triangle GSS/local-AD 源 Jacobian 和 IALMob-off 迁移率配置。

## 实现与问题修复

| 项目 | 实现/修复 | 验证 |
|---|---|---|
| 直接增广矩阵 | 保留设备稀疏 Jacobian，并加入偏压列、终端电流导数行和负载线标量方程 | 折叠模型中 Schur 失败而直接增广收敛 |
| 终端电流导数 | 由接触连续性残差及原始 Jacobian 转置构造解析梯度 | 耦合 DC sweep 集成测试通过 |
| 稀疏后端 | Eigen SparseLU，失败时依次审计 UMFPACK、SPQR 最小二范数、宽松主元和受审计静态移位 | 每个候选均按原增广矩阵残差复核 |
| 数值缩放 | 行尺度同时考虑矩阵行最大系数与该行 RHS，随后做列均衡 | 消除原先 `rhs_inf=9.48e275` 的浮点放大；修复后 `rhs_inf<=1` |
| 两点预测 | 用两个已收敛负载线状态按外压割线预测，并复用设备更新限制器统一缩放状态与内压增量 | 真实检查点无需外压回退即收敛 |
| 非输运语义 | local-AD 雪崩单元循环跳过无 `ni`/迁移率的 SiO2、Nitride、PolySilicon 单元，保持绝缘体准费米单位 gauge | 最大伪更新由 `8.33e6` 降至 `0.156`；混合 Si/SiO2 回归通过 |
| 后端部署 | SuiteSparse DLL 仅在对应 pkg-config 目标存在时暂存 | SuiteSparse 开启和 `VELA_ENABLE_UMFPACK=OFF` 两种 CMake 配置均成功 |

## 真实检查点结果

| 外压 V | 内漏极电压 V | 漏极电流 A/um | 负载线残差 V | 增广 Newton 迭代 |
|---:|---:|---:|---:|---:|
| 1187.82348632813 | 15.7209505709223 | 1.1721025357572052e-9 | -3.410605131648481e-12 | 0 |
| 1188.82348632813 | 15.721531914507386 | 1.1731019544135857e-9 | -3.865352482534945e-11 | 41 |

第二点满足 `Vout = Vin + R*Id`，闭合误差远低于 `1e-6 V` 配置容差。

## 回归结果

| 测试 | 结果 |
|---|---:|
| `test_coupled_load_line` | 37 assertions / 5 cases 全部通过 |
| `test_dc_sweep [coupled_newton]` | 42 assertions / 1 case 全部通过 |
| `test_mos_mixed_material` | 1327 assertions / 8 cases 全部通过 |
| 真实 `simulation_local_ad_coupled_step.json` | `converged=true`, 2 points |
| `git diff --check` | 无 whitespace error（仅 Windows 行尾提示） |

## 尚未完成的范围

本次完成的是直接增广算法及一个真实高压步的闭环验证，不等同于完整
IALMob-off/on BVDS 曲线已经跑完。下一阶段应在固定相同外接电阻、步长策略
和检查点策略下分别完成 off/on 扫描，越过 `Id=1e-7 A/um`，再对内压交点
插值并报告对应外压和 IALMob 引起的 BVDS 偏移。
