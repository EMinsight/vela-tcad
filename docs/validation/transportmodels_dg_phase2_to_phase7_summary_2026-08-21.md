# TransportModels DG 阶段 2–7 验证总结

日期：2026-08-21
对比基准：Sentaurus 2022 TransportModels 电子 DG MOS 算例
工作点：空间与单点审计使用 `Vg=1 V, Vd=2 V`；曲线回归包含 21 个 Id–Vg 点和 21 个 Id–Vd 点。

## 结论

阶段 2–7 已全部执行并形成可复核的数据、图表和哈希报告。经材料带边驱动修正、`sentaurus_box` 离散、连续中性界面条件以及现有隐式 Lombardi 迁移率组合后，完整 42 点回归全部收敛。

该组合显著改善了 Id–Vd 和 2 V 端点：Id–Vd 最大相对误差为 `1.081%`，中位误差为 `0.576%`，2 V 端点误差为 `0.960%`，均达到预设门槛。但 Id–Vg 仍未完全达标：过渡区最大对数误差为 `0.1780 dex`（门槛 `0.15 dex`），导通区最大相对误差为 `14.942%`（门槛 `10%`）。因此最终状态为：**执行完成，验收部分通过**。

## 阶段结果

| 阶段 | 任务 | 主要结果 | 决策 | 状态 |
|---:|---|---|---|---|
| 2 | DG 固定状态残差审计 | 自由节点最大原始残差 `21010.7`；全局 L1 中反应项占 `73.53%`，基板区占 `69.04%` | 误差主因位于体方程/材料驱动，不是显式界面项 | 完成 |
| 3 | DG 参数、质量、单位与带边驱动核对 | TDR 材料契约使固定状态残差 L1 降低 `41.67%`；Si/PolySi 亲和势需 `+22.740 mV`，SiO2 需 `-50.000 mV`；BGN 份额保持 `0.5` | 采用 Sentaurus 2022 材料带边契约，不用非物理 BGN 拟合 | 完成 |
| 4 | Si/SiO2 界面条件比较 | 连续中性条件 L1=`405275`；半跳变与旧 SingleDevice 仿射条件分别高 `0.13%`、`0.27%` | 采用连续中性界面；旧网格的仿射系数不迁移 | 完成 |
| 5 | DG 离散与自洽验证 | `sentaurus_box` 自洽收敛，Id=`7.122957e-4 A/um`，2 V 误差 `0.9596%`；P1 direct 控制未收敛 | 采用 `sentaurus_box`；不把不可比的原始算子残差当作电压误差 | 完成 |
| 6 | 表面迁移率空间审计 | Frozen-Q 下现有隐式 Lombardi 电流误差 `1.3156%`，优于显式沟道 `2.0399%`、关闭 Enormal `17.9055%`、关闭高场饱和 `170.800%` | 保留现有隐式 Lombardi，避免单端点拟合 | 完成 |
| 7 | 完整 Id–Vg/Id–Vd 回归 | 42/42 点全部收敛；Id–Vd 达标；Id–Vg 两项门槛未达标 | 固定当前组合为新的诊断基线，下一步只针对 Id–Vg 剩余差异 | 完成，部分通过 |

## 阶段 7 验收矩阵

| 指标 | 结果 | 门槛 | 判定 |
|---|---:|---:|---|
| DG Id–Vd 最大相对误差 | `1.0809%` | `≤5%` | 通过 |
| DG Id–Vd 2 V 端点相对误差 | `0.9596%` | `≤3%` | 通过 |
| DG Id–Vg 过渡区最大对数误差 | `0.1780 dex` | `≤0.15 dex` | 未通过 |
| DG Id–Vg 导通区最大相对误差 | `14.9418%` | `≤10%` | 未通过 |

深关断区前三点的最大误差为 `8.4005 dex`。该区 Sentaurus 电流约停留在 `10^-15 A/um` 量级，而 Vela 继续下降到约 `10^-23 A/um`；它主要反映电流下限/极低泄漏模型差异，不参与本轮过渡区和导通区门槛。

## 严格端点空间场验证

空间比较固定使用阶段 5 的严格 `1 mV` 外层收敛状态，而不是反向 Id–Vd 扫描最终的 `Vd=0 V` 状态。对 Si/SiO2 表面 59 个节点：

| 物理量 | 中位绝对误差 | p95 | 最大值 |
|---|---:|---:|---:|
| 电子量子势 | `11.935 mV` | `14.786 mV` | `15.537 mV` |
| 电子浓度对数误差 | `0.00534 dex` | `0.08288 dex` | `0.17746 dex` |

量子势 p95 已低于原定 `20 mV` 空间门槛；电子浓度在绝大多数表面节点也保持良好一致。剩余 Id–Vg 误差更可能来自偏压相关的量子势/迁移率耦合、阈值附近的载流子统计或极低电流处理，而不是 2 V 单端点的整体 DG 离散错误。

## 后续建议

1. 冻结当前材料、`sentaurus_box`、连续中性界面和隐式 Lombardi 组合，不再同时改变多项模型。
2. 只对 Id–Vg `-0.4 V` 至 `0.2 V` 的过渡区逐点保存 Qn、电子浓度、Enormal、迁移率和准费米势，做偏压相关误差归因。
3. 对导通区最大误差点做 Frozen-Q 与 Frozen-mobility 双向试验，区分 DG 偏差和表面迁移率偏差。
4. 单独定义深关断电流下限或泄漏机制的验收规则，避免把数值地板差异混入主要输运模型验收。

## 证据与复现入口

- [阶段 2 固定状态残差审计](transportmodels_dg_fixed_state_residual_audit_2026-08-20.md)
- [阶段 3 参数与单位审计](transportmodels_dg_parameter_sweep_2026-08-21.md)
- [阶段 3 带边驱动审计](transportmodels_dg_band_drive_audit_2026-08-21.md)
- [阶段 4 界面条件比较](transportmodels_dg_interface_sweep_2026-08-21.md)
- [阶段 5 离散格式审计](transportmodels_dg_discretization_audit_2026-08-21.md)
- [阶段 5 自洽端点验证](transportmodels_dg_discretization_self_consistent_2026-08-21.md)
- [阶段 6 表面迁移率审计](transportmodels_dg_surface_mobility_audit_2026-08-21.md)
- [阶段 7 完整回归报告](transportmodels_dg_phase7_regression_2026-08-21.md)

复现脚本为 `scripts/run_transportmodels_dg_fixed_state_residual_audit.py` 至 `scripts/run_transportmodels_dg_phase7_regression.py`。阶段 7 脚本会显式把 `D:\msys64\ucrt64\bin` 加入子进程 PATH，以保证 `libgcc_s_seh-1.dll`、`libstdc++-6.dll` 和 `libwinpthread-1.dll` 可被加载。
