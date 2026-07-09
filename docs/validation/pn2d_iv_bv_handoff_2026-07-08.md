# PN2D IV/BV Debug Handoff - 2026-07-08

## 目的

这份文档用于在另一台电脑继续当前 PN2D IV/BV 调试任务。当前分支：

`codex-pn2d-sentaurus2018-calibration`

当前任务主线：

- PN2D IV：继续定位正向 IV 电流比旧基线/Sentaurus 高 1-2 个数量级的问题。
- PN2D BV：继续定位雪崩击穿源项偏小的问题；本轮已修掉 SG avalanche 场强/alpha 单位缩放错误，剩余问题集中到 current proxy/continuation。

## 本轮已完成的代码修复

### 1. SG edge avalanche helper 接入 fieldFactor

修改文件：

- `include/vela/equation/AssemblerUtils.h`
- `src/equation/DDAssembler.cpp`
- `src/equation/CoupledDDAssembler.cpp`
- `src/solver/GummelSolver.cpp`
- `src/simulation/DCSweep.cpp`

核心变化：

- `detail::sgEdgeCurrentAvalancheSourceRecords(...)`
- `detail::sgEdgeCurrentAvalancheSourceComponentIntegrals(...)`
- `detail::sgEdgeCurrentAvalancheSourceIntegrals(...)`

都新增可选参数：

```cpp
Real fieldFactor = 1.0
```

并在以下量上使用：

- edge electric field
- electron/hole quasi-Fermi driving field
- current-aligned signed electric field
- SG continuity flux coefficient

实际求解器调用方已传入：

```cpp
scaling.unitSystem().fieldFromCoordinateDeltaFactor()
```

这使 BV 残差源项和诊断使用同一套 unit-scaling 场强。

### 2. SG avalanche edge 诊断 CSV 改为物理单位

修改文件：

- `src/simulation/DCSweep.cpp`
- `tests/test_dc_sweep.cpp`

修复前，`sg_avalanche_edges.csv` 的列名写着 `*_V_per_m`、`*_m_inv`，但实际部分列仍是 TCAD internal 单位。

修复后：

- `x0_um/y0_um/x1_um/y1_um` 输出物理 um。
- `edge_length_m` 输出物理 m。
- `edge_area_proxy_m2` 输出物理 m^2。
- `electric_field_V_per_m`、`electron_impact_field_V_per_m`、`hole_impact_field_V_per_m` 输出物理 V/m。
- `electron_alpha_m_inv`、`hole_alpha_m_inv` 输出物理 1/m。
- `electron_mobility_m2_V_s`、`hole_mobility_m2_V_s` 输出物理 m^2/(V*s)。

## 新增测试

### 1. Impact-ionization helper 缩放测试

文件：

`tests/test_impact_ionization.cpp`

新增 case：

```text
SG edge avalanche source records apply coordinate field scaling
```

覆盖：

- `fieldFactor=1e6` 时，edge electric field、electron/hole impact field、electron/hole raw flux proxy 按比例放大。

### 2. DCSweep SG edge 诊断单位测试

文件：

`tests/test_dc_sweep.cpp`

增强 case：

```text
DCSweep: SG avalanche edge diagnostics write assembled source rows
```

新增断言：

- `sg_avalanche_edges.csv` 中最大 `electric_field_V_per_m` 应等于 `result.points.front().maxElectricField * 100.0`。
- 这里 `result.points.front().maxElectricField` 是 internal electric field，`*100` 转为物理 V/m。

## 已验证命令

在 `D:\code-repo\vela-tcad`，PowerShell 下先设置 MSYS2 UCRT64：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
```

本轮已跑过：

```powershell
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
git diff --check
```

结果：

- build 通过。
- `ctest` 结果：`437/437` 通过。
- `git diff --check` 通过。

## BV 当前定位结论

### 已修复的问题：SG avalanche 场强/alpha 缩放

复跑粗网格配置：

```powershell
.\build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_previous_full20_sentaurus_math_20260708\simulation_coarse_previous_full20.json
```

关键输出路径：

- `build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_previous_full20_sentaurus_math_20260708\coarse_previous_full20.csv`
- `build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_previous_full20_sentaurus_math_20260708\sg_avalanche_edges.csv`

修复后的关键数值，最终点 `-20 V`：

- 主 CSV `max_electric_field_V_per_m = 4.174020311437548e7`
- SG edge CSV 最大 `electric_field_V_per_m = 4.174020311437548e7`
- SG edge 最大 `electron_alpha_m_inv ~= 3.302963e6`

结论：

- 之前 SG edge 诊断中最大场强只有约 `4.17e1` 或 internal/物理单位混用的假象，已经修正。
- avalanche alpha 不再是 `~1e-298` 的假零。

### 仍未解决的问题：coarse `cell_reconstructed` source proxy 压塌

同一粗网格 `-20 V` 最终点仍然：

- `breakdown_detected = 0`
- `sg_sum_source_integral ~= 1.525e-62`
- 最大 source edge 的 `electron_final_over_raw_flux_proxy ~= 3.05e-80`
- 最大 source edge 的 `hole_final_over_raw_flux_proxy ~= 1.16e-81`

这说明剩余 BV 源项偏小不是场强/alpha 了，而是 `current_approximation = cell_reconstructed` 的 flux proxy 被重构中点密度压塌。

对照实验：

将同一粗网格配置临时改为：

```json
"current_approximation": "density_gradient"
```

配置输出：

`simulation_coarse_previous_full20_density_gradient_probe.json`

该实验未跑到 `-20 V`，但最后收敛点已经足够说明问题：

- last converged bias: `-5.24321892849 V`
- `sg_sum_source_integral ~= 6.54908394089131e13`
- 最大 source edge `electron_flux_proxy == electron_raw_flux_proxy`
- `electron_final_over_raw_flux_proxy = 1`
- 后续在约 `-5.39331658584 V` 因 `line_search_non_decrease` 失败。

结论：

- `density_gradient` 能释放 avalanche source，但强反馈使 continuation 提前失败。
- `cell_reconstructed` 更稳定但会把 flux proxy 压到近零。

### BV 下一步建议

优先级最高：

1. 不要再追 SG field/alpha 单位；它已经修正。
2. 继续比较 `cell_reconstructed`、`density_gradient`、`grad_qf`、`conserved_total_current` 对源项和 continuation 的影响。
3. 对 `cell_reconstructed` 检查 `cellReconstructedAvalancheMidpointDensity(...)` 的 Bernoulli midpoint 逻辑，在反偏耗尽区是否过度选择极低载流子密度。
4. 可考虑新的 guarded proxy：保留 `density_gradient` raw SG flux 量级，但对 avalanche feedback 做 continuation damping / source limiting，而不是用 midpoint density 直接压没源项。
5. 继续遵守已有记忆里的 gate：到达 `-20 V` 或 Newton residual 收敛不代表恢复击穿分支；必须检查 multiplication-current/source order。

## IV 当前定位结论

当前 IV 诊断配置：

`build-release\reference_tcad\pn2d_sentaurus2018\reports\iv_rerun_20260708\simulation_iv_current_gap_diag.json`

复跑命令：

```powershell
.\build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\iv_rerun_20260708\simulation_iv_current_gap_diag.json
```

关键输出：

- `pn2d_iv_current_gap_diag.csv`
- `pn2d_iv_current_gap_diag_terminal_compare.csv`
- `pn2d_iv_current_gap_diag_contact_edge.csv`
- `pn2d_iv_current_gap_diag_newton_history.csv`

关键数值：

0.3 V：

- `current_total_A_per_um = -2.4084716164176403e-12`
- `current_electron_A_per_um = -2.399772232218817e-12`
- `current_hole_A_per_um = 8.6993841988234887e-15`
- `newton_convergence_reason = poisson_line_search_stall_floor`

终端电流方法对照，0.3 V Cathode：

- `I_sgflux_A_per_um = -2.4084716164176403e-12`
- `I_residual_A_per_um = -2.4084716164176403e-12`
- `I_sgflux_with_qf_floor_A_per_um = -2.4084716164176403e-12`
- `sg_avalanche_source_integral_total = 0`

结论：

- IV 偏差不是 terminal current extraction 问题。
- 不是 BV avalanche 代码造成的，因为 IV avalanche source 为 0。
- SG flux 和 residual current 一致，说明偏差在求解状态/接触边离散/边界重构。

### 已排除的 IV 假设

1. `masetti` vs `masetti_field`

临时配置：

`simulation_iv_current_gap_diag_masetti_field.json`

结果：

- 0.1 V：约 `0.9995x`
- 0.2 V：约 `0.9971x`
- 0.3 V：约 `0.9808x`

结论：迁移率模型差异只有约 2%，不是数量级偏差主因。

2. 严格 Poisson line-search floor

临时配置：

`simulation_iv_current_gap_diag_strict_floor.json`

设置：

- `stall_residual_floor = 1e-10`
- `poisson_line_search_stall_residual_floor = 1e-10`
- `poisson_line_search_stall_carrier_residual_floor = 1e-10`

结果：

- 只收敛 0 V。
- 在 `0.000625 V` 因 `line_search_non_decrease` 失败。

结论：当前严格 Newton 不能简单通过降低 floor 进入 0.3 V；之前的 floor 是稳定性门槛，但不直接解释高电流。

3. 终端电流抽取

0.3 V 下 SG flux、residual current、QF-floor current 完全一致。

结论：不是后处理抽取路径错误。

### IV 接触边定位

0.3 V 下 `pn2d_iv_current_gap_diag_contact_edge.csv` 显示：

- 高电流主要来自 Cathode 接触边多数电子项。
- 典型边 `edge_id=5579`：
  - `current_electron ~= -1.499857645e-7` internal
  - `current_total_A_per_um ~= -1.505294792e-13`
  - 多条类似接触边累加得到总 `-2.408e-12 A/um`
  - `phip0 ~= 0.1952894 V`, `phip1 = 0`
  - `n0/n1 ~= 1e17`

下一步应检查：

1. Cathode 接触边的 contact/internal node 关系是否符合 Sentaurus 的 ohmic 边界处理。
2. `contact_boundary_reconstruction = dominant_signed_contact_mean` 对 IV 低/中正偏是否引入接触边多数载流子漂移项。
3. `ContactCurrent.cpp` 和 `CoupledDDAssembler` 在接触边 SG flux 约定是否与 residual/current extraction 一致但物理边界过强。
4. 对比旧 stale baseline `build-release\reference_tcad\pn2d_sentaurus2018\vela\pn2d_sentaurus2018_iv.csv`，旧文件在 0.2/0.3 V 与 Sentaurus 很接近，但当前 rerun 高很多；需要用 git history 或临时 worktree 找到改变电流状态的提交。

## 回家后建议执行顺序

1. 更新/切到本地分支：

```powershell
Set-Location "D:\code-repo\vela-tcad"
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git status --short
git log --oneline -5
```

2. 先确认提交中的测试仍然通过：

```powershell
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

3. 若继续 BV：

先不要再改 field scaling。优先做 current proxy matrix：

- `cell_reconstructed`
- `density_gradient`
- `grad_qf`
- `conserved_total_current`

比较每种的：

- last converged bias
- `sg_sum_source_integral`
- `electron/hole_final_over_raw_flux_proxy`
- terminal current order
- failure reason

4. 若继续 IV：

优先做 contact-edge/old-baseline bisect：

- 固定 `simulation_iv_current_gap_diag.json`
- 对比旧 baseline 与当前 rerun 的 0.1/0.2/0.3 V
- 使用 git history 或临时 worktree 找“接触边电流从 Sentaurus 量级跳到当前量级”的提交
- 重点看 `ContactCurrent.cpp`、`CoupledDDAssembler.cpp`、接触边界重构、carrier row recovery、Newton handoff 相关提交

## 当前未提交前的工作树范围

预期修改文件：

- `include/vela/equation/AssemblerUtils.h`
- `src/equation/CoupledDDAssembler.cpp`
- `src/equation/DDAssembler.cpp`
- `src/simulation/DCSweep.cpp`
- `src/solver/GummelSolver.cpp`
- `tests/test_dc_sweep.cpp`
- `tests/test_impact_ionization.cpp`
- 本文档

临时仿真输出位于 `build-release/...`，一般被忽略，不需要提交。

## 2026-07-09 续查进展

### BV：cell_reconstructed midpoint 权重缺陷

新增定位发现 `cell_reconstructed` 源项压塌的直接原因在 `detail::avalancheMidpointAux2(x)` 的负分支。该函数文档语义是 `1/(1+exp(x))`，但旧实现对 `x < 0` 返回 `exp(x)/(1+exp(x))`，导致大负 potential drop 下两个 Bernoulli midpoint 权重都趋近 0，载流子中点密度和 avalanche flux proxy 被压低约 80 个数量级。

已修改：

- `include/vela/equation/AssemblerUtils.h`
- `tests/test_cell_reconstructed_avalanche.cpp`

新增测试：

```text
Bernoulli avalanche midpoint weights stay normalized at large potential drops
```

阶段性复跑 coarse7x3 proxy matrix 曾显示：

- `cell_reconstructed` 不再以假低源项跑到 `-20 V`。
- last converged bias 从旧的 `-20 V` 变为 `-5.24321892849 V`。
- `sg_sum_source_integral` 从 `1.525e-62` 提升到 `5.557e13`。
- 最大源项边的 `electron_final_over_raw_flux_proxy` 从 `3.05e-80` 提升到 `0.848`。

### IV：Poisson stall 接触边 QF drop 门控

继续定位后确认，旧基线分叉点被夹在以下提交之间：

- `de7b8ba Enable reference TCAD sweep initialization`：完全复现旧 IV 基线，0.3 V 为 `-6.9321746619913786e-14 A/um`，`reltol` 收敛。
- `4fd40bd Implement TCAD internal unit system`：只到 0 V，0.000625 V 因 `line_search_non_decrease` 失败。
- `b12f342 Stabilize PN2D scaling and Newton stalls`：跑到 0.3 V，但 0.3 V 为 `-2.4084716164176403e-12 A/um`，由 `poisson_line_search_stall_floor` 放行。

同一 Cathode 主导边 `edge_id=5608` 在 0.3 V 的对比：

- `de7b8ba`：`current_total_A_per_um ~= -4.3326e-15`，`current_electron ~= -3.0203e-09`，`phin0 ~= 2.59e-12`。
- b12/current-stall path：`current_total_A_per_um ~= -1.5053e-13`，`current_electron ~= -1.4999e-07`，`phin0 ~= 1.29e-10`。

本次新增 `poisson_line_search_stall_contact_majority_qf_drop_limit_V`，默认 `5e-11 V`，只约束 `poisson_line_search_stall_floor` 接受路径。当接触边最大多数载流子 quasi-Fermi drop 超过该值时，不再把 Poisson residual floor/line-search 停滞声明为收敛。

### 2026-07-09 后续修复：unit_scaling coupled Poisson residual

严格 floor deck 的 0.000625 V 失败诊断显示：

- `failure_reason = line_search_non_decrease`
- `residual_norm = 2.9337775388330295e-07`
- block residuals：`psi = 2.9337644336137227e-07`，`phin = 8.768996355699269e-10`，`phip = 2.537923058723764e-10`
- `max_contact_majority_qf_drop_V = 1.3444649039906498e-14`

这说明微小偏置停滞不是接触 QF drop 放行问题，而是 `unit_scaling` 下 coupled Newton 的 Poisson residual 尺度问题。新增红灯测试 `CoupledDDAssembler unit scaling Poisson residual matches DDAssembler` 后，修复前同一状态下 coupled Poisson residual 为 `-7.976646082218936e8`，而 `DDAssembler` nonlinear residual 为 `0.41805291134105771`。

修复内容：

- `src/equation/CoupledDDAssembler.cpp`：Poisson residual 的体电荷项乘以 `chargeVolumeFactor`。
- `src/equation/CoupledDDAssembler.cpp`：Poisson Jacobian 中载流子/势导数同步乘以 `chargeVolumeFactor`。
- `tests/test_sg_flux.cpp`：新增 coupled/DD Poisson residual 等价测试。

复跑结果：

- `simulation_iv_current_gap_diag_strict_floor.json`：`converged=true`，16 points，0.3 V `current_total_A_per_um = -1.0777065624368375e-10`，`newton_convergence_reason = reltol`。
- `simulation_iv_current_gap_diag.json`：`converged=true`，16 points，0.3 V 同为 `-1.0777065624368375e-10 A/um`，`reltol`。
- Sentaurus reference 0.3 V 为 `+6.93205320128e-14 A/um`，旧 Vela baseline 为 `-6.9321746619913786e-14 A/um`。因此 line-search/non-decrease 已解决，但 IV 电流量级仍偏高约 `1.55e3x`，后续应继续定位 transport/current scaling 或已收敛状态差异。
- 0.3 V terminal current direct/residual/mixed 三种提取一致：Cathode `-1.0777065624368375e-10 A/um`，Anode `+1.0777065635302712e-10 A/um`，说明当前偏差不是终端电流后处理方法之间的不一致。
- `simulation_cell_reconstructed.json`：coarse7x3 BV 现在 `converged=true`，51 points，到达 `-20 V`；末点 `current_total_A_per_um = -3.2326146014545371e-12`，`max_electric_field_V_per_cm = 104068.09472457141`，`carrier_product_max_np_over_ni2 = 3.1262979253679117`。

当前结论：

- 微小偏置严格 Newton 的 `line_search_non_decrease` 根因已收窄并修复为 coupled Poisson `unit_scaling` 体电荷/Jacobian 漏乘 `chargeVolumeFactor`。
- BV 的 `cell_reconstructed` 源项压塌和 continuation 早停均已有明显改善，当前 coarse 算例可跑到 `-20 V`，但击穿判据仍未触发，雪崩强度/电流增长还需继续和 Sentaurus 对齐。
- IV 已从“阻止坏状态被接受”推进到“严格 Newton 可收敛”，但收敛解仍不在 Sentaurus/旧基线电流量级；下一步应定位 SG transport/current scaling、contact-edge drift/diffusion cancellation、以及 current internal-unit conversion。

## 2026-07-09 提交前最新进展：0.3 V 接触边状态级对齐

本轮继续执行“对 0.3 V 主导接触边做状态级对齐”的计划。当前 evidence 已基本排除基础 terminal current scaling 漏系数，问题进一步收窄到当前收敛解中的接触边 electron quasi-Fermi drop / drift-diffusion cancellation。

### 复跑和对齐数据

旧基线 worktree：

```text
C:\tmp\vela-iv-bisect
commit de7b8ba Enable reference TCAD sweep initialization
```

旧基线使用当前 `unit_scaling` deck 复跑，并把输出写入：

```text
build-release\reference_tcad\pn2d_sentaurus2018\reports\iv_state_alignment_20260709\
```

关键产物：

- `old_de7b8ba_iv.csv`
- `old_de7b8ba_iv_contact_edge.csv`
- `pn2d_iv_0p3_curve_state_alignment_summary.csv`
- `pn2d_iv_0p3_contact_edge_state_alignment_paired.csv`
- `pn2d_iv_0p3_contact_edge_state_alignment_physical_normalized.csv`
- `pn2d_iv_0p3_contact_edge_state_alignment_physical_normalized.md`

注意：这些位于 `build-release/...` 的分析产物是本地临时输出，一般不提交。回家后若 build 目录未同步，需要重新运行对应 deck 或从当前机器拷贝。

### 曲线级结论

0.3 V 电流：

- Sentaurus reference curve：`+6.93205320128e-14 A/um`
- 旧 Vela `de7b8ba`：`-6.9321746619913786e-14 A/um`，与 Sentaurus magnitude 比值 `1.000018`
- 当前分支 `11b3c04`：`-1.0777065624368375e-10 A/um`，比 Sentaurus/旧基线高约 `1554.6x`

本地未找到可直接用于 0.3 V 状态场逐点比较的 Sentaurus field export，因此 Sentaurus 当前只作为曲线电流锚点；状态级比较是旧 Vela 基线 vs 当前 Vela 分支。

### 接触边状态级结论

matched cathode contact edges 上的代表性结果：

| edge_id | old_rank | current_rank | edge current ratio | electron QF drop ratio | normalized n0 ratio | normalized n1 ratio | normalized mun ratio | cancellation old/current |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5608 | 1 | 15 | `1554.519x` | `2612.0x` | `0.733x` | `1.000x` | `1.000x` | `2.98 / 2.407e6` |
| 5435 | 15 | 1 | `1554.770x` | `2612.4x` | `0.733x` | `1.000x` | `1.000x` | `2.98 / 2.407e6` |
| 5455 | 8 | 8 | `1554.644x` | `2612.2x` | `0.733x` | `1.000x` | `1.000x` | `2.98 / 2.407e6` |

归一化说明：

- 旧基线 diagnostic 的 `n*_cm3` 列名沿用了 cm^-3，但旧 `unit_scaling` 路径实际写出 SI m^-3 数值；对比时需除以 `1e6`。
- 旧基线 diagnostic 的 `mun_cm2_V_s` 列名沿用了 cm2/V/s，但实际写出 SI m2/V/s 数值；对比时需乘以 `1e4`。
- 归一化后，majority electron mobility 完全一致，contact-side `n1` 完全一致，`n0` 仅为旧基线的约 `0.733x`。

因此当前 IV 电流 `1554x` 偏高不是 mobility/density 物理单位整体错位，也不是 `ContactCurrent` 后处理抽取路径问题。更强 evidence 指向：

- 当前收敛状态在 cathode contact edge 上的 `phin0 - phin1` 比旧基线大约 `2.6e3x`。
- 当前解的 electron drift/diffusion 大数抵消极强，cancellation ratio 从旧基线约 `3` 增至约 `2.4e6`。
- direct SG flux、residual current、QF-floor current 三种 terminal current 仍完全一致，所以偏差在求解状态/接触边 quasi-Fermi 边界附近，而不是后处理电流方法之间的不一致。

### 本次代码修改范围

本次提交前的主要源码/测试修改：

- `src/equation/CoupledDDAssembler.cpp`：修正 coupled Poisson residual/Jacobian 的 `chargeVolumeFactor` 缩放，并补齐 unit-scaling continuity/current 相关一致性。
- `src/equation/DDAssembler.cpp`：同步 SG avalanche/current helper 的 scaling 参数。
- `include/vela/equation/AssemblerUtils.h`：修正 avalanche midpoint 权重和 SG avalanche helper 缩放接口。
- `src/post/ContactCurrent.cpp`：修正 contact edge diagnostic 的 `edge_length_m`、`edge_couple_m` 物理单位输出。
- `src/solver/NewtonSolver.cpp`、`include/vela/solver/NewtonSolver.h`、`src/simulation/DCSweep.cpp`：补充 strict Newton/Poisson stall/QF drop guard 相关诊断与控制。
- `tests/test_cell_reconstructed_avalanche.cpp`：覆盖 midpoint 权重归一化。
- `tests/test_newton_solver.cpp`：覆盖 Newton stall/floor 行为。
- `tests/test_sg_flux.cpp`：新增 coupled/DD Poisson residual 等价、continuity residual 等价、ContactCurrent legacy SI vs unit_scaling terminal current 等价测试。
- `docs/config_schema.md`：补充新增 solver 配置项。
- 本文档：更新 IV/BV 交接记录。

### 已执行验证

提交前本轮已执行：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
git diff --check
```

最近一次全量测试结果：`442/442` passed。随后新增的是文档交接内容；提交前仍需再跑一次 `git diff --check` 和必要的状态检查。

### 回家后继续执行建议

优先方向：不要再从 terminal current 后处理或 mobility/density 整体单位系数入手。下一步应直接查当前收敛解为什么在 contact edge 上产生过大的 electron quasi-Fermi drop。

建议顺序：

1. 先确认本地分支和提交：

```powershell
Set-Location "D:\code-repo\vela-tcad"
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git status --short
git log --oneline -5
```

2. 复跑当前 IV deck，确认 0.3 V 仍为当前状态：

```powershell
.\build-release\vela_example_runner.exe --config build-release\reference_tcad\pn2d_sentaurus2018\reports\iv_rerun_20260708\simulation_iv_current_gap_diag.json
```

预期当前 0.3 V：`current_total_A_per_um = -1.0777065624368375e-10`。

3. 对 contact edge 做更细状态比较：

- 固定 edge `5608`、`5435`、`5455`。
- 比较 old/current 的 `psi0/psi1`、`phin0/phin1`、`n0/n1`、Bernoulli argument、electron drift term、electron diffusion term。
- 重点看 `CoupledDDAssembler` 对 contact quasi-Fermi Dirichlet row 的设置是否与旧基线一致。
- 检查 carrier row recovery 是否在接触邻边把 `phin` 或 density 推到新的收敛分支。
- 检查 strict Newton 在低偏置 continuation 的初始状态继承是否让 contact-side `phin` drop 逐步积累。

4. 如果需要 Sentaurus 状态级锚点，需要在 Sentaurus 端导出 0.3 V field，包括：

- `Potential`
- `eQuasiFermiPotential`
- `hQuasiFermiPotential`
- `ElectronDensity`
- `HoleDensity`
- contact voltage/current scalar

当前仓库本地未发现可直接用于 0.3 V 状态逐点对齐的 Sentaurus field export。

5. BV 方向暂时作为第二优先级：

- coarse `cell_reconstructed` 可跑到 `-20 V`，但 `breakdown_detected` 仍未触发。
- 继续比较 `density_gradient`、`grad_qf`、`conserved_total_current` 的 source/current proxy 和 continuation 稳定性。
- 不要回退已修的 SG field/alpha/midpoint 权重问题。
