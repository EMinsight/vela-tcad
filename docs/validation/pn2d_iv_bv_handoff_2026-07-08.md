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
