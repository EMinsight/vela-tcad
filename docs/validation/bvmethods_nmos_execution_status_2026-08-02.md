# Sentaurus BVmethods NMOS：Vela 分层执行状态（2026-08-02）

## 结论

官方 Sentaurus 2018 BVmethods 六种方法已经完整运行并形成参考曲线；Vela
已经修通混合材料 0 V 平衡态，并完成
`constant -> SRH -> OldSlotboom -> Masetti -> Masetti high-field` 五层物理开启。

非零漏极偏压尚未通过物理门禁。Vela 在 10 mV 后处理算例中得到
`3.895659e-5 A/um`，而 Sentaurus IIC 是 `8.563450e-10 A/um`，比值
`4.54917e4`（`4.65793 decade`）。严格自洽雪崩在 0 V 收敛，但 10 mV
在五分钟内未收敛。因此没有继续生成 1/2/4/5/6 V 或 0--7 V 的伪 BV
曲线。

## 已完成项

### 1. 官方 Sentaurus 参考

参考目录：
`build-release/reference_tcad/bvmethods_sentaurus2018/run01`

提取结果：

| 方法 | BV (V) |
| --- | ---: |
| ABA Poisson | 5.305525633 |
| ABA coupled / IIC | 6.377494278 |
| external resistor | 6.379791636 |
| voltage-to-current | 6.383184201 |
| continuation | 6.383727169 |
| transient | 6.378835044 |

### 2. 混合材料修复

- Si/绝缘层界面节点的本征浓度现在优先使用可输运材料，避免氧化层单元排序
  把共享 Si 界面节点错误设为 `ni=0`。
- Newton 载流子合法性检查允许纯绝缘节点的合法 `n=p=0`，但仍要求所有
  半导体节点载流子严格为正且有限。
- Newton 失败诊断使用相同语义，不再把 pinned 绝缘节点计为
  `nonpositive carrier`。
- 线性系统非有限量和零行诊断会输出首批具体行列位置。

### 3. 0 V 分层物理结果

执行脚本：`scripts/run_bvmethods_nmos_equilibrium_stages.py`

结果目录：
`build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/equilibrium_stages`

| 阶段 | 收敛 | Newton 次数 | Id (A/um) | max(E) (V/cm) |
| --- | ---: | ---: | ---: | ---: |
| constant, no SRH, no BGN | 是 | 65 | 3.58489e-17 | 5.01484e6 |
| + SRH | 是 | 0 | 3.58213e-17 | 5.01484e6 |
| + OldSlotboom | 是 | 8 | -8.74944e-19 | 4.16437e6 |
| + Masetti | 是 | 6 | -8.82850e-22 | 4.16437e6 |
| + high-field QF-gradient | 是 | 0 | -8.84631e-22 | 4.16437e6 |

每层使用上一层的状态重启，并固定到同一冷启动残差尺度。否则重启状态会被
重新用作相对残差基准，造成同一物理误差被人为要求再下降八个数量级。

### 4. 雪崩门禁

执行脚本：`scripts/run_bvmethods_nmos_bias_validation.py`

- 严格 self-consistent 0 V：3 次 Newton，最终归一化残差
  `3.20998e-13`，终端电流约零。
- postprocess 10 mV：数值收敛且 QF 包络检查通过，但漏极电流比 Sentaurus
  高 `4.65793 decade`，判定为错误输运分支。
- self-consistent 10 mV：在五分钟限制内未收敛。电子连续性块主导，载流子
  保持有限且为正，不是原先的绝缘节点 `carrier_invalid` 问题。

## 当前阻塞点

1. 需要核对 Sentaurus `Barrier=-0.55` 与 Vela
   `flatband_voltage=-0.55` 的符号、参考能级和 metal-gate DD 边界语义。
   当前 Vela 会把 metal-gate 有效电势也写入绝缘节点的准费米边界；这些节点
   无输运，但该值会污染部分接触/QF 诊断。
2. 10 mV 漏极电流由约 `-0.415` 与 `+0.425 A/um` 的漂移/扩散项相消后得到，
   对平衡保持和 SG 电流精度非常敏感。应逐边对照 Sentaurus 的 QF、电流和
   密度场，优先定位 drain/source 接触邻边与沟道边。
3. Vela 尚未覆盖官方脚本中的 Fermi 统计、`Enormal` 表面迁移率和
   `Band2Band(E2)`。在确认边界和低偏压电流闭合前，不应通过调雪崩系数弥补
   这些上游差异。
4. Gummel 高级物理路径在节点 1172/1176 的非有限连续性项仍需独立定位；
   当前成功的 0 V 路径使用 coupled Newton。

## 下一执行顺序

1. 为 metal-gate 增加“不施加载流子准费米 Dirichlet”的显式边界语义和
   单元测试，同时让 QF bounds/接触 QF-drop 只检查输运接触。
2. 在 0、1、5、10 mV 导出 Vela/Sentaurus 的 `psi/phin/phip/n/p` 与接触邻边
   SG 通量，定位 `4.66 decade` 电流差的第一条边。
3. 修复低偏压分支后，重新执行 postprocess 1/2/4/5/6 V；只有电流、QF
   包络和连续性闭合同时通过，才计算 alpha、离化积分和积分雪崩电流。
4. 最后执行 self-consistent 0--7 V，并以 Sentaurus IIC
   `6.377494278 V` 为主要 BV 基准。

## 验证

- `test_newton_solver`: 1035 assertions / 70 cases passed
- `test_mos_mixed_material`: 659 assertions / 5 cases passed
- `test_linear_solver`: 19 assertions / 5 cases passed
- `test_analyze_sentaurus_bvmethods`: 3 tests passed
- 三个 Python 执行/分析脚本均通过 `py_compile`
