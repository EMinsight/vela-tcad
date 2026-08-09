# BV 仿真程序热点与计算瓶颈定位、改进计划

## 目标

在不调整迁移率、雪崩系数、带间隧穿参数或经验电流缩放的前提下，建立可重复的
Vela BV 性能测量体系，准确回答“时间花在哪里、为什么 Newton 次数多、为什么
单次更新慢”，并按证据依次降低完整 DD 评估次数、Newton 更新数和单次更新成本。

本计划只优化数值算法、线性代数、装配实现、状态预测和诊断开销；所有已有
Sentaurus/Vela 物理验收结果必须保持。

## 已知基线

| 指标 | Voltage-to-Current | 外接电阻 |
|---|---:|---:|
| 当前恢复段墙钟时间 | `1032.37 s` | `3885.10 s` |
| 完整 DD 边界评估 | 2 | 7 |
| Newton 更新 | 279 | 1356 |
| 有效时间/Newton | `3.700 s` | `2.865 s` |
| Vela BV | `6.395904175 V` | `6.395887866 V` |
| Sentaurus BV | `6.383184201 V` | `6.379791636 V` |

Sentaurus 完整运行为 `103.24 s`、370 次更新、约 `0.279 s/update`，但运行范围、
硬件和线程环境不同，只作为方向性参照。

## 总体执行顺序

1. 固定受控基准和正确性门槛。
2. 先用 `gprof` 获取函数级 flat profile 和 call graph，形成热点函数候选清单。
3. 再加入低开销、分层、可聚合的阶段计时，解释热点所处的求解上下文。
4. 用三个最小复现定位热点和迭代放大来源。
5. 先减少完整 DD 调用和 Newton 更新，再优化单次更新。
6. 每项优化独立提交、独立 A/B 测量，未通过正确性门槛立即回退。

## 阶段 0：建立受控基准

### 分支与构建

- 从执行时最新本地 `main` 创建 `codex/bv-performance-profiling`。
- 使用 MSYS2 UCRT64 Release 构建，固定编译器、优化级别、线程数和环境变量。
- 在基准报告中记录 CPU 型号、逻辑核数、内存、操作系统、编译器版本、CMake
  配置、Git 提交和命令行。
- 每个短基准至少运行 3 次，报告中位数和最大偏差；长时间完整验收可先运行 1 次，
  优化候选通过后再补 3 次。

### 三个基准场景

1. `6.08709 -> 6.09959 V` 高场最小复现：定位单个困难 Newton 轨迹。
2. Voltage-to-Current 最终 `1e-4 A/um`：定位两次完整 DD 评估的构成。
3. 外接电阻最终 `1206 V` 外部目标：定位 7 次完整 DD 评估和括区预测成本。

另保留一个完整 `0 -> 5.9 V -> boundary` 端到端场景，只用于最终验收，不用于
每次开发迭代。

### 基准输出

- `performance_run.json`：环境、配置哈希、Git 提交、总时间和计数。
- `performance_phases.csv`：每偏压、边界目标、DD 评估、Newton 迭代和阶段耗时。
- `performance_summary.md`：基线/候选比较、正确性门槛和结论。

## 阶段 1A：使用 gprof 定位函数级热点

本机已确认存在 `D:\msys64\ucrt64\bin\gprof.exe`，版本为 GNU Binutils
2.47；默认编译器为 UCRT64 GCC 16.1。第一轮热点定位使用 `gprof`，不先凭代码
阅读猜测热点。

### Profiling 构建

- 新增独立 `build-profile`，不能覆盖 Debug 或 Release 构建目录。
- 使用接近生产的优化级别并保留符号，建议 `RelWithDebInfo` 或
  `Release + -g`。
- 编译和链接必须同时加入 `-pg`；建议同时加入
  `-fno-omit-frame-pointer`，方便后续采样工具复核。
- 第一轮保持正常优化和内联，以反映真实优化构建；如果大量时间只能归入匿名或
  内联模板，再增加一个仅用于归因的 `-fno-inline-functions` 对照构建，不能把该
  对照构建的总时间当作生产性能。
- 推荐新增显式 CMake 选项 `VELA_ENABLE_GPROF=ON`，只对 Vela 目标增加 profiling
  编译/链接选项，默认关闭。

示例配置方式：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake -S . -B build-profile -G Ninja `
  -DCMAKE_BUILD_TYPE=RelWithDebInfo `
  -DVELA_ENABLE_GPROF=ON
cmake --build build-profile --parallel
```

### gprof 运行方式

每个基准使用独立工作目录，因为程序正常退出时会在当前目录生成 `gmon.out`。
不要强制终止被测进程，否则 profile 可能不会刷新。

```powershell
Set-Location runs/profile/high_field_transition
D:\code-repo\vela-tcad\build-profile\vela_example_runner.exe `
  --config simulation.json
D:\msys64\ucrt64\bin\gprof.exe -b -p -q `
  D:\code-repo\vela-tcad\build-profile\vela_example_runner.exe `
  gmon.out > gprof_report.txt
```

对三个基准分别保存：

- `gprof_flat.txt`：按 self time 排序，回答 CPU 时间直接消耗在哪些函数；
- `gprof_callgraph.txt`：按 cumulative time 和调用关系定位上游放大来源；
- `gmon.out`：原始 profile，连同可执行文件提交哈希和构建参数一起归档；
- `gprof_hotspots.csv`：提取函数、self time、cumulative time、calls、self/call。

第一轮候选门槛：self time 超过 `5%`、cumulative time 超过 `10%`，或调用次数
异常高的函数进入热点清单。重点观察但不预设结论：

- `CoupledDDAssembler::residual`、`assembleJacobian` 及雪崩/BTBT 子路径；
- `LinearSolver::solve`、Eigen SparseLU `factorize`/三角求解；
- `continuityRowWeights`、全局连续性闭合和稀疏行缩放；
- `BacktrackingLineSearch::search` 及候选残差计算；
- DCSweep 边界评估、状态预测、检查点和诊断输出。

### gprof 的边界

- `gprof` 提供函数级 CPU profile，但不能直接区分同一 `residual()` 是 Newton 初始
  残差、线搜索候选还是诊断调用。
- `-pg` 会插入函数调用计数，短函数和高频调用的耗时会受到扰动；因此只用它排序
  热点，不直接用 profile 构建的总时间作为最终加速比。
- Eigen 模板和被内联函数可能被归并到调用者；必要时用对照构建或外部采样 profiler
  复核。
- 当前 Vela/Eigen SparseLU 路径基本为单线程，适合 gprof；若后续引入多线程，需改用
  Windows Performance Recorder、Intel VTune 或 AMD uProf 等线程感知采样工具。
- I/O 等待和按偏压/边界目标划分的上下文仍由下一阶段内部计时补足。

## 阶段 1B：增加低开销分阶段计时

### 计时架构

新增默认关闭的 `solver.performance_profiling` 配置。使用单调时钟和 RAII scoped
timer，记录纳秒整数和调用次数；热循环内不写文件，只在偏压点结束或运行结束时
聚合输出。关闭时编译路径或运行分支的额外开销目标低于 `0.5%`，开启时低于 `2%`。

### 必须拆分的阶段

| 层级 | 计时项 | 主要代码位置 |
|---|---|---|
| DCSweep | 预偏置、边界目标、检查点读写、CSV/VTK 诊断 | `src/simulation/DCSweep.cpp` |
| BoundaryControl | 括区恢复、预测、DD 评估、根更新 | `src/simulation/BoundaryControl.cpp`、`DCSweep.cpp` |
| Newton | 初始残差、行权重、连续性闭合、Jacobian、行缩放、线性求解、线搜索、验收诊断 | `src/solver/NewtonSolver.cpp` |
| Assembler | 残差总计、Jacobian 总计、输运、复合、BTBT、雪崩源和雪崩 Jacobian | `src/equation/CoupledDDAssembler.cpp` |
| LinearSolver | pattern analyze、numeric factorize、solve | `src/solver/LinearSolver.cpp` |

### 同时记录的计数器

- Newton 更新数、残差求值数、Jacobian 装配数、线搜索尝试数。
- 稀疏矩阵维度、`nonZeros()`、符号分析次数、数值分解次数、solve 次数。
- 边界目标数、完整 DD 评估数、括区扩展数、割线预测命中/回退数。
- 每次装配的单元数、边数和活跃雪崩支持数。
- 诊断 CSV 行数、检查点字节数和 I/O 次数。

### 测试

- 关闭 profiling 时求解结果和现有 CSV 完全不变。
- 开启时所有阶段时间非负，父阶段时间不小于已计子阶段时间之和的合理范围。
- 调用计数与已知的 Newton/线搜索测试严格一致。
- 用伪时钟或容差断言测试聚合逻辑，不对真实墙钟时间做脆弱断言。

## 阶段 2：形成热点证据

对三个基准场景生成以下排名：

1. 总墙钟时间占比前 10 的阶段。
2. 单次 Newton 的平均/P50/P95 时间。
3. 每个边界 DD 评估的 Newton 数、初始残差、最终残差和状态预测来源。
4. 残差/Jacobian 次数是否因线搜索或诊断被放大。
5. SparseLU `factorize` 与 `solve` 的时间、矩阵规模和非零元变化。
6. 雪崩与 BTBT 装配占残差/Jacobian 的比例。

只有 gprof 与分阶段计时相互印证后，才决定单次迭代优化的内部优先级。当前“雪崩装配慢”或
“SparseLU 最慢”都只能视为候选假设。

## 阶段 3：优先减少完整 DD 评估

### 3.1 加强现有嵌套边界算法

- 复核跨目标括区端点恢复，确保正、负端点和对应状态均被复用。
- 使用最近两到三个已收敛边界状态构造受保护割线/切线预测。
- 预测同时作用于 InnerVoltage 和 `psi/phin/phip`，失败后回退常量 warm start。
- 根据局部 `dI/dV` 和负载线斜率调整括区步长，但不改变验收残差。

第一里程碑：

- 外接电阻每个新外部目标不超过 3 次完整 DD 评估；
- Voltage-to-Current 最终目标保持不超过 2 次完整 DD 评估；
- 所有负载线和电流边界残差不劣于当前基线。

### 3.2 联立边界约束

在计时和加强预测完成后，设计分块/Schur 边界未知量：

- 将接触 InnerVoltage 作为额外未知量；
- 外接电阻增加 `Vouter - Vinner - direction*I*R = 0`；
- 电流控制增加 `direction*I - Itarget = 0`；
- 在 Newton Jacobian 中加入接触电流对器件状态和接触电压的导数；
- 用 Schur 补或增广稀疏系统一次联立求解，消除外层嵌套标量根。

先保留原嵌套实现作为配置回退和结果对照。联立实现必须用单元测试验证符号、单位、
有限差分 Jacobian 和两种电流方向。

## 阶段 4：降低高场 Newton 更新数

按以下顺序逐项 A/B：

1. 保持连续 Newton 轨迹，避免同一偏压失败后从过旧检查点重启。
2. 使用电压/电流控制变量一致的二阶或受限割线状态预测。
3. 根据上一步收敛率自适应 InnerVoltage 步长和 continuation 步长。
4. 复核连续性行缩放在高场近零行上的条件数，优化数值尺度，不改变物理残差。
5. 将 carrier-row 保持为诊断，全局电子/空穴连续性继续作为硬验收。
6. 最后再评估伪弧长 continuation；不得以放宽残差或调物理参数代替收敛改进。

阶段性目标以实测为准，初始建议门槛：外接电阻 Newton 更新从 1356 至少下降
`50%`；Voltage-to-Current 从 279 至少下降 `25%`。

## 阶段 5：优化单次 Newton 成本

仅对阶段 2 证明的热点实施：

### 若装配占主导

- 缓存不随状态变化的几何、材料、边/单元拓扑和散射索引。
- 预分配 residual/Jacobian 工作区和 triplet 容量，减少热循环内分配。
- 将恒定 Jacobian 块与状态相关块分离，避免重复生成常量项。
- 合并残差和 Jacobian 中重复的迁移率、场、载流子统计和雪崩中间量计算。
- 确认诊断没有在生产求解路径重复重构同一物理量。

### 若稀疏数值分解占主导

- 保留当前符号 pattern 复用，并确认 pattern 是否实际稳定。
- 消除行缩放产生的不必要稀疏矩阵整份复制。
- 比较 Eigen SparseLU 排序策略和可用的稀疏直接求解后端。
- 评估带预条件的迭代求解只作为独立实验；必须保持残差和连续性门槛。
- 不错误复用已失效的数值分解；跨 Newton 复用必须有可证明的矩阵不变条件。

### 若线搜索占主导

- 复用已计算的候选残差和中间量。
- 使用上一迭代接受步长作为初值。
- 在不改变 Armijo/验收条件的前提下减少重复全量诊断。

当前最终接受状态绝大多数只进行一次线搜索尝试，因此线搜索暂列次要候选。

## 正确性硬门槛

每个优化提交必须同时满足：

- 外接电阻 BV 相对 `6.379791636 V` 误差小于 `3%`。
- Voltage-to-Current BV 相对 `6.383184201 V` 误差小于 `3%`。
- 与当前 Vela 基线的 BV 漂移建议不超过 `1 mV`；超过时必须解释数值路径变化。
- 电流边界残差不高于 `1e-10 A/um`，负载线残差不高于 `1e-4 V`。
- 全局电子/空穴连续性比率均不高于 `1e-2`。
- 电势 P95 误差、高场相关系数、载流子和准费米势比较不得显著退化。
- 配置中迁移率、雪崩、BTBT 参数及其哈希与基线一致。
- 原有 C++、模板和 BVmethods 回归测试全部通过。

## 性能验收门槛

分两级执行，避免先设不现实的最终目标：

### 第一里程碑

- profiling 关闭开销 `<0.5%`，开启开销 `<2%`；
- 外接电阻完整 DD 评估 `7 -> <=3`；
- 外接电阻恢复段墙钟时间至少降低 `40%`；
- Voltage-to-Current 恢复段墙钟时间至少降低 `20%`。

### 第二里程碑

- 根据热点数据，使单次 Newton 墙钟时间至少降低 `20%`；
- 完整端到端 BV 流程至少获得 `2x` 加速；
- 报告同机、同线程、同预偏置范围的 Sentaurus/Vela 对比。

若未达到门槛，报告必须区分“评估次数未降”“Newton 数未降”和“每更新成本未降”，
不得用总时间单一数字掩盖原因。

## 提交拆分建议

1. `perf(solver): add hierarchical BV profiling counters`
2. `test(perf): add reproducible BV hotspot benchmarks`
3. `perf(boundary): reduce repeated DD bracket evaluations`
4. `feat(boundary): add coupled contact constraint prototype`
5. `perf(assembler): optimize measured BV assembly hotspot`
6. `perf(linear): optimize measured sparse solve hotspot`
7. `docs(validation): report BV performance gains and correctness gates`

每个提交只解决一个可测问题，并附基线/候选 JSON；不要把 instrumentation、算法变化
和物理变化混在同一提交。

## 新聊天启动指令

```text
请在 D:\code-repo\vela-tcad 中执行 BV 仿真程序热点与计算瓶颈定位和改进任务。

先检查当前本地 main、工作区和最近提交，然后阅读：
docs/superpowers/plans/2026-08-08-bv-performance-hotspot-localization-and-improvement.md
docs/validation/bvmethods_nmos_field_iv_performance_analysis_2026-08-08.md
docs/validation/bvmethods_nmos_boundary_methods_implementation_2026-08-06.md

从最新本地 main 创建 codex/bv-performance-profiling 分支。严格按计划顺序执行：
先固定受控基准，再建立带 -pg 的独立 RelWithDebInfo profiling 构建。分别对
6.08709 -> 6.09959 V、Voltage-to-Current 最终目标和外接电阻 1206 V 三个最小
基准运行 gprof，保存 gmon.out、flat profile、call graph 和热点 CSV。根据 gprof
找出 self time、cumulative time 和调用次数最高的函数后，再加入低开销分阶段计时；
在获得 Vela 内部阶段数据前，不直接优化推测热点。对同样三个最小基准输出
每阶段时间、调用次数、Newton 数、残差/Jacobian 次数、
SparseLU analyze/factorize/solve 时间及边界 DD 评估次数。

依据数据，先减少完整 DD 评估和 Newton 更新，再优化单次更新成本。所有修改不得调整
迁移率、雪崩、BTBT 参数或经验电流缩放；外接电阻与 Voltage-to-Current 的 BV、
边界残差、全局连续性和场量比较必须持续通过既有门槛。每项优化独立提交，并保存
可机器比较的 baseline/candidate JSON。请先制定本次执行计划，然后直接开始检查和实现。
```
