# PN2D Minimal6 Operator Audit 交接文档（2026-07-13）

## 1. 本次暂停原因与当前结论

用户要求暂停当前开发，回家后在另一台电脑继续。本次暂停发生在：

- Task 5 最终独立复审已完成，结论为双 `FAIL`。
- Task 5 复审修复已完成 RED 阶段并进入 GREEN 实现中段。
- 修复代理已被明确停止。
- Task 6 真实六状态 C++ 重放尚未开始。
- 暂停时没有 Python、CMake、Ninja 或 minimal6 审计进程在运行。

最重要的状态判断：当前工作区是一个可保留、但尚未验证的 Task 5 WIP
checkpoint。不要把它当作已完成实现，也不要直接从 Task 6 开始。

## 2. 仓库与分支状态

工作区：

```text
D:\code-repo\vela-tcad
```

分支：

```text
codex-pn2d-minimal6-operator-audit
```

暂停前已提交的最新提交：

```text
9546826 Fix PN2D minimal6 fixed-state audit provenance
26d66b4 Add PN2D minimal6 fixed-state reports
37a9545 Add Vela fixed-state operator audit
e33fba7 Fix PN2D minimal6 Sentaurus log recovery
4c6ebee Tighten PN2D minimal6 state bias validation
719becd Add PN2D minimal6 exact state exports
a7e9603 Add PN2D minimal6 Sentaurus topology gate
```

暂停时相对远端同名分支为 `ahead 10`。本交接只创建本地 checkpoint，
不会自动推送；换电脑前必须自行推送此分支或以其他方式传输 checkpoint。

本次 checkpoint 包含以下 WIP 文件：

```text
scripts/audit_pn2d_minimal6_fixed_state.py
tests/regression/test_pn2d_minimal6_fixed_state_audit.py
docs/validation/pn2d_minimal6_handoff_2026-07-13.md
```

## 3. Tasks 1-4 与真实状态基线

Tasks 1-4 已实现、提交并通过各自的实现后复审流程：

- Task 1：显式 sketch/mirror 六节点拓扑和 DF-ISE 输入。
- Task 2：Sentaurus 无重网格拓扑门禁。
- Task 3：两拓扑在 `0/-12/-19 V` 的精确状态导出。
- Task 4：生产 C++ 固定态 operator audit，可输出 node/edge/triangle CSV，
  不运行 Newton、Gummel 或 continuation。

真实 Task 3 状态根位于：

```text
build-release/reference_tcad/pn2d_sentaurus2018_minimal6/
  state_exports/minimal6_states_live_20260713_v2
```

暂停前已确认：

- manifest schema 为 `vela.pn2d_minimal6_states.v1`。
- `outputs_complete=true`。
- sketch/mirror × `0/-12/-19 V` 共六个状态全部为 `passed`。
- requested bias 与 actual bias 精确相等。
- 每个 export 已有 `nodes.csv`、`elements.csv`、`contacts.csv`、
  `doping.csv`、`state.csv`、`field_manifest.json` 和原始 fields。
- 真实 export 尚无 `mesh.json`、`audit.json`、`vela_node_state.csv`、
  `vela_edge_audit.csv`、`vela_triangle_audit.csv`。

Task 4 独立复审报告：

```text
.superpowers/sdd/task-4-minimal6-report.md
```

## 4. Task 5 在最终复审前的证据

Task 5 初版提交为 `26d66b4`，第一轮修复提交为 `9546826`。修复后曾得到：

```text
Focused Python tests: 25/25 OK
Task 1 + Task 3 Python regression grouping: 36 tests OK
Task 4 C++ test_fixed_state_operator_audit: 63 assertions, 3 cases OK
Synthetic CLI: PASS, node=36 edge=54 triangle=24 figures=14
Max state parity hybrid error: 0
Max C++/Python formula hybrid error: 1.7898855730891228e-13
```

这些结果只能证明 `9546826` 当时的既有测试通过，不能覆盖最终复审发现的
门禁绕过。Task 5 最终独立复审范围为：

```text
37a9545..9546826
```

## 5. Task 5 最终独立复审结论

最终复审明确给出：

```text
SPEC COMPLIANCE: FAIL
CODE/DATA QUALITY: FAIL
Task 6 may not start.
```

### 5.1 Critical：Task 4 provenance 记录了但未强制执行

旧的主 CLI 只检查 manifest 中是否有六个 `exit_code=0`，没有在生成报告前
强制调用 replay verifier。对抗性复现已经证明：把 `producer_sha256` 和输入
哈希改成伪造值后，CLI 仍返回 0、打印 PASS，并写出 `gate_status=PASS`。

修复要求：任何报告文件落盘前必须验证：

- producer 名称和 producer SHA-256。
- Task 4 source commit。
- 精确的六个 topology/bias replay identity，且无缺失、重复或额外记录。
- 精确参数和命令记录。
- 六个 recorded exit code 均为 0。
- 每个 replay 的全部输入哈希。
- 当前已提交的 18 个 C++ 输出哈希。
- 使用目标 executable 重新运行后的 18 个输出哈希。

### 5.2 Critical：Vela source 汇总列由 Python 值填充

旧实现把 Python 重建的局部源项累加后写入
`vela_*_source_integral_per_m_s`，图和外部诊断也使用这些值。这违反
“Python 不能生成 Vela production columns”的要求。

修复要求：

- `vela_*` 汇总必须只累加 Task 4 CSV 的 raw local-edge source。
- `python_*` 使用独立公式重建并单独汇总。
- C++ 与 Python 汇总之间必须有严格 formula gate。
- 图和 Sentaurus/Vela diagnostic 必须使用真正的 C++ `vela_*` 汇总。

### 5.3 Important：几何零归一化范围过宽

旧的 `1e-285` 小源项归一化只判断两个源项都很小，没有确认对应 partial
volume 是否为数学零。因此非零体积上的明显相对误差也能被抹成 0。

修复要求：只有 C++ 与 Python partial volume 都落在几何零阈值内时，才允许
对相应源项做 tiny-zero 归一化；非零体积必须使用普通 hybrid gate。

### 5.4 Important：source oracle 循环复用 C++ mobility/alpha

旧 Python source oracle 直接读取 Task 4 CSV 中的 mobility 和 alpha 来构造
expected source。内部一致但错误的 C++ 因子因此可能逃过检查。

修复要求：从不可变 `audit.json`、Silicon 默认材料参数和独立
Van-Overstraeten 公式重建 mobility/alpha，或先对这些输入做独立 gate，再用于
期望值。任何 expected formula 都不能由 Vela-labeled 数值直接填充。

### 5.5 Important：拓扑与原始掺杂仍有绕过路径

- inline topology contract 旧实现只比较排序后的 node set，反向 tuple 可通过。
- 原始 Task 3 donor/acceptor fields 没有与 canonical doping / `doping.csv`
  逐节点比较。

修复要求：验证每个 triangle 是批准的 CCW tuple，并逐节点验证原始 donor 和
acceptor 场的单位换算及语义。

### 5.6 Important：对抗性测试和图件 QA 不足

必须加入：

- 缺失、重复、伪造 provenance。
- producer/input/committed-output/fresh-output hash 篡改。
- Vela/Python aggregate 分离。
- 非零体积 tiny source 拒绝。
- mobility/alpha 篡改。
- raw donor/acceptor 篡改。
- 反向 inline topology tuple。
- 精确 14 个图件文件名、PNG/PDF signature、可解码性、非空白像素和
  fixed-state/not-BV 文案。

## 6. 当前 WIP 的 RED 与中断点

修复代理按 TDD 先写测试。生产代码修改前的 RED 结果为：

```text
Ran 32 tests in 25.282s
FAILED (8 failures, 3 errors)
```

RED 覆盖了预期缺口：

- 坏 producer hash 未阻止报告生成。
- provenance 缺失、重复、伪造记录未失败。
- committed output tamper 未失败。
- fresh replay hash tamper 缺少可区分错误。
- Vela/Python aggregate 不存在或混用。
- `geometric_source_gate` 不存在。
- mobility/alpha corruption 只通过循环的 generic formula 偶然失败。
- raw Task 3 donor corruption 通过。
- 图件增强测试最初还暴露过 PIL import 漏写，WIP 中已补上 import。

RED 后已开始修改生产代码，但暂停前没有运行 GREEN。因此当前 32 个测试的
实际通过/失败状态未知。

## 7. 当前 WIP 已写入但尚未验证的内容

当前 `scripts/audit_pn2d_minimal6_fixed_state.py` 大致新增/修改了：

- 固定的 Task 4 producer/source-commit 常量和哈希验证。
- `geometric_source_gate()`，尝试把 tiny-source 归一化限制到几何零。
- 独立 `van_overstraeten_alpha()`。
- 固定 Silicon constant mobility 的独立 gate。
- `load_audit_model()`，要求每个 export 有 immutable `mesh.json` 和
  `audit.json`。
- raw Task 3 donor/acceptor 对 `doping.csv` 的逐节点 gate。
- C++ `vela_*` 与 Python `python_*` 三角源项的分离汇总。
- C++/Python aggregate source formula gate。
- 更严格的 `verify_task4_replay()`：六状态 identity、参数、命令、输入哈希、
  committed outputs、fresh outputs。
- `write_report()` 在创建输出目录前调用 provenance verification。

当前测试文件已增加 7 个左右的对抗性测试组，并加强 14 图件 QA。

这些改动是 WIP，不能直接视为正确实现。

## 8. 恢复后首先检查的高风险点

### 8.1 不要错误要求 `elements.csv` 的行顺序等于 canonical tuple 顺序

当前 WIP 的 `validate_topology()` 按 element `id` 排序后，直接要求 tuple list
与 `TRIS[topology]` 完全同序。这很可能过严。

真实 sketch export 的 element 行顺序为：

```text
0: (1,5,2)
1: (2,4,3)
2: (2,6,4)
3: (5,6,2)
```

它包含四个批准的 CCW triangle，但行顺序不同于 canonical list。恢复时应验证：

- 每个 tuple 本身精确为批准的 CCW tuple，不能反向。
- 集合完整且无重复。
- canonical report key/order 由明确的 canonicalization 产生。
- 不应仅因 Sentaurus/export 行顺序不同而拒绝真实数据。

先修正这一点，再跑 Task 5 GREEN，否则真实 Task 6 state root 很可能被误拒绝。

### 8.2 独立 alpha 必须逐项对照生产 C++

检查：

```text
src/physics/ImpactIonizationModel.cpp
include/vela/physics/ImpactIonizationModel.h
include/vela/equation/AssemblerUtils.h
```

确认 A/B 参数的 SI 单位、switch field、temperature gamma、minimum field 和
exponent 处理与 production 完全一致。不要只凭 fixture 当前可通过来确认公式。

### 8.3 provenance 失败必须在任何 artifact 创建之前发生

`write_report()` 应先验证 provenance，失败时 `out-dir` 不应存在。CLI 的 PASS
必须只在 `write_report()` 完成且 summary/gate 为 PASS 后打印。

### 8.4 真实 Task 6 必须先补 replay provenance

真实 Task 3 manifest 原本没有 Task 4 provenance。Task 6 逐状态 C++ 重放时要在
忽略目录中生成/记录：

- topology-specific `mesh.json`。
- 与 committed fixture 相同的 `audit.json`。
- 三类 Vela CSV。
- 六条完整 replay command、exit、input/output hashes 和 producer hash。

否则新的严格 Task 5 CLI 应按设计拒绝真实 report。

## 9. 恢复执行顺序

### Step A：获取 checkpoint 并确认没有丢改动

在家里电脑执行：

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Set-Location D:\code-repo\vela-tcad
git fetch
git switch codex-pn2d-minimal6-operator-audit
git status --short --branch
git log -15 --oneline
```

如果本次 checkpoint 尚未推送，先在公司电脑推送该分支，或用 `git bundle` /
其他安全方式传输。不要在家里从 `9546826` 重新写一遍而忽略 checkpoint。

### Step B：只读审查 WIP diff

```powershell
git show --stat --oneline HEAD
git show --check HEAD
git diff HEAD^ -- `
  scripts/audit_pn2d_minimal6_fixed_state.py `
  tests/regression/test_pn2d_minimal6_fixed_state_audit.py
```

优先修复第 8.1 节的 element 行顺序问题，并逐项核对第 5 节所有复审发现。

### Step C：完成 Task 5 GREEN

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_fixed_state_audit -v
```

然后运行 synthetic production CLI：

```powershell
D:\msys64\ucrt64\bin\python.exe `
  scripts\audit_pn2d_minimal6_fixed_state.py `
  --state-root tests\fixtures\pn2d_minimal6_synthetic `
  --out-dir build-release\reference_tcad\pn2d_minimal6_synthetic_audit_resume
```

要求：

- 32 个聚焦测试全部通过（如果继续增加测试，数量可以更高）。
- CLI provenance replay 实际执行，不能跳过。
- 输出精确 `36/54/24`。
- 14 个图件 QA 全部通过。
- 最大 C++/Python formula error `<5e-12`。
- 最大 state parity error `<1e-12`。

### Step D：运行 Tasks 1、3、4 回归

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_pn2d_minimal6_sentaurus_gate `
  tests.regression.test_pn2d_minimal6_state_export `
  tests.regression.test_pn2d_minimal6_fixed_state_audit -v

build-release\test_fixed_state_operator_audit.exe
git diff --check
```

### Step E：提交 Task 5 修复并重新独立复审

只提交 Task 5 所属文件，更新：

```text
.superpowers/sdd/task-5-report.md
```

报告必须记录 RED、GREEN、CLI、replay 和全部复审问题关闭证据。建议提交信息：

```text
Enforce PN2D minimal6 audit provenance
```

重新生成完整 Task 5 review package，范围为：

```text
37a9545..新的 Task 5 HEAD
```

必须重新取得两个明确结论：

```text
SPEC COMPLIANCE: PASS
CODE/DATA QUALITY: PASS
```

在此之前禁止启动 Task 6。

### Step F：双 PASS 后执行 Task 6

Task 6 brief：

```text
.superpowers/sdd/task-6-minimal6-brief.md
```

真实状态根：

```text
build-release/reference_tcad/pn2d_sentaurus2018_minimal6/
  state_exports/minimal6_states_live_20260713_v2
```

顺序：

1. 为六个真实 export 生成 topology-matched `mesh.json` 并复制相同 committed
   `audit.json`。
2. 对六个状态分别执行 `build-release/pn2d_minimal6_operator_audit.exe`。
3. 要求六个命令成功，并记录完整 provenance/hashes。
4. 用 `--state-root` 生成真实 joined audit，要求 `36/54/24`。
5. 检查 7 PNG + 7 PDF，包括人工视觉检查。
6. 更新：
   - `docs/validation/pn2d_bv_validation.md`
   - `docs/validation/pn2d_bv_current_progress_summary.md`
7. 明确声明：没有 Vela nonlinear solve，也没有物理 BV curve。
8. 运行 Task 6 完整聚焦验证、提交并独立复审。

Task 6 计划提交信息：

```text
Document PN2D minimal6 operator audit
```

## 10. Task 6 完整验证命令

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"

cmake --build build-release --target `
  test_fixed_state_operator_audit `
  test_impact_ionization `
  test_cell_reconstructed_avalanche `
  pn2d_minimal6_operator_audit --parallel 2

build-release\test_fixed_state_operator_audit.exe
build-release\test_impact_ionization.exe
build-release\test_cell_reconstructed_avalanche.exe

D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_pn2d_minimal6_sentaurus_gate `
  tests.regression.test_pn2d_minimal6_state_export `
  tests.regression.test_pn2d_minimal6_fixed_state_audit -v

git diff --check
```

## 11. 暂停状态总表

```text
Task 1: committed, reviewed
Task 2: committed, live topology gate passed
Task 3: committed, six exact live states exported
Task 4: committed, reviewed
Task 5 implementation at 9546826: tests passed but final review FAIL/FAIL
Task 5 review-fix RED: captured (32 tests, 8 failures, 3 errors)
Task 5 review-fix GREEN implementation: WIP, interrupted before verification
Task 5 fresh re-review: pending
Task 6 real six-state replay: not started
Task 6 validation docs: not updated
```

恢复后的最高优先级不是运行真实六状态，而是把当前 Task 5 WIP 做到 GREEN，
修正 element-order 过严风险，并取得 Task 5 的独立双 PASS。
