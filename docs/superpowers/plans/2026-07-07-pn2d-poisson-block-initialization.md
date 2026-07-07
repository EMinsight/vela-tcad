# PN2D Poisson Block Initialization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an explicit, configurable Poisson block initialization pass before the first coupled Newton solve in PN2D BV sweeps, replacing ad-hoc external restart generation and avoiding the unusable `gummel_max_iter=1` handoff proxy.

**Architecture:** Reuse the existing Newton cold-start construction and `NewtonSolver::evaluateBlockStep(..., "poisson_only")` machinery inside production code. DCSweep will optionally build an in-memory `DDSolution` for the first bias point, feed it to the existing `initial` path, and record diagnostics showing the cold residual and Poisson-block residual. Default behavior remains unchanged; PN2D BV recommended configs explicitly enable the new mode.

**Tech Stack:** C++20, CMake/Ninja, Catch2 tests, nlohmann/json, existing Vela `DDSolution`, `NewtonSolver`, and `DCSweep` APIs.

## Global Constraints

- Windows development uses MSYS2 UCRT64: prepend `D:\msys64\ucrt64\bin;D:\msys64\usr\bin` to `PATH`.
- Keep generated simulation outputs under ignored `build*` or report directories; do not commit generated sweep CSV/VTK unless explicitly requested.
- Default off: existing regression and numerical results must be unchanged unless `sweep.initialization.mode` is set.
- `initial_state_file` remains supported; if both file restart and Poisson-block initialization are requested, fail loudly.
- The Poisson-block initializer is only a first-bias initializer, not a per-step continuation method.
- Use existing PN2D settings: `solver.handoff.gummel_max_iter=0`, `solver.handoff.require_gummel_convergence=false`, `solver.handoff.newton_max_iter=220`, `solver.warm_start=true`, `solver.quasi_fermi_update_limit_V=0.1`.

---

## File Structure

- Modify `include/vela/solver/NewtonSolver.h`
  - Add a public method that returns the Newton cold initial guess after one Poisson block step.
  - Keep the existing private assembler-aware `buildInitialGuess(...)` helper unchanged.

- Modify `src/solver/NewtonSolver.cpp`
  - Implement the public initializer by building the same assembler/BCs as `solve()`, constructing the cold initial guess, evaluating residuals, solving the Poisson block, and returning the trial solution.
  - Return a small diagnostic struct or expose enough values for DCSweep to report cold and trial block norms.

- Modify `include/vela/simulation/DCSweep.h`
  - Add `SweepInitializationConfig` and a member on `DCSweepConfig`.
  - Add optional point/summary fields for initialization mode and residual reduction.

- Modify `src/simulation/DCSweep.cpp`
  - Parse `sweep.initialization`.
  - Build the Poisson-block initial state before the first point solve.
  - Feed that state through the existing `initialState` pointer path.
  - Write diagnostics to summary CSV and optional diagnostic CSV.

- Modify `tests/test_newton_solver.cpp`
  - Add a focused unit test for the new NewtonSolver initializer.

- Modify `tests/test_dc_sweep.cpp`
  - Add config parsing and DCSweep integration tests.

- Modify PN2D config-generation scripts that produce the coarse7x3 recommended BV deck:
  - `scripts/run_pn2d_coarse7x3_previous_full20_compare.py`
  - `scripts/rebaseline_pn2d_post_contact_handoff_fix.py`
  - Any newer task-specific config generator that emits `fixed_sweep_strict_eps1e3` or avalanche A/B configs.

- Optional docs update:
  - `docs/validation/pn2d_bv_validation.md`, only if this branch normally records solver strategy changes there.

---

### Task 1: Add NewtonSolver Poisson Block Initializer

**Files:**
- Modify: `include/vela/solver/NewtonSolver.h`
- Modify: `src/solver/NewtonSolver.cpp`
- Test: `tests/test_newton_solver.cpp`

**Interfaces:**
- Consumes: existing private `NewtonSolver::buildInitialGuess(const CoupledDDAssembler&, const CoupledDDBoundaryConditions&) const`
- Consumes: existing `NewtonSolver::evaluateBlockStep(const DDSolution&, const std::string&) const`
- Produces:
  ```cpp
  struct NewtonPoissonBlockInitialization {
      DDSolution coldInitial;
      DDSolution poissonBlockInitial;
      NewtonBlockResidualInfo coldBlockResiduals;
      NewtonBlockResidualInfo poissonBlockResiduals;
      Real rawStepNorm = 0.0;
      Real stepNorm = 0.0;
  };

  NewtonPoissonBlockInitialization buildPoissonBlockInitialization() const;
  ```

- [ ] **Step 1: Declare the public result type and method**

In `include/vela/solver/NewtonSolver.h`, near the existing `NewtonBlockStepEvaluation` declaration, add:

```cpp
struct NewtonPoissonBlockInitialization {
    DDSolution coldInitial;
    DDSolution poissonBlockInitial;
    NewtonBlockResidualInfo coldBlockResiduals;
    NewtonBlockResidualInfo poissonBlockResiduals;
    Real rawStepNorm = 0.0;
    Real stepNorm = 0.0;
};
```

In the public `NewtonSolver` API, after `NewtonResult solve() const;`, add:

```cpp
NewtonPoissonBlockInitialization buildPoissonBlockInitialization() const;
```

- [ ] **Step 2: Write the failing unit test**

In `tests/test_newton_solver.cpp`, add a test next to existing Newton step/probe tests. Use the existing PN junction fixture helpers in that file; if no exact helper exists, reuse the smallest mesh/doping setup already used by Newton tests rather than creating a new fixture style.

The assertion shape should be:

```cpp
TEST_CASE("NewtonSolver builds a Poisson-block initialized cold state",
          "[newton][poisson_block_initialization]")
{
    auto fixture = makeCoupledPnJunctionFixture();
    NewtonConfig cfg;
    cfg.inputScaling.mode = "unit_scaling";
    cfg.warmStart = true;
    cfg.maxUpdate = 0.0;

    NewtonSolver solver(
        fixture.mesh,
        fixture.matdb,
        fixture.doping,
        fixture.contactBiases,
        cfg);

    const auto init = solver.buildPoissonBlockInitialization();
    const auto coldResidual = solver.evaluateResidual(init.coldInitial);
    const auto poissonResidual = solver.evaluateResidual(init.poissonBlockInitial);

    REQUIRE(init.coldInitial.psi.size() == init.poissonBlockInitial.psi.size());
    REQUIRE(init.rawStepNorm > 0.0);
    REQUIRE(poissonResidual.blockNorms.psi < coldResidual.blockNorms.psi);
    REQUIRE(poissonResidual.blockNorms.phin == Catch::Approx(coldResidual.blockNorms.phin));
    REQUIRE(poissonResidual.blockNorms.phip == Catch::Approx(coldResidual.blockNorms.phip));
}
```

If `makeCoupledPnJunctionFixture()` is not the local helper name, use the local fixture function already present in `tests/test_newton_solver.cpp` and keep the assertions unchanged.

- [ ] **Step 3: Run the focused test and confirm failure**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build-release --output-on-failure -R test_newton_solver
```

Expected: compile failure because `NewtonPoissonBlockInitialization` or `buildPoissonBlockInitialization()` does not exist.

- [ ] **Step 4: Implement the initializer**

In `src/solver/NewtonSolver.cpp`, implement the method immediately after `NewtonSolver::solve() const`:

```cpp
NewtonPoissonBlockInitialization NewtonSolver::buildPoissonBlockInitialization() const
{
    const double Vt = thermalVoltage(cfg_.temperature_K);
    const MobilityModelConfig mobilityConfig = cfg_.mobility;
    RecombinationModelConfig recombinationConfig =
        recombinationModelConfig(cfg_.recombination, cfg_.taun, cfg_.taup);
    recombinationConfig.augerCn = cfg_.augerCn;
    recombinationConfig.augerCp = cfg_.augerCp;
    const DDScalingSpec scaling = buildScalingSpec();
    CoupledDDAssembler assembler(
        mesh_,
        matdb_,
        doping_,
        Vt,
        mobilityConfig,
        recombinationConfig,
        cfg_.bandgapNarrowing,
        cfg_.impactIonization,
        fixedCharges_,
        sheetCharges_,
        scaling,
        cfg_.carrierDiagonalFloor);
    const CoupledDDBoundaryConditions bcs = buildBoundaryConditions(assembler);

    NewtonPoissonBlockInitialization out;
    out.coldInitial = buildInitialGuess(assembler, bcs);
    const NewtonBlockStepEvaluation step = evaluateBlockStep(out.coldInitial, "poisson_only");
    out.poissonBlockInitial = step.trialSolution;
    out.coldBlockResiduals = step.residual.blockNorms;
    out.poissonBlockResiduals = step.trialResidual.blockNorms;
    out.rawStepNorm = step.rawStepNorm;
    out.stepNorm = step.stepNorm;
    return out;
}
```

- [ ] **Step 5: Run the focused test and confirm pass**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure -R test_newton_solver
```

Expected: `test_newton_solver` passes.

- [ ] **Step 6: Commit**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add include/vela/solver/NewtonSolver.h src/solver/NewtonSolver.cpp tests/test_newton_solver.cpp
git commit -m "Add Newton Poisson block initializer"
```

---

### Task 2: Add DCSweep JSON Configuration and Parsing

**Files:**
- Modify: `include/vela/simulation/DCSweep.h`
- Modify: `src/simulation/DCSweep.cpp`
- Test: `tests/test_dc_sweep.cpp`

**Interfaces:**
- Consumes: `DCSweepConfig`
- Produces:
  ```cpp
  struct SweepInitializationConfig {
      std::string mode = "none";
      std::string diagnosticCsv;
      std::string writeStateFile;
  };
  ```

Supported JSON:

```json
"sweep": {
  "initialization": {
    "mode": "poisson_block",
    "diagnostic_csv": "poisson_block_initialization.csv",
    "write_state_file": "poisson_block_initial_state.csv"
  }
}
```

- [ ] **Step 1: Add config struct**

In `include/vela/simulation/DCSweep.h`, before `struct DCSweepConfig`, add:

```cpp
struct SweepInitializationConfig {
    std::string mode = "none";
    std::string diagnosticCsv;
    std::string writeStateFile;
};
```

Inside `DCSweepConfig`, add:

```cpp
SweepInitializationConfig initialization;
```

- [ ] **Step 2: Write parsing tests**

In `tests/test_dc_sweep.cpp`, add a test that writes a minimal sweep config with:

```json
"initialization": {
  "mode": "poisson_block",
  "diagnostic_csv": "init.csv",
  "write_state_file": "init_state.csv"
}
```

Then run `DCSweepSimulation` on a tiny existing PN fixture and assert:

```cpp
REQUIRE(std::filesystem::exists(dir / "init.csv"));
REQUIRE(std::filesystem::exists(dir / "init_state.csv"));
```

Add a second test where both `"initial_state_file": "restart.csv"` and `"initialization": {"mode": "poisson_block"}` are present. Assert that `runWithResult()` throws and the message contains:

```text
DCSweep: sweep.initialization.mode='poisson_block' cannot be combined with initial_state_file
```

- [ ] **Step 3: Implement parsing**

In `src/simulation/DCSweep.cpp`, inside `dcSweepConfigFromJson(...)` after `sweep.initialStateFile` is parsed, add:

```cpp
if (j.contains("initialization")) {
    const auto& init = j.at("initialization");
    if (!init.is_object())
        throw std::invalid_argument("DCSweep: sweep.initialization must be an object.");
    sweep.initialization.mode = init.value("mode", std::string("none"));
    if (sweep.initialization.mode != "none" &&
        sweep.initialization.mode != "poisson_block") {
        throw std::invalid_argument(
            "DCSweep: sweep.initialization.mode must be 'none' or 'poisson_block'.");
    }
    sweep.initialization.diagnosticCsv =
        init.value("diagnostic_csv", std::string{});
    sweep.initialization.writeStateFile =
        init.value("write_state_file", std::string{});
}
if (!sweep.initialStateFile.empty() &&
    sweep.initialization.mode == "poisson_block") {
    throw std::invalid_argument(
        "DCSweep: sweep.initialization.mode='poisson_block' cannot be combined with initial_state_file.");
}
```

In the config path resolution block near existing `sweep.initialStateFile = resolve(...)`, add:

```cpp
if (!sweep.initialization.diagnosticCsv.empty())
    sweep.initialization.diagnosticCsv = resolve(sweep.initialization.diagnosticCsv);
if (!sweep.initialization.writeStateFile.empty())
    sweep.initialization.writeStateFile = resolve(sweep.initialization.writeStateFile);
```

- [ ] **Step 4: Run focused tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure -R test_dc_sweep
```

Expected: new parsing tests pass.

- [ ] **Step 5: Commit**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add include/vela/simulation/DCSweep.h src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp
git commit -m "Parse DCSweep Poisson block initialization"
```

---

### Task 3: Wire Poisson Block Initialization Into DCSweep Execution

**Files:**
- Modify: `src/simulation/DCSweep.cpp`
- Test: `tests/test_dc_sweep.cpp`

**Interfaces:**
- Consumes: `NewtonSolver::buildPoissonBlockInitialization() const`
- Produces:
  - First sweep point receives the Poisson-block state through the existing `initial` pointer path.
  - Diagnostic CSV columns:
    ```csv
    bias_V,cold_psi,cold_phin,cold_phip,poisson_psi,poisson_phin,poisson_phip,raw_step_norm,step_norm
    ```

- [ ] **Step 1: Add diagnostic writer helper**

In `src/simulation/DCSweep.cpp`, near other CSV writer helpers, add:

```cpp
void writePoissonBlockInitializationCsv(
    const std::filesystem::path& path,
    Real bias,
    const NewtonPoissonBlockInitialization& init)
{
    if (path.empty())
        return;
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out.is_open())
        throw std::runtime_error(
            "DCSweep: cannot write Poisson block initialization diagnostics: " + path.string());
    out << "bias_V,cold_psi,cold_phin,cold_phip,"
        << "poisson_psi,poisson_phin,poisson_phip,"
        << "raw_step_norm,step_norm\n";
    out << std::setprecision(17)
        << bias << ','
        << init.coldBlockResiduals.psi << ','
        << init.coldBlockResiduals.phin << ','
        << init.coldBlockResiduals.phip << ','
        << init.poissonBlockResiduals.psi << ','
        << init.poissonBlockResiduals.phin << ','
        << init.poissonBlockResiduals.phip << ','
        << init.rawStepNorm << ','
        << init.stepNorm << '\n';
}
```

- [ ] **Step 2: Build the initialization state before the bias loop**

In `runDCSweep(...)`, replace the existing initial-state setup:

```cpp
std::unique_ptr<DDSolution> initialState;
if (!sweep.initialStateFile.empty()) {
    initialState = std::make_unique<DDSolution>(
        readDDSolutionStateCsv(sweep.initialStateFile, mesh.numNodes()));
}
```

with:

```cpp
std::unique_ptr<DDSolution> initialState;
if (!sweep.initialStateFile.empty()) {
    initialState = std::make_unique<DDSolution>(
        readDDSolutionStateCsv(sweep.initialStateFile, mesh.numNodes()));
} else if (sweep.initialization.mode == "poisson_block") {
    Real initializationBias = sweep.start;
    if (!sweep.biasPoints.empty())
        initializationBias = sweep.biasPoints.front();
    auto initializationBiases = baseBiases;
    initializationBiases[sweep.contact] = initializationBias;

    NewtonConfig initializationNewton = newton;
    initializationNewton.warmStart = true;
    NewtonSolver initializationSolver(
        mesh,
        matdb,
        doping,
        initializationBiases,
        initializationNewton,
        fixedChargeSpecs,
        sheetChargeSpecs);
    const NewtonPoissonBlockInitialization init =
        initializationSolver.buildPoissonBlockInitialization();
    initialState = std::make_unique<DDSolution>(init.poissonBlockInitial);

    writePoissonBlockInitializationCsv(
        std::filesystem::path(sweep.initialization.diagnosticCsv),
        initializationBias,
        init);
    if (!sweep.initialization.writeStateFile.empty())
        writeDDSolutionStateCsv(sweep.initialization.writeStateFile, *initialState);
}
```

- [ ] **Step 3: Mark first point handoff stage**

In the first-point solve path, after `attempt = solvePointWithContinuation(...)` and before `savePoint(...)`, append the initialization marker only for the first point:

```cpp
if (sweep.initialization.mode == "poisson_block" && points.empty() &&
    attempt.handoffStage == "newton") {
    attempt.handoffStage = "poisson_block_newton";
}
```

Do not alter continuation points; they should keep normal `newton` or `newton_failed` stages.

- [ ] **Step 4: Extend the DCSweep integration test**

In the test from Task 2, read `init.csv` and assert:

```cpp
const auto rows = readCsvRows(dir / "init.csv");
REQUIRE(rows.size() == 1);
REQUIRE(std::stod(rows[0].at("cold_psi")) > std::stod(rows[0].at("poisson_psi")));
REQUIRE(std::stod(rows[0].at("raw_step_norm")) > 0.0);
```

Then read the sweep output CSV and assert the first row has:

```cpp
REQUIRE(rows[0].at("handoff_stage") == "poisson_block_newton");
```

- [ ] **Step 5: Run focused tests**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure -R "test_dc_sweep|test_newton_solver"
```

Expected: both focused test binaries pass.

- [ ] **Step 6: Commit**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add src/simulation/DCSweep.cpp tests/test_dc_sweep.cpp
git commit -m "Use Poisson block initialization in DCSweep"
```

---

### Task 4: Update PN2D BV Recommended Config Generation

**Files:**
- Modify: `scripts/run_pn2d_coarse7x3_previous_full20_compare.py`
- Modify: `scripts/rebaseline_pn2d_post_contact_handoff_fix.py`
- Modify: any current PN2D BV config generator that writes `fixed_sweep_strict_eps1e3` or `avalanche_gain_ab` configs

**Interfaces:**
- Produces JSON fragment:
  ```json
  "sweep": {
    "initialization": {
      "mode": "poisson_block",
      "diagnostic_csv": "<report_dir>/poisson_block_initialization.csv",
      "write_state_file": "<report_dir>/poisson_block_initial_state.csv"
    }
  }
  ```

- [ ] **Step 1: Locate config generation helpers**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
rg -n "gummel_max_iter|fixed_sweep_strict|avalanche_gain|carrier_row_convergence|output_csv|write_state_file" scripts
```

Expected: identify every script that writes the current PN2D BV coarse7x3 solver/sweep JSON.

- [ ] **Step 2: Add explicit Poisson-block initialization to each recommended config**

Where a script constructs `config["sweep"]`, add:

```python
config["sweep"]["initialization"] = {
    "mode": "poisson_block",
    "diagnostic_csv": str(report_dir / "poisson_block_initialization.csv"),
    "write_state_file": str(report_dir / "poisson_block_initial_state.csv"),
}
```

Keep:

```python
config["solver"]["handoff"]["gummel_max_iter"] = 0
config["solver"]["handoff"]["require_gummel_convergence"] = False
```

Do not set `gummel_max_iter=1`; the controlled experiment showed that path fails at 0V with `linear_solve_failed`.

- [ ] **Step 3: Add or update regression tests for script output**

In `tests/regression/test_reference_tcad_tools.py`, add assertions for each generated recommended PN2D BV config:

```python
self.assertEqual(config["sweep"]["initialization"]["mode"], "poisson_block")
self.assertTrue(config["sweep"]["initialization"]["diagnostic_csv"].endswith("poisson_block_initialization.csv"))
self.assertTrue(config["sweep"]["initialization"]["write_state_file"].endswith("poisson_block_initial_state.csv"))
self.assertEqual(config["solver"]["handoff"]["gummel_max_iter"], 0)
```

- [ ] **Step 4: Run regression tests for config generation**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
python -m unittest tests.regression.test_reference_tcad_tools -v
```

Expected: regression tests pass.

- [ ] **Step 5: Commit**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add scripts tests/regression/test_reference_tcad_tools.py
git commit -m "Enable Poisson block initialization in PN2D BV configs"
```

---

### Task 5: Validation on coarse7x3 BV Sweep

**Files:**
- No source edits expected
- Generated reports under:
  `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/poisson_block_initialization_validation_20260707/`

**Interfaces:**
- Consumes: updated PN2D BV JSON config with `sweep.initialization.mode="poisson_block"`
- Produces:
  - baseline vs Poisson-block table
  - `poisson_block_initialization.csv`
  - sweep summary CSV

- [ ] **Step 1: Generate two controlled configs**

Use the same base as the current strict coarse7x3 BV config:

- Baseline: `sweep.initialization.mode="none"`
- Candidate: `sweep.initialization.mode="poisson_block"`

Both must retain:

```json
"solver": {
  "method": "gummel_newton",
  "handoff": {
    "fallback": "none",
    "require_gummel_convergence": false,
    "gummel_max_iter": 0,
    "newton_max_iter": 220
  },
  "warm_start": true,
  "quasi_fermi_update_limit_V": 0.1
}
```

- [ ] **Step 2: Run both sweeps**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\vela_example_runner.exe --config <baseline_config.json>
.\build-release\vela_example_runner.exe --config <poisson_block_config.json>
```

Expected:

- Baseline cell-reconstructed strict config remains 55/55.
- Poisson-block config remains 55/55.
- `poisson_block_initialization.csv` reports `poisson_psi < cold_psi`.

- [ ] **Step 3: Report initial residual and early iteration comparison**

Extract:

- `0V` first Newton residual from `newton_history`.
- Newton iterations at biases `0`, `-0.002`, `-0.004386666667`, `-0.007250666667`, `-0.01067792`.
- Current at `-15.2432189285`, `-18`, `-20`.

Expected based on the controlled experiment:

| metric | baseline | poisson_block |
|---|---:|---:|
| 0V first Newton residual | about `2.04e-2` | about `2.42e-4` |
| 0V Newton iterations | 3 | 2 |
| -20V current, cell_reconstructed | about `-5.16e-17 A/um` | about `-5.16e-17 A/um` |

- [ ] **Step 4: Re-run grad_qf failure control**

Run baseline and Poisson-block versions for `current_approximation="grad_qf"`.

Expected based on the controlled experiment:

| case | last converged bias | failed bias | failure |
|---|---:|---:|---|
| grad_qf baseline | `-15.2432189285` | `-15.5691757933` | `carrier_row_convergence_line_search_rejected` |
| grad_qf poisson_block | `-15.2432189285` | `-15.5691757933` | `carrier_row_convergence_line_search_rejected` |

If the failure point changes, preserve the new diagnostics and do not call the experiment failed; treat it as useful evidence that initialization now interacts with the high-field source branch.

- [ ] **Step 5: Run full test suite**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build-release --output-on-failure
```

Expected: all tests pass.

- [ ] **Step 6: Commit validation docs if requested**

Only commit generated validation summaries if the user asks to preserve them. Otherwise leave generated CSV/report outputs uncommitted.

---

### Task 6: Final Documentation and Safety Review

**Files:**
- Modify: `docs/validation/pn2d_bv_validation.md` if this branch is already tracking PN2D solver strategy there
- No generated CSV committed unless requested

**Interfaces:**
- Produces a short documented solver strategy:
  - `gummel_max_iter=0` remains the recommended handoff.
  - `sweep.initialization.mode="poisson_block"` replaces manual restart/probe CSV generation.
  - `gummel_max_iter=1` is explicitly not recommended for this case because it fails at 0V in the controlled experiment.

- [ ] **Step 1: Add a validation note**

Append a section:

```markdown
## PN2D Explicit Poisson Block Initialization

The coarse7x3 BV deck now uses `sweep.initialization.mode="poisson_block"` for the first bias point.
This performs one Newton Poisson block solve from the same cold initial state used by the coupled Newton solver, then hands the trial state to the existing coupled Newton path.

Controlled experiment summary:
- Cold initial 0V combined block residual: `2.704409e+1`.
- Poisson-block initial 0V combined block residual: `5.510267e-1`.
- First recorded 0V Newton residual improved from about `2.04e-2` to about `2.42e-4`.
- cell_reconstructed strict BV remains 55/55 to `-20V`.
- grad_qf still fails at `-15.5691757933V`, so that blocker is not caused by 0V initialization.
```

- [ ] **Step 2: Check default-off numerical neutrality**

Run one known baseline config without `sweep.initialization` and compare summary CSV against the pre-change run for:

- converged point count
- first five `newton_iterations`
- `current_total_A_per_um` at `0`, `-15.2432189285`, `-18`, `-20`

Expected: unchanged within existing floating-point output stability.

- [ ] **Step 3: Commit docs**

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git add docs/validation/pn2d_bv_validation.md
git commit -m "Document PN2D Poisson block initialization"
```

Skip this commit if no docs update is made.

---

## Self-Review

- Spec coverage: The plan covers code, config, tests, validation, and default-off safety.
- Placeholder scan: No task uses TBD/TODO; each implementation step names exact files, API names, JSON keys, commands, and expected outcomes.
- Type consistency: `NewtonPoissonBlockInitialization`, `SweepInitializationConfig`, and `sweep.initialization.mode` are defined before use.
- Risk note: The public Newton initializer duplicates assembler setup from `solve()`. This is acceptable for the first implementation because it preserves exact solver semantics; a later refactor can share an internal `makeCoupledAssembler()` helper after tests are green.
