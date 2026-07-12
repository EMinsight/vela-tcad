# PN2D Minimal6 Fixed-State Operator Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an exact six-node, four-triangle Sentaurus/Vela fixed-state operator audit for the sketch and vertically mirrored PN2D topologies at 0 V, -12 V, and -19 V.

**Architecture:** A canonical JSON topology generates explicit DF-ISE `.grd/.dat` files that Sentaurus must accept without remeshing. Sentaurus supplies immutable nodal states; a dedicated Vela C++ audit API evaluates production SG, gradient, current-proxy, alpha, geometry, and source operators without solving, while Python independently recomputes the documented mathematics and joins the results.

**Tech Stack:** Python 3 standard library, C++20, Eigen, nlohmann-json, Catch2, CMake/Ninja, Sentaurus Device 2018 in the SSH VM, DF-ISE grid/dataset text, existing Vela TDR import and CSV utilities, matplotlib/Pillow for static figures.

## Global Constraints

- The fixture is diagnostic-only and must not contribute a physical BV curve or replace an existing PN2D reference.
- Geometry is exactly `2.0 um x 0.5 um`; the junction is exactly `x=1.0 um`.
- Nodes 2 and 6 have `Donor=Acceptor=1e17 cm^-3`; no legacy p-side rewrite is allowed.
- Required topologies are `sketch` and `mirror`; each has 6 nodes, 4 counter-clockwise triangles, 9 unique edges, 6 boundary edges, 3 interior edges, and 2 contact edges.
- Required exact Sentaurus biases are `0`, `-12`, and `-19 V`; nearest-bias substitution is forbidden.
- Vela must use the Sentaurus TDR and imported state without Newton, Gummel, continuation, or carrier/QF modification.
- Coordinate interpolation, nearest-neighbor matching, automatic remeshing fallback, QF clamps, source limiters, alpha scaling, contact fallback, and solver tuning are forbidden.
- Production defaults such as `density_gradient` and existing triangle-GSS behavior remain unchanged.
- Manual edits use ASCII and generated simulation outputs remain under ignored `build-release/` directories.

---

### Task 1: Canonical Topologies And DF-ISE Writer

**Files:**
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/source/minimal6_topologies.json`
- Create: `scripts/pn2d_minimal6_topology.py`
- Create: `tests/regression/test_pn2d_minimal6_topology.py`

**Interfaces:**
- Consumes: no earlier task output.
- Produces: `load_topology(path: Path, topology_id: str) -> Topology`; `validate_topology(topology: Topology) -> TopologySummary`; `write_dfise_grid(topology: Topology, path: Path) -> None`; `write_dfise_doping(topology: Topology, path: Path) -> None`; `validate_dfise_roundtrip(topology: Topology, grd: Path, dat: Path) -> dict[str, object]`.

- [ ] **Step 1: Write failing topology and round-trip tests**

```python
def test_sketch_and_mirror_have_exact_contract(self):
    sketch = module.load_topology(FIXTURE, "sketch")
    mirror = module.load_topology(FIXTURE, "mirror")
    self.assertEqual(sketch.triangles, [(1, 5, 2), (5, 6, 2), (2, 6, 4), (2, 4, 3)])
    self.assertEqual(mirror.triangles, [(1, 5, 6), (1, 6, 2), (2, 6, 3), (6, 4, 3)])
    for topology in (sketch, mirror):
        summary = module.validate_topology(topology)
        self.assertEqual((summary.nodes, summary.triangles, summary.edges), (6, 4, 9))
        self.assertEqual(summary.contact_edges, {"Anode": (1, 5), "Cathode": (3, 4)})

def test_dfise_roundtrip_preserves_connectivity_contacts_and_doping(self):
    topology = module.load_topology(FIXTURE, "sketch")
    with tempfile.TemporaryDirectory() as tmp:
        grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
        module.write_dfise_grid(topology, grd)
        module.write_dfise_doping(topology, dat)
        report = module.validate_dfise_roundtrip(topology, grd, dat)
    self.assertTrue(report["passed"])
    self.assertEqual(report["silicon_triangle_count"], 4)
    self.assertEqual(report["contact_element_counts"], {"Anode": 1, "Cathode": 1})
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_topology -v
```

Expected: import failure because `scripts/pn2d_minimal6_topology.py` does not exist.

- [ ] **Step 3: Add the canonical topology JSON**

The JSON must contain these exact canonical values:

```json
{
  "schema": "vela.pn2d_minimal6_topologies.v1",
  "length_unit": "um",
  "doping_unit": "cm^-3",
  "nodes": {
    "1": [0.0, 0.5], "2": [1.0, 0.5], "3": [2.0, 0.5],
    "4": [2.0, 0.0], "5": [0.0, 0.0], "6": [1.0, 0.0]
  },
  "contacts": {"Anode": [1, 5], "Cathode": [3, 4]},
  "acceptors_cm3": {"1": 1e17, "2": 1e17, "3": 0.0, "4": 0.0, "5": 1e17, "6": 1e17},
  "donors_cm3": {"1": 0.0, "2": 1e17, "3": 1e17, "4": 1e17, "5": 0.0, "6": 1e17},
  "topologies": {
    "sketch": [[1,5,2], [5,6,2], [2,6,4], [2,4,3]],
    "mirror": [[1,5,6], [1,6,2], [2,6,3], [6,4,3]]
  }
}
```

- [ ] **Step 4: Implement immutable topology types and validation**

```python
@dataclass(frozen=True)
class Topology:
    topology_id: str
    nodes: dict[int, tuple[float, float]]
    triangles: list[tuple[int, int, int]]
    contacts: dict[str, tuple[int, int]]
    acceptors_cm3: dict[int, float]
    donors_cm3: dict[int, float]

def canonical_edges(triangles):
    return sorted({tuple(sorted(edge)) for tri in triangles
                   for edge in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0]))})
```

`validate_topology` must reject non-CCW triangles, any count mismatch, any coordinate other than the six approved coordinates, a changed contact, and non-compensated nodes 2 or 6. Mirror validation must reflect `y -> 0.5-y` and map labels `1<->5`, `2<->6`, `3<->4`.

- [ ] **Step 5: Implement deterministic DF-ISE `.grd` and `.dat` output**

Enumerate sorted canonical edges once. Encode a forward edge as its zero-based ID and a reversed edge as `-(id+1)`. The `.grd` must declare 6 vertices, 9 edges, 6 elements, and regions `R.Si`, `Cathode`, `Anode`; elements 0-3 are triangles, element 4 is contact segment 3-4, and element 5 is contact segment 1-5. Emit `Locations` as `i` for the three interior edges and `e` for the six boundary edges.

The `.dat` must emit six vertex values for each dataset:

```text
DopingConcentration = donor - acceptor
PhosphorusActiveConcentration = donor
BoronActiveConcentration = acceptor
```

Use `scripts.compare_sentaurus_tdr_tdx.parse_grd` and `parse_dat` in `validate_dfise_roundtrip`; reconstruct triangle vertex connectivity from signed edge IDs before comparison.

- [ ] **Step 6: Run focused and existing DF-ISE tests**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_sentaurus_import_tools -v
```

Expected: all tests pass; the new module reports two valid topology IDs.

- [ ] **Step 7: Commit Task 1**

```powershell
D:\msys64\ucrt64\bin\git.exe add reference_tcad/pn2d_sentaurus2018_minimal6/source/minimal6_topologies.json scripts/pn2d_minimal6_topology.py tests/regression/test_pn2d_minimal6_topology.py
D:\msys64\ucrt64\bin\git.exe commit -m "Add PN2D minimal6 explicit topology fixtures"
```

---

### Task 2: Sentaurus Explicit-Grid Compatibility Gate

**Files:**
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/pn2d_sentaurus2018_minimal6_reference.json`
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par`
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_gate_sdevice.cmd`
- Create: `scripts/run_pn2d_minimal6_sentaurus_gate.py`
- Create: `tests/regression/test_pn2d_minimal6_sentaurus_gate.py`

**Interfaces:**
- Consumes: all Task 1 topology interfaces and generated `.grd/.dat` files.
- Produces: `build_gate_bundle(topology_id: str, output_dir: Path) -> GateBundle`; `validate_returned_tdr(topology: Topology, neutral_export: Path) -> dict[str, object]`; CLI manifest schema `vela.pn2d_minimal6_sentaurus_gate.v1`.

- [ ] **Step 1: Write failing dry-run, command-array, and returned-topology tests**

```python
def test_dry_run_builds_two_explicit_grid_commands_without_sde(self):
    manifest = gate.prepare_gate(...)
    self.assertEqual([run["topology_id"] for run in manifest["runs"]], ["sketch", "mirror"])
    for run in manifest["runs"]:
        self.assertNotIn("sde", " ".join(run["remote_commands"]))
        self.assertIn("sdevice pn2d_minimal6_gate_sdevice.cmd", run["remote_commands"][0])

def test_returned_tdr_gate_rejects_added_node(self):
    with self.assertRaisesRegex(ValueError, "expected 6 nodes"):
        gate.validate_returned_tdr(SKETCH, export_with_seven_nodes)
```

- [ ] **Step 2: Run the tests and verify RED**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_sentaurus_gate -v
```

Expected: module import failure.

- [ ] **Step 3: Add the reference descriptor and gate deck**

The reference descriptor names `minimal6_topologies.json`, the generated grid and doping files, `models.par`, the gate deck, both topology IDs, and the coordinate tolerance `1e-12 um`. Copy the existing coarse7x3 `models.par` byte-for-byte and record its SHA-256 in the generated manifest.

The gate deck uses the staged names `pn2d_minimal6.grd` and `pn2d_minimal6.dat`:

```text
File {
  Grid="pn2d_minimal6.grd"
  Doping="pn2d_minimal6.dat"
  Parameter="models.par"
  Plot="pn2d_minimal6_gate_des.tdr"
  Output="pn2d_minimal6_gate_des.log"
}
Electrode {
  { Name="Anode" Voltage=0.0 }
  { Name="Cathode" Voltage=0.0 }
}
Physics { Mobility(DopingDep HighFieldSaturation) Recombination(SRH Auger Avalanche) EffectiveIntrinsicDensity(OldSlotboom) }
Solve { Coupled(Iterations=100 LineSearchDamping=1e-4) { Poisson Electron Hole } }
```

- [ ] **Step 4: Implement dry-run and live VM orchestration**

Reuse `default_windows_openssh`, `run_checked`, and `write_manifest` from `scripts/run_sentaurus_vm_reference.py`. Build argv arrays; do not construct local shell strings. For each topology, stage exactly the generated `.grd/.dat`, deck, and `models.par`, then run SDevice and copy back `.tdr`, `.log`, `.grd`, `.dat`, and stdout.

CLI:

```text
python scripts/run_pn2d_minimal6_sentaurus_gate.py
  --topologies sketch,mirror
  --ssh-target sentaurus
  --run-id <stable-id>
  --output-dir <ignored-build-path>
  [--dry-run]
```

- [ ] **Step 5: Implement the strict returned-TDR gate**

Run `build-release/sentaurus_import.exe export-neutral` on each returned TDR. Canonically map coordinates with absolute error `<1e-12 um`, reconstruct the triangle sets, contacts, and nodal donor/acceptor fields, and compare all values against Task 1. Emit a per-topology `topology_gate.json`; any mismatch returns nonzero.

- [ ] **Step 6: Run dry-run tests and generate a local gate bundle**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_sentaurus_gate -v
D:\msys64\ucrt64\bin\python.exe scripts\run_pn2d_minimal6_sentaurus_gate.py --dry-run --run-id minimal6_gate_prepare_20260712
```

Expected: tests pass; the manifest lists two SDevice-only runs and no SDE command.

- [ ] **Step 7: Execute the live compatibility gate and enforce the stop condition**

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\run_pn2d_minimal6_sentaurus_gate.py --run-id minimal6_gate_live_20260712
```

Expected: both topology gates report exactly `6/4/9`, exact contact edges, exact compensated doping, and unchanged connectivity. If either topology fails, stop the implementation here, preserve logs and manifests, document the incompatibility, and do not start Task 3.

- [ ] **Step 8: Commit Task 2 after the live gate passes**

```powershell
D:\msys64\ucrt64\bin\git.exe add reference_tcad/pn2d_sentaurus2018_minimal6 scripts/run_pn2d_minimal6_sentaurus_gate.py tests/regression/test_pn2d_minimal6_sentaurus_gate.py
D:\msys64\ucrt64\bin\git.exe commit -m "Add PN2D minimal6 Sentaurus topology gate"
```

---

### Task 3: Exact-Bias Sentaurus State Export Contract

**Files:**
- Create: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd`
- Create: `scripts/export_pn2d_minimal6_states.py`
- Create: `tests/regression/test_pn2d_minimal6_state_export.py`

**Interfaces:**
- Consumes: Task 2 live-compatible grid bundles.
- Produces: six export directories keyed by `(topology_id, bias_V)`; `field_manifest.json`; `state.csv` in `DDSolutionCsv` columns; manifest schema `vela.pn2d_minimal6_states.v1`.

- [ ] **Step 1: Write failing exact-bias and field-contract tests**

```python
def test_required_state_matrix_is_exact_and_complete(self):
    matrix = export.validate_state_matrix(exports)
    self.assertEqual(set(matrix), {(t, v) for t in ("sketch", "mirror") for v in (0.0, -12.0, -19.0)})

def test_vector_fields_require_two_components_and_cm_units(self):
    with self.assertRaisesRegex(ValueError, "eCurrentDensity.*A\\*cm\\^-2"):
        export.validate_field_manifest(manifest_with_wrong_current_unit)
```

- [ ] **Step 2: Run the tests and verify RED**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_state_export -v
```

- [ ] **Step 3: Add an exact-bias state deck**

Use the gate-compatible explicit grid and the same physical model configuration as the coarse7x3 BV deck. Generate one deck invocation per exact target bias so the final plot filename encodes the requested value. Ramp from 0 V with bounded steps but accept only a final contact voltage matching the request within `1e-12 V`. Plot every field required by the design, including vector ElectricField and vector carrier CurrentDensity.

- [ ] **Step 4: Implement state export and manifest validation**

Use `sentaurus_import.exe export-neutral` for each final TDR. Require vector field entries with `region=0`, `components=2`, `mapping_status=complete`, and `global_node_mapping=global_vertex_order`. Require current units `A*cm^-2`, mobility units `cm^2*V^-1*s^-1`, density units `cm^-3`, alpha units `cm^-1`, and QF/potential units `V`.

Create `state.csv` with canonical node order and SI values:

```text
node_id,psi_V,phin_V,phip_V,n_m3,p_m3
```

Do not derive missing `phin/phip` from density; absence is fatal.

- [ ] **Step 5: Run tests and a mocked six-state export**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_state_export -v
```

Expected: exact six-state fixture passes; nearest-bias, scalar-vector, unit, and missing-QF fixtures fail with explicit messages.

- [ ] **Step 6: Run all six real Sentaurus states**

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\export_pn2d_minimal6_states.py --run-id minimal6_states_20260712 --biases=0,-12,-19
```

Expected: six exact converged exports. If a bias fails, retain a partial manifest marked `outputs_complete=false`; do not substitute another bias and do not continue to final report generation.

- [ ] **Step 7: Commit Task 3**

```powershell
D:\msys64\ucrt64\bin\git.exe add reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd scripts/export_pn2d_minimal6_states.py tests/regression/test_pn2d_minimal6_state_export.py
D:\msys64\ucrt64\bin\git.exe commit -m "Add PN2D minimal6 exact state exports"
```

---

### Task 4: Vela C++ Fixed-State Operator Evaluator

**Files:**
- Create: `include/vela/equation/FixedStateOperatorAudit.h`
- Create: `src/equation/FixedStateOperatorAudit.cpp`
- Create: `src/tools/pn2d_minimal6_operator_audit.cpp`
- Create: `tests/test_fixed_state_operator_audit.cpp`
- Modify: `CMakeLists.txt`

**Interfaces:**
- Consumes: Task 3 mesh JSON, doping CSV, `state.csv`, and an audit JSON containing the existing impact-ionization and mobility configuration.
- Produces: `FixedStateOperatorAuditResult evaluateFixedStateOperators(const DeviceMesh&, const VectorXd& doping, const DDSolution&, const GummelConfig&)`; CLI CSV files `vela_node_state.csv`, `vela_edge_audit.csv`, and `vela_triangle_audit.csv`.

- [ ] **Step 1: Write failing C++ tests for one triangle and the full sketch topology**

```cpp
TEST_CASE("fixed-state audit preserves supplied carrier state", "[minimal6][fixed-state]") {
    const auto result = evaluateFixedStateOperators(mesh, doping, state, config);
    REQUIRE(result.nodes.size() == 6);
    REQUIRE(result.nodes[1].psi == Catch::Approx(state.psi[1]).margin(0.0));
    REQUIRE(result.nodes[1].n == Catch::Approx(state.n[1]).margin(0.0));
}

TEST_CASE("fixed-state audit enumerates nine edges and four triangles", "[minimal6][fixed-state]") {
    const auto result = evaluateFixedStateOperators(mesh, doping, state, config);
    REQUIRE(result.edges.size() == 9);
    REQUIRE(result.triangles.size() == 4);
    REQUIRE(result.triangles[0].signedDoubleArea > 0.0);
}
```

- [ ] **Step 2: Add the test target and verify RED**

```powershell
D:\msys64\ucrt64\bin\cmake.exe --build build-release --target test_fixed_state_operator_audit --parallel 2
```

Expected: compile failure because `FixedStateOperatorAudit.h` does not exist.

- [ ] **Step 3: Define explicit result records**

```cpp
struct FixedStateNodeRecord { Index nodeId; Real psi, phin, phip, n, p; };
struct FixedStateEdgeRecord {
    Index edgeId, node0, node1;
    Real length, electronRawSignedFlux, holeRawSignedFlux;
    Real electronMidpointDensity, holeMidpointDensity;
    Real electronImpactField, holeImpactField;
    Real electronAlpha, holeAlpha, edgeArea;
};
struct FixedStateTriangleRecord {
    Index cellId;
    std::array<Index, 3> nodes;
    Real signedDoubleArea;
    Point2 gradPsi, gradPhin, gradPhip;
    std::vector<TriangleGssAvalancheSourceRecord> localEdges;
};
```

The public names above are stable contracts for Task 5. Add only fields required by the design; do not expose solver mutation APIs.

- [ ] **Step 4: Implement production-helper evaluation without solving**

Build edge/cell/node adjacency using `AssemblerUtils.h`. Call production SG quasi-Fermi flux helpers, `cellReconstructedAvalancheMidpointDensity`, `computeCellScalarGradientCache`, `edgeQuasiFermiCoefficientField`, `selectAvalancheCurrentFluxProxy`, `avalancheSourceEdgeArea`, and `triangleGssAvalancheSourceRecords`. Do not duplicate these formulas in the new C++ source.

Copy the supplied `DDSolution` into node records before evaluation and assert it is byte-for-byte unchanged afterward. Reject non-Tri3 cells, non-six-node meshes, non-four-cell meshes, non-finite state, and any state vector with the wrong size.

- [ ] **Step 5: Add the CLI and deterministic CSV serialization**

CLI:

```text
pn2d_minimal6_operator_audit
  --mesh mesh.json
  --doping doping.csv
  --state state.csv
  --config audit.json
  --node-out vela_node_state.csv
  --edge-out vela_edge_audit.csv
  --triangle-out vela_triangle_audit.csv
```

Use `CsvUtils` and canonical node-pair keys. Include native signed values and SI units in column names. The CLI performs no solve and returns nonzero if counts are not exactly 6/9/4.

- [ ] **Step 6: Run focused and avalanche regression tests**

```powershell
D:\msys64\ucrt64\bin\cmake.exe --build build-release --target test_fixed_state_operator_audit test_impact_ionization test_cell_reconstructed_avalanche pn2d_minimal6_operator_audit --parallel 2
build-release\test_fixed_state_operator_audit.exe
build-release\test_impact_ionization.exe
build-release\test_cell_reconstructed_avalanche.exe
```

Expected: all tests pass; existing avalanche assertion counts do not regress.

- [ ] **Step 7: Commit Task 4**

```powershell
D:\msys64\ucrt64\bin\git.exe add CMakeLists.txt include/vela/equation/FixedStateOperatorAudit.h src/equation/FixedStateOperatorAudit.cpp src/tools/pn2d_minimal6_operator_audit.cpp tests/test_fixed_state_operator_audit.cpp
D:\msys64\ucrt64\bin\git.exe commit -m "Add Vela fixed-state operator audit"
```

---

### Task 5: Independent Mathematical Reference And Joined Reports

**Files:**
- Create: `scripts/audit_pn2d_minimal6_fixed_state.py`
- Create: `tests/regression/test_pn2d_minimal6_fixed_state_audit.py`
- Create: `tests/fixtures/pn2d_minimal6_synthetic/`

**Interfaces:**
- Consumes: Task 1 canonical topology, Task 3 Sentaurus exports, and Task 4 Vela CSV files.
- Produces: exact 36-row `node_state.csv`, 54-row `edge_audit.csv`, 24-row `triangle_audit.csv`, `summary.json`, `summary.md`, and seven PNG/PDF figure pairs; schema `vela.pn2d_minimal6_fixed_state_audit.v1`.

- [ ] **Step 1: Write failing Python-reference and matrix tests**

```python
def test_linear_triangle_gradient_matches_closed_form(self):
    grad = audit.triangle_gradient([(0,0), (1,0), (0,1)], [2, 5, 7])
    self.assertEqual(grad, (3.0, 5.0))

def test_complete_matrix_has_exact_row_counts_and_unique_keys(self):
    report = audit.build_report(two_topology_six_state_fixture)
    self.assertEqual(len(report.node_rows), 36)
    self.assertEqual(len(report.edge_rows), 54)
    self.assertEqual(len(report.triangle_rows), 24)
    audit.require_unique_keys(report.edge_rows, ("topology_id", "bias_V", "node0", "node1"))
```

- [ ] **Step 2: Run tests and verify RED**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_fixed_state_audit -v
```

- [ ] **Step 3: Implement the independent reference formulas**

Implement triangle shape gradients from the inverse `[[1,x,y], ...]` matrix, stable Bernoulli limits, production sign conventions for electron and hole SG flux, GSS logistic midpoint density, canonical vector projection, and the documented Genius-truncated partial volume. These functions may import constants but must not import Vela-produced CSV values as their result.

Use the hybrid error norm:

```python
def hybrid_error(actual, expected, abs_floor=1e-300):
    return abs(actual - expected) / max(abs(actual), abs(expected), abs_floor)
```

Classify two zero values as `both_zero`; do not assign a log10 ratio.

- [ ] **Step 4: Implement strict joins and gates**

Canonicalize every node by exact coordinate, every edge by sorted canonical node pair, and every triangle by the approved CCW tuple. Require field-manifest units before converting. Fail on missing or duplicate keys, inexact bias, partial state, wrong topology, non-finite nonzero values, C++/Python formula error `>=5e-12`, state parity error `>=1e-12`, or row-count mismatch.

Sentaurus versus Vela current and source errors are diagnostic columns only. Do not add a required improvement threshold.

- [ ] **Step 5: Implement orientation sensitivity and static figures**

For each bias and canonical mirrored quantity, compute `mirror/sketch`, signed difference, absolute log10 ratio when both values are nonzero, and zero classification. Generate:

```text
minimal6-topologies.{png,pdf}
minimal6-edge-current-audit-{0v,minus12v,minus19v}.{png,pdf}
minimal6-triangle-source-audit-{0v,minus12v,minus19v}.{png,pdf}
```

The topology figure labels all nodes, edges, triangles, contacts, and compensated nodes. Current/source figures use scientific axes and state explicitly that they are fixed-state operator audits, not BV curves.

- [ ] **Step 6: Run Python tests and a synthetic full report**

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_fixed_state_audit -v
D:\msys64\ucrt64\bin\python.exe scripts\audit_pn2d_minimal6_fixed_state.py --fixture tests/fixtures/pn2d_minimal6_synthetic --out-dir build-release/reference_tcad/pn2d_minimal6_synthetic_audit
```

Expected: all gates pass; output row counts are exactly 36/54/24 and every PNG/PDF is nonempty.

- [ ] **Step 7: Commit Task 5**

```powershell
D:\msys64\ucrt64\bin\git.exe add scripts/audit_pn2d_minimal6_fixed_state.py tests/regression/test_pn2d_minimal6_fixed_state_audit.py tests/fixtures/pn2d_minimal6_synthetic
D:\msys64\ucrt64\bin\git.exe commit -m "Add PN2D minimal6 fixed-state reports"
```

---

### Task 6: Real Six-State Audit, Documentation, And Final Verification

**Files:**
- Modify: `docs/validation/pn2d_bv_validation.md`
- Modify: `docs/validation/pn2d_bv_current_progress_summary.md`
- Modify: `docs/superpowers/specs/2026-07-12-pn2d-minimal6-operator-audit-design.md` only if the live compatibility gate reveals a documented Sentaurus file-format constraint without changing approved behavior.

**Interfaces:**
- Consumes: all outputs and gates from Tasks 1-5.
- Produces: ignored real artifact root `build-release/reference_tcad/pn2d_sentaurus2018_minimal6/reports/minimal6_fixed_state_audit_<date>` and durable validation conclusions.

- [ ] **Step 1: Run the real C++ replay for all six states**

For each exact export, run `pn2d_minimal6_operator_audit` with the same committed audit config and write topology/bias-specific Vela CSV files. Require six successful command records before joining.

- [ ] **Step 2: Generate the real joined audit**

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\audit_pn2d_minimal6_fixed_state.py `
  --state-root build-release/reference_tcad/pn2d_sentaurus2018_minimal6/states/minimal6_states_20260712 `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018_minimal6/reports/minimal6_fixed_state_audit_20260712
```

Expected: topology/state/formula gates pass; 36/54/24 rows are generated; Sentaurus/Vela differences and sketch/mirror sensitivity are reported without a physical BV claim.

- [ ] **Step 3: Inspect all static figures**

Check PNG dimensions and nonblank pixel range, then visually inspect topology labels, legends, scientific units, zero classifications, and lack of overlap. Re-render after any presentation-only correction.

- [ ] **Step 4: Update validation documents**

Document the exact grid route, all gate statuses, formula-parity errors, per-bias current/source discrepancies, and topology-direction ratios. State prominently that no Vela solve and no physical BV curve were performed. If the live grid gate failed, document that failure instead and do not claim Tasks 3-5 completed on real Sentaurus data.

- [ ] **Step 5: Run complete focused verification**

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

D:\msys64\ucrt64\bin\git.exe diff --check
```

Expected: every focused test passes and `git diff --check` is silent.

- [ ] **Step 6: Commit Task 6**

```powershell
D:\msys64\ucrt64\bin\git.exe add docs/validation/pn2d_bv_validation.md docs/validation/pn2d_bv_current_progress_summary.md
D:\msys64\ucrt64\bin\git.exe commit -m "Document PN2D minimal6 operator audit"
```

## Plan Self-Review

- Every approved spec requirement maps to a task: topology and doping in Task 1, Sentaurus no-remesh gate in Task 2, exact biases and fields in Task 3, production C++ replay in Task 4, independent reference and 36/54/24 reports in Task 5, real artifacts and documentation in Task 6.
- All interfaces use the same stable keys: `topology_id`, exact `bias_V`, canonical node IDs, sorted canonical edge pairs, and approved CCW triangle tuples.
- Task 2 is a hard execution checkpoint. A live Sentaurus topology failure blocks later implementation and is reported rather than bypassed.
- No task changes a production physical model, runs a Vela nonlinear solve, publishes a minimal6 BV curve, or commits generated simulation artifacts.
