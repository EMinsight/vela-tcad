# PN2D Minimal6 Operator Audit Handoff (2026-07-12)

## 1. Objective And Scope

This branch implements a diagnostic-only PN2D comparison on an exact six-node,
four-triangle mesh. The purpose is to compare Sentaurus and Vela operators on
the same immutable state and topology. It must not be presented as a physical
BV curve and must not change production avalanche defaults.

Approved geometry and state matrix:

- Rectangle: `2.0 um x 0.5 um`.
- Junction: `x = 1.0 um`.
- Nodes: `(0,0.5)`, `(1,0.5)`, `(2,0.5)`, `(2,0)`, `(0,0)`, `(1,0)`.
- Contacts: `Anode=(1,5)`, `Cathode=(3,4)`.
- Nodes 2 and 6: donor and acceptor both `1e17 cm^-3`.
- Topologies: `sketch` and vertically reflected `mirror`.
- Exact Sentaurus biases: `0 V`, `-12 V`, and `-19 V`.
- Vela receives Sentaurus state and evaluates operators without Newton,
  Gummel, continuation, interpolation, carrier/QF modification, clamps, or
  source scaling.

## 2. Repository State At Pause

Workspace:

```text
D:\code-repo\vela-tcad
```

Branch:

```text
codex-pn2d-minimal6-operator-audit
```

Relevant commits, oldest first:

```text
e2c97ee Design PN2D minimal6 operator audit
943cb48 Plan PN2D minimal6 operator audit
3af7268 Add PN2D minimal6 explicit topology fixtures
```

The first two commits define the approved behavior. Commit `3af7268` implements
Task 1. No Task 2 Sentaurus gate code or live Sentaurus result exists yet.

Primary documents:

- `docs/superpowers/specs/2026-07-12-pn2d-minimal6-operator-audit-design.md`
- `docs/superpowers/plans/2026-07-12-pn2d-minimal6-operator-audit.md`
- `.superpowers/sdd/task-1-brief.md`
- `.superpowers/sdd/task-1-report.md`

Task 1 files:

- `reference_tcad/pn2d_sentaurus2018_minimal6/source/minimal6_topologies.json`
- `scripts/pn2d_minimal6_topology.py`
- `tests/regression/test_pn2d_minimal6_topology.py`

## 3. Completed Work And Evidence

Task 1 currently provides:

- Exact canonical JSON for the sketch and mirror topologies.
- Validation of node coordinates, CCW triangles, edge counts, contacts,
  compensated junction nodes, and reflected connectivity.
- Deterministic DF-ISE `.grd` and `.dat` writers modeled on the local accepted
  PN2D files.
- Signed edge encoding and triangle reconstruction.
- Roundtrip checks for vertices, nine edges, six exterior and three interior
  edge locations, four silicon triangles, region ownership, contact endpoints,
  and the three doping datasets.
- Robust dynamic module import on Python 3.14 and repo-root CLI execution.

Controller verification performed immediately before the pause:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_sentaurus_import_tools -v
```

Result:

```text
Ran 47 tests in 14.651s
OK
```

The Task 1 implementer also reported:

```text
Focused minimal6 topology suite: 7 tests, OK
git diff HEAD --check: exit code 0, no output
```

## 4. Unresolved Review Findings

Task 1 is committed and its tests pass, but independent review was interrupted
by the pause request and returned important findings. Treat Task 1 as requiring
a fix-and-review cycle before starting the live Sentaurus gate.

### 4.1 Topology containers are not deeply immutable

`Topology` uses `@dataclass(frozen=True)`, but its dictionaries and triangle
list can still be mutated. The approved interface called for immutable topology
types. Convert nested collections to immutable representations or otherwise
prevent mutation, while preserving practical read access and deterministic
ordering.

Location at pause:

```text
scripts/pn2d_minimal6_topology.py:41
```

### 4.2 Roundtrip validation relies on permissive parsers

`validate_dfise_roundtrip()` uses the repository regex parsers. Those parsers
do not prove all required DF-ISE metadata, including version/type/dimension,
coordinate-system completeness, dataset function/type/location/validity, and
some declared counts. A malformed file may therefore still report
`passed=True`.

Location at pause:

```text
scripts/pn2d_minimal6_topology.py:398
```

Add explicit metadata validation in the minimal6 validator. Keep the existing
shared parsers for structural extraction, but do not treat their successful
parse as sufficient consumer-format validation.

### 4.3 Negative tests are too narrow

The tests corrupt a contact endpoint, but do not yet prove rejection of:

- Invalid signed-edge references or disconnected triangle edge loops.
- Incorrect `Locations` exterior/interior classification.
- Incorrect `R.Si` element ownership.
- Missing or renamed doping datasets.
- Incorrect dataset counts or values.
- Missing `function`, `location=vertex`, or `validity=["R.Si"]` metadata.
- Incorrect grid version/type/dimension or declared mesh counts.

Add focused failure tests for these cases before calling Task 1 review-clean.

Important distinction: the emitted files visually follow local references and
all current tests pass, but actual Sentaurus acceptance has not been tested.

## 5. Required Resume Sequence

Do not begin Task 2 live execution immediately. Resume in this order.

### Step A: Restore and inspect the exact branch

```powershell
Set-Location D:\code-repo\vela-tcad
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
git status --short --branch
git log -5 --oneline
```

Confirm that `3af7268` is present and inspect any newer handoff commit. Do not
discard unrelated changes. If continuing on another computer, fetch or transfer
this branch before doing any implementation.

### Step B: Fix and re-review Task 1

Use test-driven changes for the three review findings. Keep the fix separate
from Task 2. Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology -v

D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_topology `
  tests.regression.test_sentaurus_import_tools -v

git diff --check
```

Request a fresh independent review against `943cb48..HEAD`. Required verdicts:

- Spec compliance: `PASS`.
- Code quality: no Critical or Important findings.

Suggested separate commit message:

```text
Harden PN2D minimal6 DF-ISE validation
```

### Step C: Implement Task 2 dry-run gate

Follow Task 2 in the implementation plan exactly. Required new files include
the reference descriptor, copied `models.par`, SDevice gate deck, orchestration
script, and regression tests.

The orchestration must:

- Generate both explicit topologies from Task 1.
- Stage `.grd/.dat`, deck, and `models.par` only.
- Run SDevice without SDE and without remeshing fallback.
- Use argv arrays locally rather than shell-built command strings.
- Export a manifest with topology and file hashes.
- Validate the returned TDR against exact coordinates, connectivity, contacts,
  and donor/acceptor values.

Dry-run commands:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_pn2d_minimal6_sentaurus_gate -v

D:\msys64\ucrt64\bin\python.exe `
  scripts\run_pn2d_minimal6_sentaurus_gate.py `
  --dry-run `
  --run-id minimal6_gate_prepare_20260713
```

### Step D: Verify the company-computer/VM environment

Before a live run, confirm:

```powershell
Test-Path D:\msys64\ucrt64\bin\python.exe
Test-Path D:\code-repo\vela-tcad\build-release\sentaurus_import.exe
C:\Windows\System32\OpenSSH\ssh.exe sentaurus `
  "test -x /usr/synopsys/sentaurus/O_2018.06-SP2/bin/sdevice && echo sdevice-ok"
```

If the SSH alias is absent on the company computer, configure the equivalent
host explicitly; do not change manifests to hide the environment difference.

### Step E: Run the live Sentaurus topology gate

```powershell
D:\msys64\ucrt64\bin\python.exe `
  scripts\run_pn2d_minimal6_sentaurus_gate.py `
  --run-id minimal6_gate_live_20260713
```

Hard pass criteria for both `sketch` and `mirror`:

- Exactly 6 nodes.
- Exactly 4 triangles.
- Exactly 9 unique edges.
- Exact six approved coordinates within `1e-12 um`.
- Exact triangle connectivity and contact edges.
- Exact donor and acceptor fields, including compensated nodes 2 and 6.
- No SDE, remesh, interpolation, or nearest-node fallback.

If either topology fails, stop after preserving logs and manifests. Document the
consumer incompatibility and do not start Task 3.

### Step F: Continue Tasks 3-6 only after the gate passes

The remaining sequence is:

1. Export exact Sentaurus states for both topologies at `0`, `-12`, and `-19 V`.
2. Add the Vela C++ fixed-state operator evaluator without a nonlinear solve.
3. Add an independent Python mathematical reference and strict joins.
4. Produce exact `36` node rows, `54` edge rows, and `24` triangle rows.
5. Generate seven PNG/PDF static figure pairs.
6. Run the real six-state audit and update BV validation documents.

Never substitute a nearest bias or fabricate complete rows after a partial
Sentaurus run.

## 6. Pause Status Summary

```text
Task 1 implementation: committed
Task 1 tests: passing
Task 1 independent review: paused with 3 unresolved findings
Task 2 implementation: not started
Sentaurus live minimal6 topology gate: not run
Tasks 3-6: not started
```

The highest-priority next action is to harden Task 1 validation and obtain a
clean independent review. The next scientific decision point is the live
Sentaurus topology gate, not a BV convergence run.
