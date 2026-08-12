# BVmethods NMOS original full-physics BV validation (2026-08-12)

## Scope and reference

This validation is limited to the original Sentaurus Applications Library
`GettingStarted/sdevice/BVmethods` fully coupled model:

- Fermi carrier statistics and Old Slotboom band-gap narrowing;
- Masetti `DopingDep` mobility and `HighFieldSaturation(GradQuasiFermi)`;
- `Enormal` (Sentaurus O-2018.06 Enhanced Lombardi);
- `SRH(DopingDep)`, `Band2Band(E2)`, and `Avalanche(Eparallel)`.

The archived Sentaurus O-2018.06-SP2 node-6 voltage-to-current reference is
`6.38318420057198 V` at `abs(Idrain)=1e-4 A/um`.  The independent node-5
series-resistor reference is `6.379791636301563 V`.

## Implemented physics

Vela now supports independent Scharfetter lifetime parameters for electrons
and holes, the total-impurity/net-doping concentration bases, and the
Scharfetter doping-dependent lifetime law in both Gummel and coupled Newton
paths.

The `masetti_field_lombardi` mobility composition implements the O-2018.06
Enhanced Lombardi acoustic and surface-roughness terms, distance damping from
the selected semiconductor/oxide interface, and the local normal electric
field.  It is composed with the PN2D element-edge SG/GSS-Laux avalanche path;
surface-mobility cases currently use the finite-difference Jacobian fallback.

## Vela A/B/C/D ablation

All four voltage-to-current cases used the same mesh, contacts, initial 6 V
state, E2 and avalanche configuration.  Each case converged at all four
points: 6 V, 40 uA/um, 60 uA/um and 100 uA/um.

| Case | SRH | Enormal | V at 100 uA/um (V) | Delta from A (mV) |
|---|---|---:|---:|---:|
| A | constant 1e-7 s | off | 6.397141298 | 0.000 |
| B | doping dependent | off | 6.399088378 | 1.947 |
| C | constant 1e-7 s | on | 6.399945324 | 2.804 |
| D | doping dependent | on | 6.401911823 | 4.771 |

The interaction term `D-B-C+A` is `0.0194 mV`, so the two model increments are
nearly additive at the target current.  Every target row satisfies global
electron and hole continuity closure far below the 1% limit.

## Full-model acceptance

| Method | Sentaurus (V) | Vela (V) | Absolute error | Relative error | Result |
|---|---:|---:|---:|---:|---:|
| Voltage to current | 6.383184201 | 6.401911823 | 18.728 mV | 0.2934% | PASS |
| 10 Mohm um series resistor | 6.379791636 | 6.401809065 | 22.017 mV | 0.3451% | PASS |

The Vela series-resistor result has `Idrain=99.9531 uA/um` and a load-line
residual of `-0.0677 V` (below its `0.1 V` acceptance limit).  It differs from
the Vela voltage-to-current result by only `0.103 mV` (`0.00161%`), confirming
that the extracted BV is independent of the boundary-control implementation.

The original full-physics BV acceptance target (relative error no greater than
2%) is therefore complete and passes.

## Reproducibility and remaining optional evidence

- Generate all Sentaurus decks and Vela configs with
  `scripts/generate_bvmethods_nmos_full_physics_ablations.py`.
- Summarize the completed Vela cases with
  `scripts/summarize_bvmethods_nmos_full_physics_ablations.py`.
- Machine-readable output is
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/sentaurus_full_physics_ablations_20260812/vela/ablation_summary.json`.

The local archive does not contain Sentaurus A/B/C ablation outputs.  Running
those generated decks would add a solver-by-solver decomposition of the SRH
and Enormal increments, but is not required for the already-passing original
full-model BV criterion.  Uploading the prepared mesh, parameter files and
decks to the VM remains intentionally blocked until the user explicitly
authorizes that payload transfer.

## Verification

- `build/test_mobility.exe`: 124 assertions in 25 test cases passed.
- `build/test_recombination.exe`: 72 assertions in 25 test cases passed.
- `build/test_mos_mixed_material.exe '[surface]'`: 11 assertions in 3 test cases passed.
- `build/test_element_edge_gss_laux_avalanche.exe`: 138 assertions in 9 test cases passed.
- Full-physics generator and summary regression tests: 6 tests passed.

An earlier `test_mobility.exe` access violation was caused by stale objects in
the incremental Debug directory after header layout changes.  A clean rebuild
of `build` removed the mixed-ABI executable; the same path now passes the full
test suite without a popup.
