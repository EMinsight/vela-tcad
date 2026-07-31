# PN2D Task 10 dose-preserving mesh gate

Date: 2026-07-31

## Superseding correction

The original `N.Window`-only junction construction below is not a valid
finite-element representation of the intended junction position. Although it
preserves the integrated total-impurity dose, it assigns a nonzero net doping
of `+1e17 cm^-3` to the `x=1 um` junction nodes. Linear interpolation therefore
moves the zero-net-doping location away from the geometric junction.

The M0 Vela stall reported below has since been reproduced and removed by the
balanced-half construction (`ND=NA=5e16 cm^-3` at the junction nodes), which
preserves both the `1e17 cm^-3` total impurity and zero net doping. See
`pn2d_task10_m0_stall_root_cause_2026-07-31.md`.

This correction resolves the M0 continuation blocker. It does not retroactively
pass the separate M1-to-M2 source/knee mesh-independence gates; those results
must be regenerated with the corrected junction contract.

The corrected balanced-junction M1/M2 rerun has now been completed. Both
levels reach `-20 V`, but M1 develops a deterministic nonmonotonic
avalanche-on branch in both simulators while M2 remains monotonic. The rerun
still returns primary `mesh_dependent_knee` and secondary
`mesh_dependent_source`. See
`pn2d_task10_balanced_mesh_independence_2026-07-31.md`.

## Decision

Task 10 stops with the typed outcome:

`mesh_dependent_knee`

The secondary failed observation is `mesh_dependent_source`. Task 11 is not
authorized, and the SG/Laux candidate remains opt-in.

## Prospective prerequisite

The frozen dual-domain M0 contract was committed as `2d1c23d`. Two fresh
independent Vela runs then completed the same 29-point off/IIC/on lattice.
All three IV hashes, all physical/non-schedule configuration hashes, and all
87 state hashes matched between the two runs. The frozen contract review
returned `passed` and
`bv_model_consistent_low_current_precision_floor_open`.

## Dose-preserving construction

The rejected mesh sequence placed the P and N constant-profile windows on the
same `x=1 um` boundary. Since SDE rectangle windows include their boundaries,
junction nodes received both active species and their mesh-dependent control
volume changed the discrete total-impurity dose.

The replacement construction assigns the junction to `N.Window` only:

- `P.Window` ends at `0.999 um`;
- `N.Window` starts at `1.000 um`;
- the `0.001 um` gap is 62.5 times smaller than the finest requested
  `0.0625 um` junction spacing;
- no generated mesh has a node inside the gap.

An initial `1e-9 um` epsilon was rejected before physical runs because it was
below the SDE O-2018.06-SP2 geometry tolerance and still double-counted every
junction node.

| Level | Nodes | Triangles | Mesh SHA-256 | Doping SHA-256 | Dose (`cm^-3 um^2`) | Relative to M0 | Double-counted nodes |
| --- | ---: | ---: | --- | --- | ---: | ---: | ---: |
| M0 | 27 | 32 | `c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e` | `fa3cb43da81456a392a6b14a5a50983b60267af7bd52a3787e05ed4c0a7eb27a` | `1.000000000000001e17` | `0` | 0 |
| M1 | 43 | 60 | `9b858d097ad5c07e04640001a4bc144ba8b9e032f479a5dbf504b4a8271f2023` | `b591a7a9a347fc69888bc6756d04ccccaf77aedf9dac7955ea1e028c6b5ff94d` | `1.000000000000001e17` | `0` | 0 |
| M2 | 115 | 191 | `feae2486632ff6f6f45191d41a62fb606cb67e91bcb919a7005995005d8796f9` | `6e5910b2d9b6449cbc83a0a51f501b7509aed7f3650ac9a24a85e6a604d98686` | `1.0000000000000008e17` | `-1.11e-16` | 0 |

The `<0.1%` dose gate and the no-double-counting gate pass.

The source manifest SHA-256 is
`a8d52c139a77f2b9cd0b8316bb457d87a826f724ba23f69e61eb67fa79e237ab`.

## Process-observer corrections

Two observation-only defects were found and fixed before scoring:

1. Remote process artifacts are now tarred once and copied with a single SCP
   connection. This changes evidence transport only.
2. The runtime Tcl probe formerly inferred applied bias from the minimum
   electron QFP over the whole device. On a strong-avalanche state, an
   interior QFP can be below the Anode contact value. The probe now binds bias
   to the mean electron QFP on the minimum-x Anode contact nodes.

The old M1 observer omitted process records at `-19.9`, `-19.95`, and `-20 V`
although all 29 solver snapshots existed. After the contact-bound observer
rerun, M1 has all 29 records. This was not a continuation failure.

The normalized Sentaurus manifests are:

| Level | Manifest SHA-256 | Field records | Aggregate records |
| --- | --- | ---: | ---: |
| M0 | `bffa89ffbf6b0e350cb7816a45e1a7ab1715472a78fcf579ae73954fad090637` | 176,668 | 1,044 |
| M1 | `e505661162d85c705afd7fb8eaca3ae8f626bdeaa313b0544d8695f669783c44` | 302,876 | 1,044 |
| M2 | `2bf863f5390def47e8eb451bd98ede53c40fbaa829464bbf4b53a96f7e89c93a` | 889,256 | 1,044 |

## Failed mesh gates

Sentaurus self-consistent avalanche-on does not preserve a common physical
branch across the two finest meshes:

| Metric | Limit | Observed | Result |
| --- | ---: | ---: | --- |
| M1-to-M2 `V_break` change | `<=0.10 V` | `0.402 V` (`-19.830` to `-19.428 V`) | fail |
| M1-to-M2 maximum integrated-source change | `<2%` | `240316%` at `-19.85 V` | fail |
| M1-to-M2 integrated-source change at `-20 V` | `<2%` | approximately `100%` | fail |
| M1 knee-lattice monotonicity | required for a physical branch | reversals at `-19.0 -> -19.25 V` and `-19.8 -> -19.85 V` | fail |
| M2 knee-lattice monotonicity | required | monotonic | pass |

The branch discontinuity is explicit in the terminal current:

| Bias | M1 Sentaurus on (`A/um`) | M2 Sentaurus on (`A/um`) |
| ---: | ---: | ---: |
| -19.80 V | `-7.16636e-16` | `-1.37520e-15` |
| -19.85 V | `-5.44133e-19` | `-1.52410e-15` |
| -19.90 V | `-3.19203e-5` | `-1.70629e-15` |
| -20.00 V | `-3.65667e-5` | `-2.22680e-15` |

The M1 slope-knee estimator is undefined because the curve does not provide a
stable monotonic knee.

Vela independently fails to establish the dose-preserving M0 on branch:

- avalanche-off: 29/29 exact points;
- IIC postprocess: 29/29 exact points;
- SG/Laux avalanche-on: completes through `-17 V`, then stalls while targeting
  `-18 V`;
- final retained parent bias: `-17.105395800637609 V`;
- final attempted target: `-17.105395832254029 V`;
- attempted step: `-3.1616419704505461e-8 V`;
- recorded continuation attempts: 4,813.

No M1/M2 Vela run is allowed after this first failed level.

## Commands

Representative commands were:

```powershell
python scripts\prepare_pn2d_bv_off_srh_mesh_matrix.py `
  --sealed-source reference_tcad\pn2d_sentaurus2018_coarse7x3\source `
  --out-root build-release\pn2d-task10-dose-preserving-sources-v2-20260731

python scripts\run_pn2d_bv_process_matrix_vm.py `
  --source-root <level-source> `
  --output-root <level-process-root> `
  --remote-root <level-remote-root> `
  --biases 0 -1 ... -19.95 -20

python scripts\run_pn2d_bv_exact_lattice_process.py `
  --runner build-release\vela_example_runner.exe `
  --base-config <frozen-M0-config> `
  --sentaurus-manifest <level-manifest> `
  --mesh-file <paired-mesh.json> `
  --doping-file <paired-doping.csv> `
  --branches avalanche_off,iic_postprocess,avalanche_on `
  --continuation-schedule standard_0p05 `
  --sg-laux-candidate
```

## Implication

The coarse historical SG/Laux agreement is not sufficient to support a
mesh-independent physical or production claim. The failure remains within
the opt-in validation scope. Production defaults are unchanged, generated
simulation outputs and `tmp/` remain untracked, and Task 11 is not entered.
