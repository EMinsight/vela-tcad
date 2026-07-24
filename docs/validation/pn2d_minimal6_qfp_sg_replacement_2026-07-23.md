# PN2D Minimal6 internal-QFP SG replacement audit

> Quantitative correction: this report used the pre-fix Vela restart-unit and
> fixed-`ni` inverse inputs. Its numerical tables are superseded by
> `pn2d_minimal6_state_unit_and_transport_reaudit_2026-07-23.md`. The original
> report is retained as historical evidence.

## Conclusion

Replacing only the Sentaurus electron and hole quasi-Fermi potentials at
internal nodes 1 and 5 materially improves the Vela directed-edge
Scharfetter-Gummel current comparison, but it does not close the current
residual.

Across the 40 exact states and 280 replacement-affected carrier-edge samples:

| Carrier | Baseline median error | Replaced median error | Paired median improvement | Replaced bounded residual | Sign agreement |
|---|---:|---:|---:|---:|---:|
| Electron | 3.50074 dex | 2.10184 dex | 1.42225 dex | 0.984305 | 100% |
| Hole | 3.52435 dex | 2.07186 dex | 1.47622 dex | 0.983192 | 100% |

The result is therefore classified as
`sign_corrected_but_magnitude_not_closed`. The substituted QFP is causal for
the direction mismatch and part of the magnitude mismatch, but another
operator or support-semantic difference still leaves a median magnitude gap
of roughly 100--130 times.

## Experimental contract

- States: both `sketch` and `mirror` topologies at every integer reverse bias
  from -1 V through -20 V.
- Replaced: Sentaurus `eQuasiFermiPotential` and/or
  `hQuasiFermiPotential` at nodes 1 and 5.
- Frozen: Vela electrostatic potential, stored electron and hole densities,
  production baseline edge mobility, mesh, temperature, and intrinsic
  density.
- SG implementation: production variable-intrinsic-density QFP form with the
  same Bernoulli thresholds and exponent clamp as the C++ code.
- Sentaurus reference: the endpoint-mean node current-density vector projected
  on the canonical edge tangent. It is a same-support proxy, not a native
  Sentaurus SG edge flux.
- Electron current is converted to Vela electron-continuity particle-flux
  convention with `-J_n/q`; hole current uses `+J_p/q`.
- Bounded residual:
  `abs(candidate-reference)/(abs(candidate)+abs(reference))`.

The strict frozen-density SG branch is a negative control. Because it consumes
only frozen `psi`, `n`, `p`, and mobility, changing QFP cannot affect it by
construction.

## Production replay gate

The offline QFP-SG evaluator was replayed against the C++ fixed-state operator
audit for all 720 baseline carrier-edge samples.

| Check | Result |
|---|---:|
| Sample count | 720 |
| Maximum relative difference | 3.07054e-16 |
| Required maximum | 5e-11 |
| Gate | passed |

The production edge mobility was inferred from each nonzero baseline C++ edge
flux using the exact linear dependence of SG flux on mobility. This calibrated
560 carrier-edge samples. The remaining 160 zero-QFP-flux samples use the C++
triangle mobility only as a control; none is incident on replacement nodes 1
or 5.

This step also exposed an earlier diagnostic assumption that must not be used
as production evidence: at about `1e7 V/m`, the Python
`vela_masetti_native_state` reconstruction gave electron mobility near
`9.86e-3 m2/(V s)`, while the production C++ local audit was near
`1.07e-4 m2/(V s)`. The roughly 92-times difference explained the first failed
replay. The triangle-local proxy then differed from the actual SG edge
mobility by up to 0.1626%, so the final experiment froze the mobility inferred
from the baseline edge operator itself.

## Residual structure

The internal vertical edge `(1,5)` responds much more strongly than the
boundary-to-interior edges:

| Edge | Electron replaced median error | Hole replaced median error |
|---|---:|---:|
| `(0,1)` | 2.72795 dex | 2.01848 dex |
| `(0,5)` | 2.67967 dex | 1.97420 dex |
| `(1,2)` | 2.06182 dex | 2.78332 dex |
| `(1,3)` | 2.01756 dex | 2.73506 dex |
| `(1,4)` | 2.67967 dex | 1.97420 dex |
| `(1,5)` | 0.775060 dex | 0.450247 dex |
| `(2,5)` | 2.01756 dex | 2.73506 dex |
| `(3,5)` | 2.06182 dex | 2.78332 dex |
| `(4,5)` | 2.72795 dex | 2.01848 dex |

This pattern is consistent with QFP replacement correcting the internal
driving-force direction while the endpoint-mean node-current projection
remains a poor magnitude surrogate for most boundary-to-interior SG edges.

The conclusion is not driven by one isolated bias:

| Carrier and branch | -1 V | -10 V | -20 V |
|---|---:|---:|---:|
| Electron baseline | 3.62922 dex | 3.49885 dex | 4.06807 dex |
| Electron QFP replaced | 1.47750 dex | 2.06114 dex | 2.12587 dex |
| Hole baseline | 3.64023 dex | 3.51720 dex | 4.08812 dex |
| Hole QFP replaced | 1.32843 dex | 2.01576 dex | 2.10248 dex |

## Evidence

- Artifact root:
  `build-release/pn2d-minimal6-qfp-sg-replacement-20260723-a`
- Full edge samples: `qfp_replacement_edge_samples.csv`
- Summary: `qfp_replacement_summary.csv`
- Independent verification: `independent_verification.json`
- Portable report: `report.html`
- Deterministic report query: `source_query.sql`

Independent verification recomputed all summary statistics from the 3,600
edge-branch rows, checked the exact 40-state support, verified the
opposite-carrier no-op and both-carrier equivalence controls, and confirmed
the non-closure classification.

No production formula under `include/` or `src/` was modified.
