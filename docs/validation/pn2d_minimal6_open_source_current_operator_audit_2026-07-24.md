# PN2D Minimal6 open-source current-operator audit

Date: 2026-07-24

## Scope

This audit uses three open-source TCAD implementations as independent
current-operator references. It does not modify the Vela production current
formula.

Reference revisions:

- Genius-TCAD-Open: `543da8452d5dfd33e6f8c457f962f6f670f0fce7`
- DEVSIM: `58a9a87083db00c6cadc0b4011c801db2cec5844`
- Charon: `7cc38745625a6011ae3584ed111ec7ee74fb890e`

The numerical comparison uses the sealed 40-state Minimal6 data set
(mirror/sketch, -1 V through -20 V).

## Code-level operator comparison

| Solver | Current support | Driving potential | Mobility support | BGN/DOS treatment |
|---|---|---|---|---|
| Vela | directed mesh edge | endpoint electrostatic and quasi-Fermi potentials | model evaluated on the edge and averaged over adjacent transport cells | variable effective intrinsic density enters the SG potential difference |
| Genius | directed finite-volume edge | effective conduction/valence band difference | endpoint model values are arithmetically averaged | effective band edge contains BGN and effective DOS |
| DEVSIM | directed mesh edge | endpoint electrostatic potential and carrier densities | user-supplied edge model | plain density SG unless BGN is added by the user |
| Charon SGCVFEM | cell edge, later assembled into sub-control volumes | effective potential | edge mobility and diffusion | symmetric BGN and intrinsic-Fermi terms enter the effective potential |
| Charon FEM | integration point | field plus carrier-density gradient | integration-point mobility | BGN is included in the carrier driving field |

Relevant source locations:

- Vela variable-\(n_i\) SG:
  `src/discretization/ScharfetterGummel.cpp:210` and
  `src/discretization/ScharfetterGummel.cpp:420`.
- Vela edge mobility:
  `include/vela/equation/AssemblerUtils.h:620`.
- Genius SG kernels:
  `D:/code-repo/Genius-TCAD-Open/include/math/jflux1.h:145` and
  `D:/code-repo/Genius-TCAD-Open/include/math/jflux1.h:155`.
- Genius effective band-edge use:
  `D:/code-repo/Genius-TCAD-Open/src/solver/ddm1/ddm1_semiconductor.cc:185`
  through `:223`.
- Genius cell quasi-Fermi-gradient mobility driver:
  `D:/code-repo/Genius-TCAD-Open/src/solver/ddm_common/mob_semiconductor.cc:108`.
- DEVSIM electron and hole SG:
  `D:/code-repo/devsim/python_packages/simple_dd.py:36` through `:75`.
- Charon FEM current:
  `D:/code-repo/tcad-charon/src/evaluators/Charon_FEM_CurrentDensity_impl.hpp:15`.
- Charon SGCVFEM effective potential and Bernoulli coefficients:
  `D:/code-repo/tcad-charon/src/evaluators/Charon_SGCVFEM_EdgeCurrDens_impl.hpp:204`
  through `:249`.

For homogeneous silicon, Boltzmann statistics, and symmetric BGN splitting,
the Genius band-edge SG and Charon effective-potential SG reduce to the same
Bernoulli structure as Vela variable-\(n_i\) SG, modulo carrier and edge
orientation conventions. DEVSIM's default expression is the useful control
with the effective-\(n_i\) drift term removed.

## 40-state numerical discrimination

### Node/edge-support counterfactual

Sentaurus node electrostatic and quasi-Fermi potentials were inserted into
Vela. Vela then recomputed carrier density, mobility, and SG current.

| Quantity | Electron median error | Hole median error | Support |
|---|---:|---:|---|
| Recomputed carrier density | 4.27482e-6 dex | 4.25236e-6 dex | all nodes |
| Recomputed mobility | 0.0348089 dex | 0.0320825 dex | identifiable edges |
| Recomputed SG current | 0.397171 dex | 0.474486 dex | Sentaurus endpoint-current proxy edges |
| Current sign agreement | 1.0 | 1.0 | nonzero endpoint-current proxy edges |

This establishes that the remaining current residual is not explained by the
node density calculation or by the median mobility difference.

### Native Sentaurus element support

Sentaurus exports current, mobility, and quasi-Fermi gradient on the same four
native elements, but does not expose native element density or native directed
edge flux.

| Check | Electron | Hole |
|---|---:|---:|
| Median current versus negative-QFP-gradient angle | 0.008824 deg | 0.011670 deg |
| p95 angle | 0.351724 deg | 0.350602 deg |
| Median orthogonal current residual | 0.0001540 | 0.0002037 |
| Median inferred-density gap to arithmetic nodal interpolation | 6.34033 dex | 5.92305 dex |
| Median inferred-density gap to geometric nodal interpolation | 4.96827 dex | 4.02945 dex |

The direction law is therefore closed on native element support. The
magnitude law cannot be closed with a simple interpolation of the exported
node densities.

One exact example is mirror, -1 V, cell 0, electron:

| Field | Value |
|---|---:|
| Native current magnitude | 26.3122 A/m2 |
| Native mobility | 0.0462759 m2/(V s) |
| Native QFP-gradient magnitude | 1.09909e5 V/m |
| Density inferred from \(J=q\mu n|\nabla\phi_F|\) | 3.22894e16 m-3 |
| Arithmetic mean of the three exported node densities | 1.50164e10 m-3 |
| Gap | 6.33249 dex |

### Open-source SG candidates against projected native element current

The native element vectors were averaged over adjacent cells and projected
onto each canonical Vela edge. The same Sentaurus endpoint state and mapped
native mobility were then evaluated with:

- DEVSIM-style plain density SG; and
- Vela/Genius/Charon-style effective-potential or variable-\(n_i\) SG.

| Carrier | Candidate | Median error | p95 error | Sign agreement |
|---|---|---:|---:|---:|
| Electron | DEVSIM plain density SG | 6.15639 dex | 7.14047 dex | 1.0 |
| Electron | effective-potential SG | 6.15631 dex | 7.14047 dex | 1.0 |
| Hole | DEVSIM plain density SG | 6.15519 dex | 7.16067 dex | 1.0 |
| Hole | effective-potential SG | 6.15527 dex | 7.16064 dex | 1.0 |

The median separation between the two SG candidates is only about
0.00008 dex. Selecting the DEVSIM, Genius, Charon, or current Vela Bernoulli
form does not resolve the projected native-element magnitude gap.

Deterministic numerical root:

`build-release/pn2d-minimal6-native-cell-sg-unitfix-open-source-20260724-a`

## Error localization

The evidence ranks the remaining causes as follows:

1. **Dominant and directly demonstrated: support/interpolation mismatch.**
   Native Sentaurus element current and QFP gradient are mutually consistent
   in direction, but the element current magnitude implies a carrier density
   four to six decades away from arithmetic or geometric interpolation of the
   exported node density.
2. **Secondary: node/edge current-support mapping and discrete integration.**
   On the node-derived endpoint proxy, a 0.40-0.47 dex SG current residual
   remains even after potential, density, and mobility replacement. This
   comparison is conditional because the Sentaurus reference is not a native
   directed-edge flux.
3. **Small at the median: mobility model/support.** The corrected mobility
   residual is about 0.03-0.05 dex, too small to explain the current residual.
   Large p95 mobility tails still deserve a separate edge-by-edge check.
4. **Rejected as the leading cause: Bernoulli/SG formula family.** Plain
   density SG and effective-potential SG differ by about 0.00008 dex in this
   experiment, while both retain the same six-decade projected-element gap.
5. **Rejected: carrier sign convention.** All evaluated nonzero current
   comparisons have full sign agreement, and native current is nearly
   collinear with negative QFP gradient.

## Decision and next evidence requirement

The three open-source implementations support the existing Vela SG structure;
they do not justify a production current-formula change.

The next decisive Sentaurus export is one of:

1. native element electron and hole density on exactly the same support as
   current, mobility, and QFP gradient; or
2. native directed-edge electron and hole flux/current for the nine Minimal6
   edges.

Without one of these fields, changing Vela's current formula would fit a
node-to-element support artifact rather than identify a transport law.
