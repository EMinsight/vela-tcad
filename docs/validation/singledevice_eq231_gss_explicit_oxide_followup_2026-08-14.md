# SingleDevice Eq. 231 GSS explicit-oxide follow-up (2026-08-14)

## Scope

This follow-up tested whether the GSS DG-DDM discretization can close the
remaining all-material SingleDevice Eq. 231 mismatch.  The imported linear
and saturation Sentaurus states were held fixed for the initial residual
audit.  A full Id-Vg curve was not admitted unless both endpoint gates passed.

The literature contracts used here are:

- GSS Eq. 9.125-9.128: a control-volume sqrt-density flux with a decaying
  exponential branch and a second-order branch in the growing direction;
- GSS explicit-oxide boundary guidance: solve SiO2 as a wide-bandgap material
  without adding a WKB truncation source;
- Sentaurus Device O-2018.06 Eq. 231: the continuous internal-interface state
  is `Phi=Ec+Phi_m+Lambda`, and the unweighted normal driving flux is
  continuous for Formula 0.

The last point is different from a continuous-Lambda GSS interpretation and
is decisive for this deck.

## Implementation

The following opt-in diagnostic operators were added without changing the
qualified default:

- `p1_lambda_direct`: expanded P1 Eq. 231 with continuous Lambda;
- `gss_density_fitted`: Eq. 9.128 flux with continuous Lambda and material-side
  band/DOS traces;
- `gss_potentiallike_fitted`: Eq. 9.128 flux with the existing continuous
  potential-like state;
- `conservative_sqrt_fitted`: exact theta=0.5 sqrt-density weak form with a
  common fixed row scale and the matching Lambda-times-sqrt-density reaction.

The GSS operator uses `expm1` in the decaying branch and the Eq. 9.128
second-order branch in the growing direction.  Its analytic Jacobian is
assembled consistently.  A fitted-edge diagnostic now writes directed edge
jumps, selected branch, stiffness, flux contribution, and Jacobian
contribution.  JSON validation rejects the physically incompatible
`include_insulators: true` plus `oxide_boundary: devsim_wkb` combination.

## Unit qualification

The targeted release-build results are:

```text
test_density_gradient_quantum_potential: 77 assertions, 20 cases passed
test_newton_solver [density_gradient]:   39 assertions, 3 cases passed
```

Manufactured tests cover a continuous Lambda across a material band step, the
small-jump linear sqrt-density limit, the conservative linear auxiliary-field
identity, and the existing stagnation/residual convergence guards.

## Fixed-state and endpoint results

| endpoint | state/operator | max fixed residual | max node | inner result | raw Lambda change |
|---|---|---:|---:|---|---:|
| linear | potential-like + expanded P1 | 2080.200 | 1848 | 500 iterations, fail | 1.3895 V |
| saturation | potential-like + expanded P1 | 2080.192 | 1848 | 500 iterations, fail | 1.4654 V |
| linear | continuous Lambda + expanded P1 | 1043984.226 | 280 | 500 iterations, fail | 0.6016 V |
| saturation | continuous Lambda + expanded P1 | 1044400.004 | 283 | 500 iterations, fail | 1.3823 V |
| linear | continuous Lambda + GSS 9.128 | 1044007.016 | 166 | 178 iterations, converged | 3.1786 V |
| saturation | continuous Lambda + GSS 9.128 | 1044385.470 | 283 | 173 iterations, converged | 3.1790 V |
| linear | potential-like + GSS 9.128 | 952.158 | 30 | 36 iterations, converged | 0.9871 V |
| saturation | potential-like + GSS 9.128 | 952.319 | 30 | 36 iterations, converged | 0.9871 V |
| linear | potential-like + conservative sqrt(n) | 9.191e12 | 2178 | 500 iterations, fail | 1.3097 V |
| saturation | potential-like + conservative sqrt(n) | 9.186e12 | 2178 | 228 iterations, fail | 1.2952 V |

The potential-like GSS variant reduces the maximum fixed-state residual by
about 54%, but it still changes the imported quantum field by nearly 1 V and
therefore fails the reference-state gate.  No self-consistent curve was run.

## Hotspot evidence

For continuous Lambda plus GSS, the linear endpoint hotspot is the
GateOx/Poly interface node 166.  Along the 166-to-167 oxide edge, the imported
Lambda changes from approximately -0.129 V to -3.385 V over about 0.0113 nm.
The fitted half jump is 62.967.  Each of the two adjacent oxide triangles
contributes about 5.219e5 through the quadratic Eq. 9.128 branch, producing the
1.044e6 row residual.

For the potential-like GSS variant, the hotspot moves to PolyReox node 30.
The imported oxide Lambda is approximately -3.256 V.  The decaying GSS flux
saturates and contributes only about -17 in total while the local reaction
contributes about -941 across the adjacent oxide cells.  This is the practical
manifestation of the nonconservative behavior noted in the GSS text.

The conservative sqrt-density weak form is not a remedy: material-side
sqrt-density is discontinuous at the Sentaurus potential-like interface, so a
direct conservative assembly creates a much larger interface residual.

## Decision

The GSS explicit-oxide formulations are rejected as drop-in Sentaurus Eq. 231
parity implementations for this SingleDevice deck.  They remain opt-in
diagnostic operators only.  The evidence supports keeping the existing
continuous potential-like state and the Formula-0 unweighted interface flux.

The next implementation target is a Sentaurus-box-equivalent discretization
of the expanded Eq. 231 quadratic term and interface control volume.  A scalar
theta correction is not sufficient: hotspot rows imply effective values near
0.257-0.275, while the global minimax scalar is about 0.352 and still leaves a
maximum fixed-state residual near 714.  Solver damping, continuous Lambda,
WKB, and ordinary affine-Tri3 CVFEM are no longer supported root-cause paths.

