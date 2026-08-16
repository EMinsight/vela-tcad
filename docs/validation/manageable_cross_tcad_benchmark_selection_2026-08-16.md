# Manageable cross-TCAD benchmark selection

## Decision

Use the two-dimensional silicon Schottky diode as the next bounded parity
case, with Charon's `n_diode` case as the canonical source and the Genius and
PISCES II examples as independent scope checks.  Do not begin with another PN
diode or MOS transistor.

The selected first slice is deliberately small: a 1 um by 1 um uniformly
doped n-type silicon slab, one Schottky contact, one Ohmic contact, 300 K,
isothermal two-carrier drift diffusion, SRH recombination, and a 0--1 V DC
sweep.  Barrier lowering, tunnelling, incomplete ionization, AC, self-heating,
and circuit resistance are excluded unless a later executed source deck makes
one of them necessary.

## Candidate script audit

| Source | Executed example | Actual scope | Vela coverage | Selection assessment |
|---|---|---|---|---|
| DEVSIM | `examples/diode/diode_2d.py` | 0.1 um long silicon PN step junction, two-carrier DD, 0--0.5 V | Already covered by Vela PN fixtures and prior Sentaurus PN2D work | Too duplicative; no Schottky example was found in the local DEVSIM example/test tree |
| DEVSIM | `testing/mos_2d.py` | Si/oxide/poly MOS equilibrium plus DD initialization | Mixed material and MOS support exist, but a useful curve would reopen oxide/interface and mobility questions | More uncertain than the selected contact-focused case |
| Genius | `examples/PN_Diode/2D/pn2d.inp` | 3 um by 3 um analytic PN diode, DD, contact resistance, 0--1 V | Mostly already supported | Too duplicative |
| Genius | `examples/MOSCAP/cap.inp` | MOS capacitor with nonuniform doping, EBML3 bias ramp, DDMAC AC sweep | Vela has only quasi-static C-V | Exact AC gap is larger than desired |
| Genius | `examples/Schottky_Diode/sdiode.inp` | Uniform n-Si, Schottky/Ohmic contacts, DD, 0--2 V; source mesh is 3-D | Vela has a 2-D Dirichlet-barrier Gummel prototype | Strong corroborating case, but dimensional reduction and 1 kohm resistance would add ambiguity |
| Charon | `test/nightlyTests/schottky_contacts/n_diode.dd.forward.inp` | 1 um by 1 um n-Si, 1e16 cm^-3, CVFEM DD, SRH, thermionic Schottky flux, 0--1 V | Geometry, material, DD, SRH, sweep and current extraction exist; Newton Schottky/thermionic flux is the isolated gap | **Selected canonical case** |
| PISCES II | `examples/pisces/test08/shot.p2` | 1 um by 1 um p-Si, Schottky equilibrium followed by two-carrier Newton forward sweep | Same contact-family gap, but 5e19 cm^-3 adds degeneracy uncertainty | Independent historical corroboration; defer as polarity/degeneracy extension |
| PISCES II | `examples/pisces/test01/mosIV.p2` | NMOS Id-Vg and Id-Vd with field-dependent mobility | Reopens MOS surface-mobility calibration | Larger and less controlled |

## Canonical Charon semantics

- Geometry: 2-D mapped quad mesh, 1 um by 1 um silicon.
- Doping: uniform donors, `1e16 cm^-3`.
- Contacts: Schottky anode with metal work function 4.75 eV; Ohmic cathode.
- Boundary flux: electron and hole Richardson constants 250 and
  130 A/(cm^2 K^2), without barrier lowering or tunnelling in the selected
  deck.
- Bulk physics: isothermal two-carrier drift diffusion and SRH only.
- Forward scan: Schottky anode 0 to 1 V, continuation step at most 0.05 V.

Charon's implemented thermionic flux is

```text
Jn = An*T^2/Nc * (n - Nc*exp(-(Wf-chi)/kT))
Jp = -Ap*T^2/Nv * (p - Nv*exp((-Eg+Wf-chi)/kT))
```

in A/cm2.  This is a Robin carrier boundary, not the fixed-carrier-density
Dirichlet approximation currently used by Vela's Schottky Gummel prototype.

## Sentaurus translation contract

The committed SDE/SDevice inputs preserve the canonical geometry, doping,
contact work function, two-carrier DD, SRH, temperature, and forward-bias
range.  Charon's Richardson constants are converted to explicit Sentaurus
electron/hole surface recombination velocities at 300 K.  No mobility,
barrier-lowering, tunnelling, high-field, quantum, AC, or thermal model is
added by inference.

The Sentaurus run must finish normally and provide the full 0--1 V current
curve before any Vela implementation work begins.

## Vela development boundary

First run the existing Gummel Dirichlet-barrier model against the frozen
Sentaurus curve.  If it cannot satisfy current sign, monotonicity, equilibrium
current conservation, and log-current shape, the only initially authorized
feature is a thermionic Robin Schottky carrier flux for the coupled Newton
path, with analytic or verified finite-difference Jacobian coverage.  Barrier
lowering, tunnelling, image force, series resistance, and AC remain out of
scope.

## Executed result and bounded acceptance

Sentaurus O-2018.06-SP2 completed the translated 0--1 V sweep in 22.87 s. The
mesh contains 697 vertices and 1280 triangles; the 1 V anode current is
`1.18824193013e-4 A/um`. The equilibrium and forward points are frozen in
`reference_tcad/schottky_charon_sentaurus2018` with source and raw-artifact
hashes.

The required ablation was decisive. Vela's existing Gummel
`dirichlet_barrier` contact failed to converge at 0 V after 150 iterations on
the imported mesh. After adding only the thermionic Robin boundary used by the
Charon and translated SDevice decks, coupled Newton converged at 0 V in five
iterations and produced the expected positive monotone forward curve.

The first-stage comparison is deliberately limited to 0.01--0.54 V:

- 14 Sentaurus reference points compared with log-current interpolation;
- current trend matches;
- maximum log-current discrepancy is 0.478652 dex, below the 0.5 dex gate;
- no current rescaling or bias offset is used.

A diagnostic extension converged through 0.5625 V and stalled at 0.563125 V
with `line_search_non_decrease` dominated by the electron continuity block.
The 0.54--1 V range is therefore a separate numerical-continuation phase. It
does not authorize image force, tunnelling, high-field mobility, AC, series
resistance, or unrelated solver features.

## Phased follow-up gates

| Phase | Scope | Gate | Development authority |
|---|---|---|---|
| S0 | Source and translation freeze | Charon hashes fixed; SDE builds; SDevice reaches 1 V | Complete; no Vela feature |
| S1 | Equilibrium | 0 V converges, finite positive carriers, terminal balance | Complete with thermionic Robin |
| S2 | Low-current branch | 0.01--0.30 V, positive monotone current, <=0.5 dex | Complete |
| S3 | Bounded forward branch | 0.01--0.54 V, 14 reference points, <=0.5 dex | Complete and current acceptance limit |
| S4 | High-forward continuation | 0.54--1 V with the same physics and mesh | Deferred; numerical diagnosis only, no new physics by inference |
| S5 | Polarity extension | Translate PISCES II p-type case after S4 or on a separately bounded nondegenerate variant | Not authorized yet |
| S6 | Resistance/dimensional extension | Genius series resistance or 3-D effects | Not authorized unless an executed deck isolates the need |

DEVSIM's local PN and MOS examples remain useful regression sources, but they
do not justify reopening PN or MOS physics in this Schottky track.
