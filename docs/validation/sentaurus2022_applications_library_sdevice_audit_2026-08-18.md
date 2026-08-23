# Sentaurus 2022 Applications Library SDevice audit

Date: 2026-08-18
Sentaurus host: `sentaurus` / `TCAD202203`
Sentaurus release: `T-2022.03-SP2`
Library root: `/atctools/Synopsys/tcad/T-2022.03/tcad/T-2022.03-SP2/Applications_Library`

## 1. Outcome

The library is available and was indexed read-only through the existing SSH
connection. It contains 9,356 files in 954 directories. There are 827 `.cmd`
files, 394 standard Workbench `*_des.cmd` SDevice files or fragments, and 267
complete command files that contain both `Electrode {}` and `Solve {}`.

Vela already overlaps the most important steady-state silicon subset:

- two-dimensional Poisson and electron/hole drift diffusion;
- Scharfetter--Gummel transport with Gummel and coupled Newton solution paths;
- Ohmic, metal-gate, fixed-barrier Schottky, and thermionic-Robin Schottky
  contacts;
- Boltzmann and Newton-path Fermi--Dirac statistics;
- Slotboom/OldSlotboom bandgap narrowing;
- constant, Caughey--Thomas, Masetti, high-field, surface, and Enhanced
  Lombardi mobility variants;
- SRH, doping/temperature-dependent SRH lifetime, Auger, local E2 BTBT, and
  Selberherr/Van Overstraeten impact ionization;
- electron density-gradient quantum potential;
- adaptive DC, reverse-BV, ABA/IIC postprocessing, external-series-resistor,
  voltage-to-current, and experimental pseudo-arclength workflows;
- quasi-static terminal charge/CV, state Save/Load, Sentaurus TDR import, and
  extensive field/current diagnostics.

The library makes the next gaps clear. Vela does not currently provide a
physical time-domain continuity solve, AC small-signal admittance, general
MixedMode circuits, carrier/lattice energy equations, explicit bulk/interface
trap occupation, incomplete ionization, radiative or optical generation,
anisotropic/piezoelectric/ferroelectric physics, or calibrated SiC/GaN/GaAs
material libraries. Its BTBT is local E2 rather than Sentaurus `NonlocalPath`,
and its quantum correction is electron density-gradient only.

The recommended order is:

1. refresh the already-supported SingleDevice and BVmethods comparisons under
   T-2022.03-SP2 without replacing the archived 2018 baselines;
2. close the existing BVmethods continuation robustness gap before adding new
   physics;
3. implement implicit AC for one device using `AC1_des.cmd` as the first new
   feature oracle;
4. add a bounded steady-state trap-recombination slice;
5. add backward-Euler device transient before attempting general MixedMode;
6. add thermodynamic self-heating before hydrodynamic transport;
7. only then open SiC/GaN projects, each with one isolated new material feature
   at a time.

## 2. Audit method and limits

The audit used only read operations on the virtual machine:

- `find` for the directory and file inventory;
- standard `*_des.cmd` naming for Workbench SDevice files/fragments;
- the simultaneous presence of `Electrode {}` and `Solve {}` to identify
  complete device command files;
- case-insensitive lexical scans to map model families;
- direct inspection of representative command files and parameter files;
- SHA-256 for a frozen shortlist of future oracle inputs.

Lexical counts below mean that the syntax occurs somewhere in a command file.
They can include a commented alternative, a `Plot` request, or a macro branch
that is inactive for one Workbench node. They are useful for sizing the corpus,
not for claiming that every counted deck executes the feature. The detailed
conclusions use inspected active `Physics`, `Math`, `System`, and `Solve`
sections instead.

No library file was copied, modified, or executed during this audit.

## 3. Corpus inventory

### 3.1 File types

| Item | Count |
| --- | ---: |
| Files | 9,356 |
| Directories | 954 |
| `.cmd` | 827 |
| standard `*_des.cmd` files/fragments | 394 |
| complete command files with `Electrode` and `Solve` | 267 |
| `.par` | 686 |
| `.tcl` | 635 |
| `.tdr` | 75 |
| `.plt` | 1,909 |
| `.project` | 274 |
| `.pdf` | 127 |

### 3.2 Complete SDevice decks by top-level application

| Application | Complete decks |
| --- | ---: |
| GettingStarted | 85 |
| Power | 73 |
| CMOS | 26 |
| Memory | 18 |
| Solar | 11 |
| Variability | 10 |
| Reliability | 8 |
| Bipolar | 8 |
| Templates | 6 |
| FinFET | 5 |
| Opto | 4 |
| Analog | 4 |
| Sensors | 3 |
| Hetero | 3 |
| AdvancedTransport | 2 |
| Backend | 1 |

The full file population is dominated by `GettingStarted` (6,117 files), then
Solar (886), Power (815), CMOS (259), FinFET (210), Memory (182), and Backend
(170). For Vela, `GettingStarted/sdevice` is the best first source because it
contains smaller, more isolated demonstrations; Power and compound-material
projects often activate several missing models at once.

### 3.3 Model-family lexical coverage in `*_des.cmd`

| Syntax or family | Files mentioning it |
| --- | ---: |
| Poisson | 338 |
| Electron / Hole | 323 / 318 |
| Fermi | 263 |
| Effective intrinsic density | 258 |
| Slotboom/OldSlotboom BGN | 274 |
| Mobility | 299 |
| Doping-dependent mobility | 226 |
| High-field saturation | 236 |
| Enormal/Lombardi | 205 |
| SRH | 273 |
| Auger | 227 |
| `eQuantumPotential` | 112 |
| `Avalanche(...)` | 53 |
| `Band2Band(...)` | 32 |
| `Traps(...)` | 106 |
| `Transient(...)` | 119 |
| `ACCoupled(...)` | 27 |
| `Hydrodynamic(...)` | 24 |
| Thermodynamic | 55 |
| `Thermode {}` | 73 |
| `System {}` | 91 |
| Continuation | 6 |
| Incomplete ionization | 35 |
| Radiative recombination | 40 |
| Piezoelectric polarization | 29 |
| Schottky | 46 |

These counts show that the unsupported families are not marginal library
features. Transient, traps, thermal contacts, and circuit/system syntax each
appear in a substantial fraction of the example corpus.

## 4. Representative script semantics

### 4.1 SingleDevice: closest full overlap

Path:
`GettingStarted/sdevice/SingleDevice/SingleDevice_des.cmd`

The deck is a four-terminal silicon MOS example. It activates electron density
gradient, OldSlotboom effective intrinsic density, doping-dependent mobility,
electron/hole high-field saturation driven by quasi-Fermi gradients, Enormal,
and doping/temperature-dependent SRH. The solve graph is:

```text
Poisson + eQuantumPotential initialization
  -> coupled Poisson/Electron/Hole/eQuantumPotential
  -> Save common state
  -> Vd = 0.1 V -> Id-Vg linear sweep
  -> Load common state
  -> Vd = 1.1 V -> Id-Vg saturation sweep
```

This is already the strongest Vela parity case. The archived 2018 workflow
passes both self-consistent branches and the exact-bias field/KCL audit. The
2022 library should therefore be used first as a release-drift oracle, not as
authorization for a new model.

### 4.2 TransportModels: one geometry, four controlled deltas

Paths:

- `GettingStarted/sdevice/TransportModels/IdVgs_des.cmd`
- `GettingStarted/sdevice/TransportModels/IdVds_des.cmd`
- `GettingStarted/sdevice/TransportModels/BV_des.cmd`

Workbench switches the same MOS device among:

- drift diffusion: Poisson/Electron/Hole;
- thermodynamic: adds lattice `Temperature` and `Thermode` surface resistance;
- hydrodynamic: adds electron temperature and lattice temperature and changes
  the high-field/avalanche driving force to `CarrierTempDrive`;
- optional electron density-gradient correction.

The DD plus electron-DG branch is substantially covered by Vela. The
thermodynamic branch isolates a missing lattice heat equation; hydrodynamic is
a later carrier-energy extension. This directory is a better thermal roadmap
than starting from a multi-physics GaN or IGBT example.

#### 4.2.1 Frozen 2022 DD/DG MOS baseline

The first execution stage completed on 2026-08-18 using the unmodified
T-2022.03-SP2 Workbench project. The common SDE mesh and four MOS sweeps are
frozen under
`build-release/reference_tcad/transportmodels_sentaurus2022/run02`:

| Branch | Sweep | Node | Samples | Final point |
| --- | --- | ---: | ---: | ---: |
| DD | Id-Vg, -1 to 2.2 V | 6 | 21 | 1.81040907111e-3 A |
| DD | Id-Vd, 0 to 2 V | 12 | 21 | 8.08310849704e-4 A |
| DD + electron DG | Id-Vg, -1 to 2.2 V | 7 | 21 | 1.68607324102e-3 A |
| DD + electron DG | Id-Vd, 0 to 2 V | 13 | 21 | 7.05525753105e-4 A |

The common mesh contains 1691 points and 3228 elements. At the final sampled
bias, electron DG reduces the drain current by 6.867831% in Id-Vg and
12.716036% in Id-Vd relative to the pure-DD branch. BV and SVisual nodes were
deliberately excluded from this first MOS baseline; one partially launched BV
node is retained only as provenance and is not a valid reference curve.

A diff of the preprocessed command files verifies a controlled experiment:
apart from node-specific output names, the only DD-to-DG changes are enabling
`eQuantumPotential` and adding that equation to every coupled solve.

Workbench node-selection semantics were also confirmed against the local
T-2022.03 Workbench manual: `~` extends a selected node to the root and `+`
forms a union. `~24+~25` therefore selects both complete branches, while
`7+13` selects only the two DG MOS jobs and reuses an existing common mesh.

### 4.3 BVmethods: solver and boundary-control oracle

Path: `GettingStarted/sdevice/BVmethods/BV_des.cmd`

The macro-selected methods are:

| Sentaurus method | Main semantics | Vela mapping |
| --- | --- | --- |
| `ABA_poisson` | Poisson-only avalanche ionization integral | implemented operator/reference workflow |
| `ABA_coupled` | coupled DD plus ionization-integral postprocess | implemented path/current IIC workflows |
| `resistor` | high-voltage source through drain resistor | implemented external series resistor |
| `voltage2current` | voltage ramp, then current-controlled contact | implemented voltage-to-current boundary control |
| `continuation` | trace folded I-V with arc continuation | implemented experimentally; NMOS arc acceptance still pending |
| `transient` | backward-Euler voltage waveform plus resistor | no physical time-domain Vela solve |

The active physics combines Fermi statistics, OldSlotboom BGN,
doping/high-field/Enormal mobility, SRH, `Band2Band(Model=NonLocalPath)`, and
avalanche driven by quasi-Fermi gradients. Vela matches the surrounding DD,
mobility, BGN, SRH, avalanche, and boundary-control families, but its BTBT is a
local E2 source rather than `NonLocalPath`.

The current frozen result accepts ABA, coupled IIC, resistor, and
voltage-to-current mappings. Pseudo-arclength remains a numerical robustness
gap; existing evidence points to state scaling/bordered-corrector robustness,
not a missing physical model.

### 4.4 AC: clean first new feature

Paths:

- `GettingStarted/sdevice/AC/AC1_des.cmd`
- `GettingStarted/sdevice/AC/AC_des.cmd`

Both use ordinary silicon DD physics and a gate DC sweep. `AC1_des.cmd` uses an
implicit AC system; `AC_des.cmd` makes the device and voltage sources explicit
in `System {}`. At each selected gate bias, `ACCoupled` evaluates the
small-signal system at 1 MHz.

Vela's finite-difference quasi-static CV is not equivalent. A faithful first
slice needs the DC Jacobian, the charge/mass derivative, a complex solve (or
real 2x2 block form) for `J + j*omega*M`, and terminal AC current/admittance
extraction. Start with `AC1_des.cmd`; only after it passes should the explicit
source/netlist form be attempted.

### 4.5 MixedMode: defer until device transient exists

Path: `GettingStarted/sdevice/MixedMode/MixedMode_des.cmd`

The example instantiates one MOS device, pulse and DC voltage sources, and a
lumped capacitor. It solves Poisson/electron/hole/contact/circuit equations in
one backward-Euler transient system and records node voltages and branch
currents.

Vela's external series-resistor and voltage-to-current controllers are scalar
outer boundary solves, not a general circuit DAE. General MixedMode therefore
should not be the first transient milestone. Implement and validate a
single-device backward-Euler continuity solve and displacement/charge KCL
before introducing circuit unknowns and netlist elements.

### 4.6 Traps: bounded steady-state model candidate

Paths:

- `GettingStarted/sdevice/Traps/TrapRecombination/sim1_des.cmd`
- `GettingStarted/sdevice/Traps/TrapRecombination/sdevice.par`

The deck is a two-contact device with a Gaussian neutral trap distribution,
electron/hole capture cross sections, forward/reverse bias, and explicit
trapped-charge and gap-state recombination outputs. It is much better for a
first trap implementation than reliability or memory projects.

Vela has fixed region/interface charge and contact surface-recombination
velocities, but neither is an explicit trap occupation/recombination model.
The first implementation should remain steady-state, single-level or Gaussian,
and local; trap dynamics, distributed interface states, degradation, and AC
trap response should remain separate phases.

### 4.7 Wide-bandgap examples: useful later, not minimal first slices

`GettingStarted/sdevice/4H-SiC_PiN/BV_des.cmd` combines:

- a 4H-SiC parameter library;
- SRH, Auger, avalanche, doping/high-field/Enormal mobility;
- split incomplete nitrogen ionization;
- optional anisotropic mobility and Poisson permittivity;
- SiC/oxide fixed interface charge;
- extended precision and a custom iterative linear solver;
- backward-Euler reverse-bias ramp and current break criterion.

`GettingStarted/sdevice/GaN_PiN_Diode/sd_fdiv_des.cmd` and
`sd_rviv_des.cmd` combine:

- GaN material parameters;
- piezoelectric polarization and anisotropic Poisson;
- incomplete ionization;
- interface donor traps;
- SRH, Auger, radiative recombination, and reverse avalanche;
- cylindrical interpretation and extended precision.

Neither is an appropriate one-feature port. A future SiC/GaN program must first
freeze a material-parameter contract, then isolate incomplete ionization,
anisotropic permittivity, polarization, radiative recombination, and avalanche
calibration in separate unit or slab cases before assembling the library deck.

### 4.8 Tunnel diode: nonlocal BTBT target

Paths:

- `Solar/TunnelDiode_Basic_IV/sdevice_des.cmd`
- `Solar/TunnelDiode_Basic_IV/sdevice.par`

The deck activates Fermi statistics, heterointerface treatment, SRH/Auger/
radiative recombination, and `Band2Band(Model=NonlocalPath)`. Vela's E2 model
uses a local electric-field generation law and cannot claim equivalence. This
case should be deferred until a one-dimensional silicon homojunction separates
the nonlocal path integral from heterojunction and radiative effects.

## 5. Vela capability matrix

Status meanings:

- **covered**: executable Vela path and direct tests or cross-TCAD evidence;
- **bounded**: implemented with an explicit limitation or incomplete parity;
- **absent**: no corresponding device equation/model path.

| SDevice capability | Vela status | Boundary of the claim |
| --- | --- | --- |
| 2-D Poisson | covered | finite-volume/box geometry on 2-D meshes |
| electron/hole drift diffusion | covered | isothermal, steady state, SG fluxes |
| Gummel and coupled Newton | covered | adaptive DC; analytic/FD Jacobian options |
| Save/Load and staged bias graph | covered | JSON/CSV state workflow, not SDevice file syntax |
| Quasistationary DC sweep | covered | IV, quasi-static CV, reverse BV |
| Fermi statistics | bounded | coupled Newton; Gummel rejects Fermi--Dirac |
| Slotboom/OldSlotboom BGN | covered | scalar silicon-oriented model |
| doping-dependent mobility | covered | Caughey--Thomas and Masetti variants |
| high-field mobility | covered | electric-field or quasi-Fermi-gradient drive |
| Enormal/Lombardi | bounded | surface selector and field reconstruction are less general than Sentaurus |
| SRH | covered | constant and doping/temperature-dependent lifetime |
| Auger | covered | silicon coefficients/configuration |
| radiative recombination | absent | no mechanism in `RecombinationModel` |
| local E2 BTBT | covered | local field law only |
| `NonlocalPath` BTBT | absent | no path-integral tunnelling model |
| avalanche | covered | Selberherr and Van Overstraeten; self-consistent/postprocess modes |
| ABA/IIC | covered | detailed path and source diagnostics available |
| external resistor | covered | scalar swept-contact boundary controller |
| voltage-to-current switching | covered | scalar outer solve |
| continuation | bounded | core and integration exist; BVmethods NMOS arc remains pending |
| electron density-gradient | covered | electron-only, Boltzmann-compatible outer/frozen coupling |
| hole density-gradient | absent | no hole quantum equation |
| Schrodinger/MLDA/SBTE/MC/NEGF | absent | outside current transport architecture |
| Ohmic/metal gate | covered | 2-D contact boundary models |
| Schottky fixed barrier | covered | prototype Gummel path |
| Schottky thermionic Robin | covered | coupled Newton; no image-force lowering/tunnelling |
| fixed region/interface charge | covered | static charge only |
| explicit trap occupation/DOS | absent | fixed charge is not a trap-state model |
| incomplete ionization | absent | dopants are fully ionized inputs |
| thermodynamic/self-heating | absent | no heat equation or Thermode |
| hydrodynamic carrier energy | absent | no carrier temperature equation |
| physical transient DD | absent | adaptive DC and pseudo-arclength are not time-domain solves |
| AC small signal | absent | quasi-static finite-difference CV only |
| general MixedMode circuit | absent | scalar boundary controllers only |
| optical generation/ray trace/TMM | absent | no optical field/generation pipeline |
| anisotropic tensors | absent | scalar permittivity/mobility material data |
| piezoelectric/ferroelectric polarization | absent | static fixed charge is not polarization dynamics |
| calibrated materials | bounded | built-in Si/SiO2; custom scalar JSON entries allowed |
| SiC/GaN/GaAs library models | absent | no calibrated formula/parameter library |
| 3-D/cylindrical geometry | absent | 2-D triangular mesh with depth scaling |
| configurable ILS/Pardiso/Super parity | absent | Vela uses its own Eigen-based solver stack |

The local source anchors for this matrix include:

- `docs/architecture.md` for solver paths and explicit implementation bounds;
- `include/vela/physics/MobilityModel.h` and
  `src/physics/MobilityModel.cpp` for mobility families;
- `include/vela/physics/RecombinationModel.h` for SRH/Auger and local BTBT;
- `include/vela/physics/BandToBandTunnelingModel.h` for the local E2 contract;
- `include/vela/physics/ImpactIonizationModel.h` for avalanche models and
  driving/source controls;
- `include/vela/physics/DensityGradientQuantumPotential.h` for electron DG;
- `include/vela/boundary/BoundaryCondition.h` for contact models;
- `include/vela/simulation/PseudoArclength.h` for continuation;
- `src/material/MaterialDatabase.cpp` for built-in Si/SiO2 and custom scalar
  material loading.

## 6. Development and comparison queue

### P0: T-2022.03-SP2 oracle refresh, no new physics

Run the 2022 library-derived cases into versioned staging under
`build-release/reference_tcad/...sentaurus2022...`; do not overwrite the
archived `reference_tcad/*sentaurus2018` sources or accepted curves.

1. **SingleDevice**
   - materialize the full Save/Load graph;
   - preserve the linear and saturation point grids;
   - compare current shape, Ion, constant-current Vth, SS, DIBL, exact-bias
     potential/density/eQP fields, and terminal KCL;
   - first determine release drift between 2018 and 2022 before changing Vela.
2. **BVmethods non-transient matrix**
   - run ABA Poisson, ABA coupled, resistor, voltage-to-current, continuation;
   - preserve active model parameters and Workbench macro values in a manifest;
   - compare breakdown/break criteria, complete curves, ionization integrals,
     terminal current, and accepted-step history.
3. **Simple silicon Schottky**
   - rerun the already translated bounded case on T-2022.03-SP2;
   - preserve the no-image-force/no-tunnelling/no-series-resistance contract.

Acceptance rule: if 2022 changes the oracle, freeze and explain that change
before any Vela calibration. Do not silently replace the 2018 accepted ledger.

### P1: close continuation numerics

The existing four non-transient BVmethods mappings pass. The remaining
continuation failure should be treated as a numerical scaling and bordered-
corrector task with physics frozen.

Gate:

- accept at least one usable point beyond the full-physics anchor;
- reproduce the Sentaurus folded branch without mobility/avalanche scaling;
- retain exact mesh, contacts, BGN, SRH, BTBT, avalanche, and source mapping;
- report state-block scaling, tangent normalization, corrector residual, and
  minimum-step termination explicitly.

### P2: implicit AC on the existing DD Jacobian

Use `AC1_des.cmd` first. Keep the device, mesh, DC physics, 1 MHz frequency,
and gate-bias sampling unchanged. Do not add circuit elements in this phase.

Implementation slice:

```text
converged DC state
  -> assemble J = dR/dx
  -> assemble charge/mass derivative M
  -> solve (J + j*w*M) dx = b_ac
  -> extract terminal complex current and Y matrix
```

Gate:

- reproduce conductance and capacitance signs;
- compare selected Y entries across the gate sweep;
- verify charge conservation/column sums and frequency limiting behavior;
- demonstrate that the zero-frequency limit approaches quasi-static charge
  derivatives where that comparison is mathematically valid.

### P3: steady-state trap recombination

Use `Traps/TrapRecombination` but reduce the first Vela port to one trap family
and a two-dimensional slab equivalent. No trap transient or AC response.

Gate:

- equilibrium neutrality including trapped charge;
- forward and reverse I-V for at least three capture-cross-section values;
- local trap occupation in `[0, 1]`;
- integrated electron/hole trap recombination consistency;
- terminal current plus integrated source balance.

### P4: backward-Euler device transient

Start without a circuit netlist. Add carrier accumulation and terminal
displacement/charge current to the coupled Newton system. Reuse the existing
state, charge, adaptive-step, and diagnostics infrastructure where possible.

Gate:

- time-step refinement;
- integrated current equals terminal charge change plus recombination/
  generation balance;
- DC limit approaches the existing steady-state solution;
- only after this passes, add one voltage source and one capacitor from the
  MixedMode example.

### P5: thermodynamic self-heating

Use the `TransportModels` Thermodynamic branch before hydrodynamic transport.
Implement lattice heat conduction, Joule/recombination heat sources, and
Thermode temperature/surface-resistance boundaries while keeping DD carrier
transport unchanged.

Gate:

- isothermal limit exactly recovers the DD curve;
- thermal-energy balance closes;
- peak temperature and Id-Vg/Id-Vd self-heating deltas track Sentaurus;
- mesh and thermal-boundary refinement are reported.

### P6: wide-bandgap material program

Create a material contract before porting complete devices:

1. 4H-SiC scalar equilibrium and DD;
2. incomplete ionization;
3. anisotropic permittivity/mobility;
4. SiC avalanche parameterization;
5. GaN radiative recombination;
6. GaN polarization and interface traps;
7. only then the complete PiN/HEMT examples.

Each stage must have a slab or diode oracle that activates exactly one new
feature. Extended precision and a different linear solver are numerical
requirements to measure, not physical models to infer.

## 7. Reproducibility contract for future library runs

For every selected project, record:

- absolute source path and the relative path under `Applications_Library`;
- Sentaurus banner `T-2022.03-SP2`;
- SHA-256 of every `.cmd`, `.par`, mesh source, and copied TDR used;
- resolved Workbench parameters and active preprocessor branches;
- executed command and environment variables;
- mesh vertex/cell counts, material/region/contact names, and area/depth factor;
- active physics and inactive/commented alternatives separately;
- solver, step, stopping, and break criteria;
- raw PLT/TDR/log artifact hashes;
- normalized current sign, current units, and terminal naming;
- exact comparison point grid and interpolation rule;
- Vela commit, configuration hash, solver path, and diagnostics settings.

The staging convention in `docs/sentaurus_vm_ssh_workflow.md` remains in force:
new T-2022.03-SP2 artifacts go under `build-release/reference_tcad/...` first,
and accepted archived inputs are updated only after explicit review.

## 8. Frozen shortlist hashes

The following hashes were read directly from the T-2022.03-SP2 virtual machine.
Paths are relative to the library root.

| SHA-256 | Relative path |
| --- | --- |
| `ea26c9856f1514b9835625bea1f4d5573122e4173da155b15ece6f10f6feabce` | `GettingStarted/sdevice/SingleDevice/SingleDevice_des.cmd` |
| `9d5fd0a1f038f191b5152baaee687d774de599905c64a6092f74611ab672b1ff` | `GettingStarted/sdevice/SingleDevice/sdevice.par` |
| `96a4be1f68fcc92c4e9c71db17b66772df027563cf2d22b1fe19447de319a2b1` | `GettingStarted/sdevice/SingleDevice/Silicon.par` |
| `5fd4ddc1f22d3abec5af53b226bfb442e3a581b4a39a47907c79255a179fc6be` | `GettingStarted/sdevice/TransportModels/IdVgs_des.cmd` |
| `09cfcebe05aab17194997db5887ba2ae85a8c4ffa3e45caaaded134cab390551` | `GettingStarted/sdevice/TransportModels/IdVds_des.cmd` |
| `6836d63570419e680ad975dd11aac3a3180614ce55e33e132bdf07139b4db70b` | `GettingStarted/sdevice/TransportModels/BV_des.cmd` |
| `a0a3254ee87318c8bf5ef5938ccf1d85cabfe0742607d299e58617ab9135c027` | `GettingStarted/sdevice/BVmethods/BV_des.cmd` |
| `5f3c2d83292c2ae020af09c5e7e73c91be92838b94a3bab56c985c3e02b0414b` | `GettingStarted/sdevice/BVmethods/sdevice.par` |
| `706867b86845bc4bda2b21774217e0ec1fb7faa66008bb4baec570f20efdc20e` | `GettingStarted/sdevice/AC/AC_des.cmd` |
| `ee9499a2d13c12d92ed3bfa3c75aad0ccdff2c068a4f0a79cdad4ee6401c2214` | `GettingStarted/sdevice/AC/AC1_des.cmd` |
| `89236a493873e4acba6615c262c1a6f8f3e5e90a7ab5d27d3def6096dcd39746` | `GettingStarted/sdevice/MixedMode/MixedMode_des.cmd` |
| `e54fe55957ca20b65a23a0d8b3102e64fe9ceeee03ed09d0d18e87e0a1b0cadc` | `GettingStarted/sdevice/Traps/TrapRecombination/sim1_des.cmd` |
| `70747cd8b5c76390538224451b65208d37d6cde23e652e96f41eea1ea27c4775` | `GettingStarted/sdevice/Traps/TrapRecombination/sdevice.par` |
| `719c10e63ff54fdaf7598687e7e00648a613f5810d61676fc4bcb16a69a4b132` | `GettingStarted/sdevice/4H-SiC_PiN/BV_des.cmd` |
| `df4137514169c57ed6581353a2e188432f414176b4d0c3d14f6893f4c767a6c4` | `GettingStarted/sdevice/4H-SiC_PiN/4HSiC.par` |
| `4d6ec0d1e90ef26fd002774c49e74ce86295be194c5e65db45b11726b098a754` | `GettingStarted/sdevice/GaN_PiN_Diode/sd_fdiv_des.cmd` |
| `70947e9ba44ad692e1f7f73d3868b8bc32768bb51038332f300ffee1562eed6b` | `GettingStarted/sdevice/GaN_PiN_Diode/sd_rviv_des.cmd` |
| `4d85e05c0f304bd6e83d3add92adc5f8228a77c428bef408f0f6665154ed13d3` | `GettingStarted/sdevice/GaN_PiN_Diode/sdevice.par` |
| `41fb7872d54e792a3b8d4d5e7faef689bba9633a5cf2b047363e120711507f6c` | `Solar/TunnelDiode_Basic_IV/sdevice_des.cmd` |
| `7bc86d56e34cd8145aed06434215e4c2be000491baf88be0e8f723c183533ae8` | `Solar/TunnelDiode_Basic_IV/sdevice.par` |

The 2022 `SingleDevice/Silicon.par` hash exactly matches the archived local
`reference_tcad/singledevice_sentaurus2018/source/Silicon.par`. This is useful
evidence that the selected silicon parameter excerpt did not drift, but the
command/parameter graph still requires an executed 2022 comparison before a
release-wide equivalence claim.

## 9. P0 execution record

### 9.1 SingleDevice release-drift result: pass

The library-derived SingleDevice case was executed on the VM with the live
SDevice banner `T-2022.03-SP2`. The archived 2018 inputs were not overwritten.
Generated files are staged under:

`build-release/reference_tcad/singledevice_sentaurus2022/sentaurus_vm_runs/singledevice_t2022_refresh_20260818`

Both Save/Load branches completed. Each extracted curve contains the same 21
gate-voltage points from -0.5 V through 2.2 V as the archived 2018 reference.

| Branch | Maximum relative 2018--2022 difference | Maximum log-current difference | Trend |
| --- | ---: | ---: | --- |
| `Vd = 0.1 V` linear | 0.167371% | 0.000727492 dex | match |
| `Vd = 1.1 V` saturation | 0.301357% | 0.00131075 dex | match |

This is small release drift, not a change that warrants Vela recalibration.
The source comparison also found that the SDE deck and `Silicon.par` are byte
identical to the archived copies. The 2022 differences are limited to:

- `sdevice.par` includes `Silicon.par` through Workbench's `@pwd@` token;
- the SDevice `Plot` list uses the newer `ImpactIonization`,
  `eImpactIonization`, and `hImpactIonization` field names in place of the
  older avalanche-generation output names.

`scripts/run_singledevice_sentaurus_vm.py` now resolves the `@pwd@` token,
checks the live Sentaurus release, and records SHA-256 values for both source
and materialized run inputs. The SDE deck already invokes TDX, so the runner no
longer repeats that licensed mesh-conversion command.

### 9.2 Remaining P0 matrix

| Case | Status on 2026-08-18 | Next evidence |
| --- | --- | --- |
| SingleDevice | pass | preserve staged raw outputs and manifest |
| BVmethods | official 2022 project complete; input evolution detected | same-input binary control |
| simple silicon Schottky | running | compare frozen two-file deck with 2018 curve |

### 9.3 BVmethods: official-project evolution, not binary-only drift

The complete 2022 Workbench project finished all 13 submitted nodes: one SDE
node, six SDevice methods, and six SVisual postprocessors. The six method PLT
files were parsed with the same voltage/current/ion-integral criteria as the
2018 archive. The 2022 release renamed the integrated source dataset from
`IntegrSemiconductor AvalancheGeneration` to
`IntegrSemiconductor ImpactIonization`; the analyzer now accepts both names.

| Method | 2018 extracted BV (V) | 2022 project BV (V) | Relative delta |
| --- | ---: | ---: | ---: |
| ABA Poisson | 5.305526 | 8.396796 | +58.265% |
| ABA coupled | 6.377494 | 6.696107 | +4.996% |
| resistor | 6.379792 | 6.673595 | +4.605% |
| voltage-to-current | 6.383184 | 6.673723 | +4.552% |
| continuation | 6.383727 | 6.674713 | +4.558% |
| transient control | 6.378835 | 6.673755 | +4.623% |

These deltas fail the 1% binary-release-drift gate, but they must not be
interpreted as a T-2022 solver regression because the official project inputs
changed materially. Direct source comparison found at least these active
differences:

- junction mesh resolution `Gpn` changed from 0.008 um to 0.002 um, and the
  meshing strategy was revised;
- BTBT changed from `Band2Band(E2)` to
  `Band2Band(Model=NonLocalPath)`;
- avalanche drive changed from `Avalanche(Eparallel)` to
  `Avalanche(GradQuasiFermi)`;
- continuation uses `NewArc`, `Vadapt=2.0`, a 0.25 V maximum voltage step, and
  a 1.3e-3 A/um current limit instead of the archived controller's 0.5 V step
  and 1.443e-3 A/um limit.

The 2022 continuation itself completed with 44 reported curve points, zero
non-smooth-curve cutbacks, and seven failed-solution cutbacks. Its extracted
1e-4 A/um threshold is 6.674713 V. This curve is the correct target when Vela
eventually ports the new official physics bundle; it is not a drop-in
replacement for the existing 2018 Vela acceptance fixture.

The required isolation control is therefore a second run of the archived 2018
mesh/decks on the T-2022.03-SP2 binary. Only that same-input run can measure
binary release drift independently of mesh and model changes.

### 9.4 Simple silicon Schottky same-input result: pass

The frozen two-file Schottky case was run unchanged on T-2022.03-SP2. Both
source and materialized bundle hashes exactly match the archived 2018 hashes:

- `schottky_n_sde.cmd`:
  `c381e21dc3bbc2fd9e56ba39db919c26ad0c8cf596d10fa64a03b631a2d56208`;
- `schottky_n_des.cmd`:
  `2d27fafd43a86501b1420ea3a6c42cf42f35c56198782684fd2162b5afd20d00`.

The complete 24-point nonzero-bias comparison from 0.01 V through 1.0 V is
monotone and agrees with the 2018 curve to a maximum
`5.59e-12 dex`. At 1 V, the 2022 current is
`1.18824193012578e-4 A/um`, versus the archived
`1.18824193013e-4 A/um` (relative difference `3.55e-12`). This is numerical
identity at the precision stored in the archived CSV, so no Schottky model or
Vela calibration change is indicated.

The reproducible runner is `scripts/run_schottky_sentaurus_vm.py`; it checks
the live release banner, records both source and bundle hashes, preserves the
bounded no-image-force/no-tunnelling physics contract, and extracts a normalized
forward-current CSV under `build-release/reference_tcad/schottky_sentaurus2022`.

## 10. P1 continuation work started

The first physics-frozen numerical defect was reproduced and corrected in the
pseudo-arclength core. Previously, the bordered corrector computed
`deltaX = a - z*deltaLambda` and only afterwards applied
`maxParameterUpdate` to `deltaLambda`. Whenever the trust-region cap was
active, the two update components no longer satisfied the same bordered
linearization:

`J*deltaX + F_lambda*deltaLambda = -F`.

A focused regression demonstrates the inconsistency: the old code evaluated a
trial state with `x = 1.0` after capping `deltaLambda = -0.05`, while the
consistent state update is `x = 0.5`. The fix caps `deltaLambda` first and then
recomputes `deltaX` from the capped value. It changes neither the residual nor
the device Jacobian or any physical coefficient.

Verification after the correction:

- pseudo-arclength: 12 test cases, 49 assertions passed;
- NewtonSolver: 84 test cases, 1,172 assertions passed;
- DCSweep: 95 test cases, 3,400 assertions passed;
- selected Python oracle/runner regressions: 17 tests passed.

The end-to-end physics-frozen NMOS rerun used the preserved adjacent accepted
states at 6.0 V and 6.049652 V. It accepted one new point at
`6.049691062499898 V` after 16 corrector iterations and eight shrink retries,
using an accepted arclength step of `3.90625e-5`. The next attempted point
failed after the step shrank below `min_step`.

This closes the linear-consistency defect and passes the first minimum gate of
one usable point beyond the full-physics anchor. It does not yet mark the full
NMOS continuation case accepted: the accepted voltage advance is only about
39 microvolts, and the second point still exhausts the shrink budget. The next
numerical slice must preserve the underlying last corrector failure and
residual instead of reporting only `arclength step shrank below min_step`, then
use that evidence to address state scaling or line-search globalization.
Ultimate acceptance still requires reproduction of the frozen folded branch.

The diagnostic preservation change was implemented immediately after this
run: a minimum-step exit now appends the termination condition to the last
corrector failure instead of overwriting it, and retains the last corrector
iteration count and residual norm. Core, NewtonSolver, and DCSweep regressions
remain green after that change.

### 10.1 Line-search globalization and state metric diagnosis

The identical two-state restart was repeated with last-corrector diagnostics.
At the second-point failure, the current combined residual was
`5.0926240e-7`. All 13 line-search trials were rejected; the best trial used
`alpha = 1.220703125e-4` and produced `5.0930032e-7`. The equation residual
decreased slightly, but the arclength constraint increased and became the
controlling norm. The raw bordered state update had infinity norm `67.6746`;
the component-wise quasi-Fermi limiter changed it by `66.9010` and left norm
`0.773635`. This proved that the limiter, not an insufficient number of
backtracking halvings, had destroyed the bordered descent direction.

The arclength-specific limiter now computes the tightest configured
quasi-Fermi cap and uniformly scales both `deltaX` and `deltaLambda`. Ordinary
Newton retains its existing component-wise limiter. The uniform scaling keeps
both bordered rows at the same first-order decrease. In the next frozen rerun,
line search accepted `alpha = 0.5`; the remaining failure was correctly
classified as exhaustion of the 20-iteration corrector budget rather than a
stale line-search error. A regression now locks that failure-reason behavior.

The restart-state metric then exposed a separate configuration defect. Between
the preserved 6.0 V and 6.049652 V states, the 8,157-component packed-state
difference has scaled L2 norm `62.8042`. With the preparation script's former
`state_weight = 1e-15`, the state contributed only `1.60e-9` of the secant norm
and `lambdaDot = 0.9999999992`; the requested pseudo-arclength controller was
effectively a voltage controller. The script now defaults to zero, which
selects the core's mesh-independent `1/N` state metric, and exposes
`--state-weight` for controlled sensitivity runs. For this restart,
`1/N = 1.22594e-4`, the state fraction is `0.99493`, and
`lambdaDot = 0.07122`.

The bounded `1/N` rerun produced 32 accepted continuation points plus the
anchor in 399 seconds, versus one accepted continuation point before this
work. It reached an apparent Vela
fold at approximately `6.049652697 V`: subsequent arclength motion was almost
entirely in state space and voltage changed by less than a nanovolt. The final
failure was no longer limiter-related (`update_limit_change_norm = 0`); its
residual was `1.0003635e-8`, only 0.036% above the configured `1e-8` tolerance,
and the best trial was 0.0083% higher. This is a near-tolerance stagnation
slice, while the much lower Vela fold relative to the frozen Sentaurus
continuation threshold (`6.383727 V`) is now the dominant physics/discretization
comparison issue.

Evidence directories:

- component-limit diagnosis:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/continuation_line_search_diagnostics_20260818`;
- uniform bordered limiter:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/continuation_uniform_bordered_limit_20260818`;
- mesh-independent state metric:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/continuation_default_state_metric_20260818`.

Verification after this slice: pseudo-arclength 13 test cases / 61 assertions,
NewtonSolver 84 / 1,172, DCSweep 95 / 3,400, and the continuation preparation
script 3 tests all pass.
