# PN2D Minimal6 Sentaurus Topology Gate (2026-07-13)

## Outcome

The live topology gate passed for both canonical six-node topologies on
Sentaurus Device `O-2018.06-SP2`.

| Topology | SDevice | Nodes | Triangles | Unique edges | Contacts | Doping |
|---|---:|---:|---:|---:|---:|---:|
| `sketch` | passed | 6 | 4 | 9 | passed | passed |
| `mirror` | passed | 6 | 4 | 9 | passed | passed |

The source DF-ISE `.grd/.dat` files remain the authoritative explicit mesh and
dataset. Sentaurus Device 2018 does not accept a DF-ISE `.grd` directly through
`File { Grid=... }`, including a `.grd` exported by the same-version `tdx`.
The validated no-remesh path is therefore:

```text
explicit pn2d_minimal6.grd + pn2d_minimal6.dat
  -> tdx -d (format conversion only)
  -> pn2d_minimal6.tdr
  -> SDevice Grid/Doping from the same TDR
```

No SDE, SNMesh, remeshing, interpolation, nearest-node mapping, or separately
generated Vela mesh was used.

## Root-cause evidence

The original direct-grid run failed at:

```text
Reading grid 'pn2d_minimal6.grd' ...
The data file "pn2d_minimal6.grd" was not found,
or its format was not recognized !
```

Same-version `tdx -d` successfully converted that exact file. A control using a
DF-ISE file produced by `tdx -dd` from an accepted TDR failed identically when
passed directly to SDevice, proving this was an SDevice input-entry limitation,
not malformed minimal6 grammar.

The first converted-TDR live retry exposed a UTF-8 BOM in the command deck;
a regression now requires the staged SDevice deck to be BOM-free. The next run
reached the topology gate and showed only one-ULP TDR doping serialization
(`1.0e17` became `1.0000000000000002e17`). Doping validation now uses a
`1.0e-15` relative tolerance with zero absolute tolerance; the regression that
changes a zero dopant value to `1.0` remains rejected.

## Passing command and evidence

```powershell
D:\msys64\ucrt64\bin\python.exe `
  scripts\run_pn2d_minimal6_sentaurus_gate.py `
  --run-id minimal6_gate_live_20260713_tdr3
```

Passing manifest:

```text
build-release/reference_tcad/pn2d_sentaurus2018_minimal6/
  sentaurus_gate_runs/minimal6_gate_live_20260713_tdr3/manifest.json
```

Per-topology gate reports and neutral exports are under:

```text
build-release/reference_tcad/pn2d_sentaurus2018_minimal6/
  sentaurus_gate_runs/minimal6_gate_live_20260713_tdr3/topologies/
    sketch/artifacts/
    mirror/artifacts/
```

Remote evidence remains under:

```text
~/sentaurus_runs/vela_oracle/minimal6_gate_live_20260713_tdr3/{sketch,mirror}
```

## Hash evidence

The uploaded and returned source hashes are identical:

| Topology | File | SHA-256 |
|---|---|---|
| `sketch` | `pn2d_minimal6.grd` | `a3c2d59fc37c80c0848bfb77dafb4239bb2de6a40619ea15a32386a61f5a006d` |
| `mirror` | `pn2d_minimal6.grd` | `4c4ddf24220d14f240a3b28575deabbc4d0da632dfdd9f4829e42e1a334ccafb` |
| both | `pn2d_minimal6.dat` | `c0dcd031db9865afe30d8b541d0d8851892c900628236eaca6c0f4c9180f1135` |

Returned solved TDR hashes:

| Topology | SHA-256 |
|---|---|
| `sketch` | `c454d54e3cd67dfbf4079a7bccdf55d94a559c571dd9f053b70284f3fb5c74e2` |
| `mirror` | `e5ae0360d1ad73b61562af06753859d33509d67787312f158c8529793c77d610` |

## Gate decision

Task 2 is complete. The hard stop blocking Task 3 is cleared: both topology
variants produced returned TDR files and passed coordinate/connectivity/contact/
doping validation. Task 3 may now export the exact `0`, `-12`, and `-19 V`
Sentaurus states on these two immutable topologies.
