# TransportModels DG Frozen-Q oracle

Work point: Vg = 1.0 V, Vd = 2.0 V

Status: **pass**

## Terminal current

| Result | Id (A/um) | Relative error versus Sentaurus |
|---|---:|---:|
| Sentaurus DG | 0.000705525753105 | 0 |
| Vela self-consistent DG baseline | 0.000774966247896 | 9.842376% |
| Vela with Sentaurus Frozen-Q | 0.000714807651367 | 1.315600% |

Frozen-Q removes **8.526775 percentage
points** or **86.63%** of the
self-consistent endpoint current error.

## Interpretation

Classification: **DG equation/field is the dominant endpoint error source**.

The Frozen-Q run changes only the DD variables while preserving the imported
Sentaurus electron quantum potential. The result therefore separates the DG
field/equation contribution from the classical transport and mobility path.
It is a diagnostic oracle, not a production configuration.

The converged Vela self-consistent 2 V state supplies the initial electrostatic
and carrier variables; only its electron quantum potential is replaced by the
Sentaurus value. This avoids attributing an initial-state representation
mismatch to the Frozen-Q experiment.

## Provenance

- Imported restart SHA-256: `1B401AF8B3A8172E4835F6603132DC7370DA0ED3416396B5BF0C8935FC5E2917`
- Hybrid restart SHA-256: `4665D7C01C4B87ABA3935B0E769909ACDA149ABBCA3E63D954DC7953971BB917`
- Config SHA-256: `878050D579CBD320A068D08D4E2AED16188B3DFBA276AF4A0784EDB501688230`
- Final state SHA-256: `080FCB75434EAD3B36125B1C4B8354EB70A5A4BD7D477E0E8E10E5B71EA886E6`
- Curve SHA-256: `A5264981F8C3D7E7913E2B3B43817CC72A933F27F6D26DC28FB88C27711B1438`
