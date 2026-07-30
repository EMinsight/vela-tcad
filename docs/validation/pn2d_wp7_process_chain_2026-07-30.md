# PN2D WP7 paired process-chain analyzer

Date: 2026-07-30

Outcome: `insufficient_observation`.

## Implemented contract

`scripts/analyze_pn2d_bv_process_chain.py` compares:

1. Sentaurus versus Vela avalanche-off;
2. Sentaurus versus Vela IIC/postprocess-only;
3. Sentaurus versus Vela avalanche-on;
4. Sentaurus IIC versus avalanche-on; and
5. Vela postprocess-only versus avalanche-on.

The fixed dependency order is state, density, drive, mobility, current,
alpha, generation, geometric source, residual/Jacobian, first Newton update,
and terminal current. The analyzer emits:

- `stage_summary.csv`;
- `support_summary.csv`;
- `hotspot_chain.csv`;
- `first_departure.json`;
- `source_terminal_closure.csv`;
- `newton_first_update.csv`;
- deterministic process-chain and hotspot SVGs; and
- `acceptance.json`.

It reports scalar signed/relative/log differences, vector component/magnitude
and angle differences, active-support overlap, hotspot, source-weighted
centroid, and cumulative 10/50/90% support.

## Focused acceptance

Synthetic fixtures inject an isolated error at each of all eleven stages.
Every injection is recovered as the first stage at both adjacent knee biases.
Row reordering and unrelated sub-floor tail values do not change the result.
A one-bias injection and a missing simulator both return
`insufficient_observation`.

A causal result also requires:

- every stage in all paired comparisons;
- matching non-implicit supports;
- source and terminal closure at `1e-12`; and
- the same earliest stage at two adjacent biases.

## Current real-data gate

The current run binds the exact 29-point Sentaurus process manifest:

`build-release/pn2d-wp21-full-lattice-a-20260729/manifest.json`

with SHA-256
`190cb08f6c128ce64bdfd9bb8dfc6242bde95238234b7492ef5740b4fa2d3d15`.

There is no matching current-code Vela process-chain manifest containing
avalanche-off, postprocess-only, and avalanche-on records on the same exact
global and knee lattices. The sealed avalanche-on solve still fails at the
request for `-19.693749999999643 V`, so the adjacent required knee states
cannot be fabricated or nearest-matched.

The generated fail-closed report is under
`build-release/pn2d-wp7-process-chain-20260730`. It contains no claimed
causal stage. WP8 and all physics/default changes remain unauthorized.
