# SingleDevice linear deep-off current diagnostic (2026-08-16)

## Scope

This short task isolates the only failing point in the self-consistent
SingleDevice Id-Vg comparison:

- `Vds = 0.1 V`;
- `Vg = -0.5 V`;
- Sentaurus drain current `1.2472150e-14 A/um`;
- Vela drain current `1.1103056e-14 A/um`;
- absolute difference `1.3690940e-15 A/um`;
- relative difference `10.9772%` and logarithmic difference `0.05050 decade`.

The Sentaurus state was regenerated with the same physics as the application
library deck and an explicit plot immediately after the drain ramp.  The
tracked command is
`reference_tcad/singledevice_sentaurus2018/source/singledevice_deep_off_sdevice.cmd`.
The TDR was imported on the original 3584-node mesh.  Vela was rerun at the
same point with all-terminal and contact-edge diagnostics, with SRH disabled,
and at `Vds = 0 V` for a numerical-floor control.

## Contact-current decomposition

| Quantity | Sentaurus | Vela | Observation |
| --- | ---: | ---: | --- |
| drain total current | `1.2472150e-14 A/um` | `1.1103056e-14 A/um` | `10.9772%` low |
| drain electron current | `1.2472141e-14 A/um` | `1.1103053e-14 A/um` | `10.9772%` low |
| drain hole current magnitude | `8.9194e-21 A/um` | `3.0584e-21 A/um` | negligible |
| Sentaurus TDR drain contact flux | `1.2429394e-14 A/um` | n/a | only `0.3428%` below its PLT value |
| Vela all-terminal KCL residual | n/a | `2.8198e-17 A/um` | `0.2540%` of drain current |

Using the independent Sentaurus TDR contact flux instead of the PLT total
reduces the Vela difference only to `10.6710%`.  The Sentaurus output convention
and the Vela contact-balance residual are therefore much smaller than the
observed gap.  The total-current error is numerically identical to the electron
current error, so the mismatch is an electron-state/transport discrepancy, not
a hole-current or displacement-current artifact.

## SRH check

The spatial SRH fields are not equal:

| Aggregate | Sentaurus | Vela |
| --- | ---: | ---: |
| maximum absolute SRH rate | `2.1085e24 m^-3 s^-1` | `1.8730e16 m^-3 s^-1` |
| mean absolute rate over all mesh nodes | `2.2078e22 m^-3 s^-1` | `5.2494e14 m^-3 s^-1` |

This is a real field-level parity gap, but it is not the cause of the drain
current failure:

- disabling SRH in Vela changes the drain current by only `2.95e-13` relative;
- the Sentaurus substrate current magnitude is `1.2329e-16 A/um`, only `0.989%`
  of its drain current;
- even treating that entire substrate current as an upper bound on the missing
  SRH contribution accounts for only about `9.0%` of the absolute
  Sentaurus-Vela drain-current difference.

SRH spatial parity can be retained as a separate field-validation item, but no
SRH development is justified solely to close this Id-Vg point.

## Minority-carrier check

In p-type Silicon, where electrons are the relevant minority carrier, Vela's
electron density is systematically lower.  Across 994 p-type nodes the median
log ratio is `-0.09835 decade`.  When weighted by the local Sentaurus electron
current-density magnitude, the geometric mean Vela/Sentaurus electron-density
ratio is `0.93256`, or `6.74%` low.  This explains most of the drain-current
direction and magnitude; the remaining few percent is consistent with the
mobility/SG flux response to that state.

The n-type-region hole population contains near-underflow values and does not
control the drain current.  Its large pointwise relative errors must not be
used as a terminal-current acceptance metric.

## Numerical-floor check

The otherwise identical Vela solve at `Vds = 0 V` gives a drain current of
`3.0015e-23 A/um`, which is `2.70e-9` of the deep-off current.  The
`1.11e-14 A/um` signal is therefore about 8.57 decades above the measured
zero-bias numerical floor.  The discrepancy is not random solver current
noise.

## Decision

The `10.9772%` point is a small but real electron-state/transport mismatch.  It
does **not** come primarily from SRH, contact-current extraction, KCL error, or
the numerical current floor.  At the same time, a strict pointwise `10%`
relative-current gate is brittle at this current level: the point exceeds it by
only `0.9772` percentage point while its absolute error is
`1.3691e-15 A/um` and its logarithmic error is only `0.05050 decade`.

Recommended acceptance policy for `|Id| < 1e-13 A/um` is a hybrid low-current
criterion rather than relative error alone.  This point passes the already
defined `0.2 decade` logarithmic criterion; an explicit absolute-current guard
should be recorded in the validation contract before declaring both curves
fully accepted.  No additional SingleDevice physics implementation is required
for this isolated point.

## Reproducibility

The reusable analyzer is `scripts/analyze_singledevice_deep_off.py`.  It writes
a compact JSON summary and a same-node Silicon CSV containing doping, electron
and hole densities, Sentaurus SRH, and Sentaurus electron-current magnitude.
Runtime artifacts from this execution are under the ignored directory:

`build-release/reference_tcad/singledevice_sentaurus2018/deep_off_20260816`
