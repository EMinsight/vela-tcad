# TransportModels DG P3 assessment

Status: **complete**. This stage evaluates whether WKB interface closure,
alternative DG operators, or mesh refinement should enter the production
TransportModels comparison.

| Candidate | Executed evidence | Result | Decision |
|---|---|---|---|
| DEVSIM/Garcia-Asenov oxide WKB | Full density-gradient unit suite: 23 cases, 97 assertions; analytic penetration-depth and MOSCAP face-source tests pass | Numerically implemented, but mutually exclusive with the current all-material `include_insulators=true` Sentaurus-box domain | Do not mix into the current baseline; retain as a separate semiconductor-only experiment |
| Neutral continuous Si/SiO2 interface | Three-contract fixed-state sweep | Lowest global L1 and hotspot residual; half-jump and transferred affine calibration are 0.13% and 0.27% worse | Keep neutral continuous |
| `sentaurus_box` | Fixed-state audit plus self-consistent Vg=1 V, Vd=2 V endpoint | Converged in 26 iterations with 0.9596% current error | Keep as production DG operator |
| `p1_direct` | Self-consistent endpoint qualification | Did not qualify within the 19-minute window; last inner residual 9770 | Do not promote |
| P1-lambda / CVFEM / GSS / conservative-sqrt | Common fixed-state manufactured/operator audit | Tests pass, but no candidate has both better residual evidence and a qualified self-consistent current result | Retain as diagnostics |
| Interface-normal mesh refinement | Exact 3315-node Sentaurus/Vela mapping; 20-point normal profiles | First spacing is 0.02114 nm, maximum spacing to 20 nm is 2.15078 nm; source/channel profiles already pass Qn/n targets | Do not refine the whole mesh before fixing the drain-end model; refinement would also break the exact same-mesh oracle |
| Mesh-quality cleanup | Runtime geometry audit | 30 negative-cotangent fallbacks and a minimum angle of 0.1123 degrees remain | A local drain-end quality study is justified later, with a paired Sentaurus remesh |

## Decision

The current production contract remains `sentaurus_box + neutral_continuous +
include_insulators=true`. WKB and alternative operators are functioning research
branches, not evidence-backed fixes for the remaining Id-Vg discrepancy. If a
new mesh study is opened, refine/remesh only the drain-end Si/SiO2/contact corner
and regenerate both Sentaurus and Vela results on the same topology.
