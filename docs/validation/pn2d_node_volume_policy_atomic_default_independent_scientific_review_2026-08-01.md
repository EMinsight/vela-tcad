# Independent scientific review: PN2D BV atomic default

Date: 2026-08-01

Verdict: **APPROVE**.

No scientific blocking findings remain. The reviewer confirmed that the v2
evidence is prospective, uses actual default renders rather than CLI SG/Laux
injection, completes duplicate 29-point off/IIC/on runs on M0 and M2, preserves
IIC/off state identity, and passes all same-grid Sentaurus curve, knee,
continuity, determinism, and forward-IV gates.

Direct box-measure evidence shows that actual M0/M2 meshes are non-obtuse and
that Sentaurus `Measure[k][j]` agrees with Vela mixed-Voronoi to numerical
precision. The `require_non_obtuse=true` runtime rejection is therefore a
scientifically conservative qualification boundary, not a claim that Vela's
local obtuse half/quarter rule reproduces Sentaurus MixAverage behavior.

Approval is restricted to PN2D BV template version 3 on the verified
non-obtuse M0/M2 Tri3 mesh family. It excludes arbitrary obtuse/non-Delaunay
meshes, other devices and materials, the PN2D IV template, and global C++
defaults.
