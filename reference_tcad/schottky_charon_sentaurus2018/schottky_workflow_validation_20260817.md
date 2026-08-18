# Simple Schottky reference closeout

Status: **PASS**

Compared 24 Sentaurus points; maximum error is `0.478652 dex`.

At 1 V, Vela is `0.000108262382 A/um` and Sentaurus is `0.000118824193 A/um` (`8.8886%`).

| Requested bias | Actual bias | Terminal KCL abs / relative | Total current | Result |
| ---: | ---: | ---: | ---: | --- |
| 0 V | 0 V | 2.38e-30 A/um / 7.83e-17 | -3.03556e-14 A/um | pass |
| 0.4 V | 0.4 V | 3.87e-23 A/um / 3.31e-17 | 1.1681e-06 A/um | pass |
| 1 V | 1.00029436 V | 3.35e-21 A/um / 3.09e-17 | 0.000108332 A/um | pass |

No optional Schottky physics is enabled.
