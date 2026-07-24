# PN2D Minimal6 high-field box-current Task 1 RED contract

Date: 2026-07-24

Status: `red_confirmed`

## Typed branch contract

- branch:
  `sentaurus_lowfield_element_electric_field`;
- reference label:
  `box_operator_reconstruction`;
- mirror Vela-triangle to Sentaurus-region-cell mapping:
  `0->0, 1->1, 2->2, 3->3`;
- sketch Vela-triangle to Sentaurus-region-cell mapping:
  `0->0, 1->3, 2->2, 3->1`; and
- geometric-zero edges remain typed and do not receive a dex value.

Each sample must record topology, bias, carrier, Vela triangle id, Sentaurus
region-cell id, low-field source hash, electric-field source hash, high-field
parameters, element-edge box coefficient, typed status, and reconstruction
label.

## Sealed input hashes

| Input | SHA-256 |
|---|---|
| native element decomposition | `1f9a34a1cf06b22b45bcc640f50aa37bdb6c03bed55a78cfe199f57b69ffadb3` |
| triangle mobility decomposition | `5d2ded369bcb51c11360bd48ae4e54245b08f35fe27eaf05739fdc77795a2d1e` |
| 40-state transport elements | `7f29dd652b4d8e2c01103ad49bf50b03cc6cfa75977f8332559681f94beccad3` |
| mirror low-field electron cells | `bd28819fd4a6d7265a351fb26b1a73b4c9e705b9ee9cdba206e36142d5657821` |
| mirror low-field hole cells | `ad749779cb77197388a5729f88480464b40303f0a67f13b2ad63bcc42f6a9546` |
| sketch low-field electron cells | `1f10bf3bbb4f35e5418ebddb48c00bc89d5f9b911c71bdc3c116ce447329852d` |
| sketch low-field hole cells | `090897138e5819d18936fb51bf4d68d5c03c84cf6ea9cb85c79a3e08ade64c75` |

## Expected RED

Command:

`python -m unittest tests.regression.test_pn2d_minimal6_highfield_box_current`

Observed failure:

`ModuleNotFoundError: No module named
'scripts.pn2d_minimal6_diagnostics.highfield_box_replay'`

The failure is the intended missing-capability RED. It is not caused by data
availability, numerical thresholds, locale, path formatting, or newline
handling.

Production files modified by Task 1: `false`.
