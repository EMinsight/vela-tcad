# PN2D Minimal6 daily-report voltage-axis design

## Objective

Generate three daily-report-only figures from the validated Task 8 comparison
evidence with applied voltage decreasing from left to right. The left side of
each plot represents `-1 V` and the right side represents `-20 V`.

## Scope

The daily-report set contains only:

- `terminal_current.png`
- `maximum_field.png`
- `source_integrals.png`

The sealed Task 8 comparison directory, its figures, report JSON, artifact
hashes, scientific conclusions, and verification contract remain unchanged.
The change is presentation-only and does not alter solver data or derived
quantities.

## Selected approach

Add a reproducible daily-report renderer that reads the already validated
`sweep_comparison.json` and writes the three figures to a separate output
directory. It reuses the established series definitions, units, palette,
legends, disclaimer, dimensions, and accepted-checkpoint data. After plotting,
it explicitly sets the voltage-axis limits in descending display order rather
than changing or negating the underlying voltage values.

The main Task 8 comparison renderer is not modified. Existing figures are not
overwritten or post-processed.

## Data flow

1. Load and validate the Task 8 `sweep_comparison.json` with the existing
   standalone comparison verifier.
2. Select the accepted Vela and Sentaurus records used by the corresponding
   Task 8 figures.
3. Render the three daily-report figures from the original signed bias values.
4. Set each voltage axis so its displayed order is `-1, -2, ..., -20 V` from
   left to right.
5. Write the PNG files and a small manifest to a new daily-report output root.

## Figure contract

Each manifest entry records:

- source comparison report path and SHA-256;
- output PNG SHA-256 and dimensions;
- `x_quantity: applied_bias_V`;
- `x_axis_order: decreasing_left_to_right`;
- finite displayed limits with the left limit greater than the right limit;
- the series identities used by the figure.

The renderer must fail rather than silently emit a misleading figure when the
source report is invalid, the required accepted records are absent, a plotted
value is non-finite, or the resulting axis direction is not decreasing from
left to right.

## Verification

Tests establish the behavior before implementation:

1. A RED contract test requires exactly the three daily-report figure names
   and the decreasing-axis metadata.
2. A rendering test checks that all PNGs decode, have the expected dimensions,
   and carry descending numeric x-axis limits.
3. A source-integrity test confirms the Task 8 comparison report and original
   figure hashes are unchanged.
4. A deterministic rerender test compares pixels from two independent output
   directories.
5. Visual QA confirms that labels, legends, annotations, and the disclaimer
   remain readable at the daily-report image size.

## Acceptance criteria

- The three daily-report PNGs show `-1 V` at the left and `-20 V` at the right.
- No Task 8 artifact is overwritten or rehashed.
- All three figures retain their original data, units, series, and scientific
  meaning.
- Contract, integrity, deterministic-render, and relevant regression tests
  pass.
