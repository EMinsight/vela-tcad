#!/usr/bin/env python3
"""Verify the M2 Sentaurus-state SG/Laux frozen replay contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


EXPECTED_BIASES = (-18.0, -19.5, -19.7, -20.0)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--source-comparison", type=Path, required=True)
    parser.add_argument("--determinism", type=Path, required=True)
    args = parser.parse_args()

    result = json.loads(args.result.read_text(encoding="utf-8-sig"))
    comparison = rows(args.source_comparison)
    determinism = rows(args.determinism)

    require(
        result["schema"] == "vela.pn2d_bv_m2_sentaurus_frozen_sg_laux.v1",
        "unexpected result schema",
    )
    require(result["observation_only"] is True, "experiment is not observation-only")
    require(result["state_advanced"] is False, "frozen state was marked advanced")
    require(
        result["continuity_feedback_enabled"] is False,
        "continuity feedback was marked enabled",
    )
    require(result["coupling_mode"] == "postprocess_only", "wrong coupling mode")
    require(tuple(result["biases_V"]) == EXPECTED_BIASES, "wrong bias lattice")
    require(
        result["verdict"]["typed_outcome"] == "state_feedback_dominant",
        "typed outcome is not state_feedback_dominant",
    )

    mapping = result["state_mapping"]["per_bias"]
    for bias in EXPECTED_BIASES:
        record = mapping[f"{bias:.12g}"]
        require(record["vela_node_count"] == 115, f"wrong node count at {bias:g} V")
        require(
            record["sentaurus_physical_record_count"] == 122,
            f"wrong Sentaurus physical record count at {bias:g} V",
        )
        require(
            record["excluded_extra_node_count"] == 7
            and record["excluded_coordinate_duplicate_count"] == 7,
            f"unexpected contact-support exclusion at {bias:g} V",
        )
        require(
            record["maximum_coordinate_mismatch_um"] <= 1.0e-12,
            f"coordinate mismatch at {bias:g} V",
        )

        runs = result["runs"][f"{bias:.12g}"]
        require(len(runs) == 2, f"expected two runs at {bias:g} V")
        for run in runs:
            require(run["state_node_count"] == 115, "roundtrip node count changed")
            require(
                run["maximum_state_relative_error"] == 0.0,
                f"state roundtrip changed at {bias:g} V",
            )
            require(
                run["solver_coupled_record_count"] == 0,
                f"solver coupling found at {bias:g} V",
            )
            require(
                run["nonzero_residual_feedback_record_count"] == 0,
                f"residual feedback found at {bias:g} V",
            )
            require(
                run["maximum_qG_closure_relative_error"] <= 1.0e-12,
                f"qG closure failed at {bias:g} V",
            )

    require(len(comparison) == 12, "expected 12 carrier/bias comparison rows")
    total_rows = [row for row in comparison if row["carrier"] == "total"]
    require(len(total_rows) == 4, "expected four total-source rows")
    self_errors = []
    for row in total_rows:
        frozen_ratio = float(row["frozen_to_sentaurus_ratio"])
        frozen_error = float(row["frozen_abs_log10_error_dex"])
        self_error = float(row["self_consistent_abs_log10_error_dex"])
        require(
            math.isclose(frozen_ratio, 1.0, rel_tol=0.005, abs_tol=0.0),
            f"frozen replay is not within 0.5% at {row['bias_V']} V",
        )
        require(frozen_error <= 0.005, "frozen source error exceeds 0.005 dex")
        require(
            self_error - frozen_error >= 0.02,
            "frozen replay did not reduce source error by at least 0.02 dex",
        )
        self_errors.append(self_error)
    require(
        all(b > a for a, b in zip(self_errors, self_errors[1:])),
        "self-consistent source error does not grow monotonically toward -20 V",
    )

    require(len(determinism) == 20, "expected 20 determinism records")
    require(
        all(row["byte_identical"] == "1" for row in determinism),
        "not all repeated artifacts are byte-identical",
    )
    require(
        result["determinism"]["all_artifacts_byte_identical"] is True,
        "result determinism flag is false",
    )
    print(
        json.dumps(
            {
                "status": "pass",
                "typed_outcome": result["verdict"]["typed_outcome"],
                "mean_frozen_abs_log10_error_dex": result["verdict"][
                    "mean_frozen_abs_log10_error_dex"
                ],
                "mean_self_consistent_abs_log10_error_dex": result["verdict"][
                    "mean_self_consistent_abs_log10_error_dex"
                ],
                "mean_error_reduction_dex": result["verdict"][
                    "mean_error_reduction_dex"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
