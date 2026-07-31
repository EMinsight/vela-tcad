#!/usr/bin/env python3
"""Verify the M2 single-family substitution and first-Newton-step contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


EXPECTED_BIASES = (-18.0, -19.5, -19.7, -20.0)
EXPECTED_SOURCE_VARIANTS = {
    "vela_baseline",
    "sent_psi_only",
    "sent_qfp_only",
    "sent_density_only",
    "sent_all",
}
EXPECTED_NEWTON_VARIANTS = {
    "vela_baseline",
    "sent_psi_only",
    "sent_qfp_only",
    "sent_all",
    "feedback_baseline",
    "feedback_density_only",
    "feedback_qfp_only",
    "feedback_density_qfp",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def key(row: dict[str, str]) -> tuple[float, str]:
    return float(row["bias_V"]), row["variant"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", type=Path, required=True)
    parser.add_argument("--source-substitution", type=Path, required=True)
    parser.add_argument("--newton-first-update", type=Path, required=True)
    parser.add_argument("--determinism", type=Path, required=True)
    args = parser.parse_args()

    result = json.loads(args.result.read_text(encoding="utf-8-sig"))
    source = read_rows(args.source_substitution)
    newton = read_rows(args.newton_first_update)
    determinism = read_rows(args.determinism)

    require(
        result["schema"] == "vela.pn2d_bv_m2_single_family_state_substitution.v1",
        "unexpected result schema",
    )
    require(result["status"] == "passed", "experiment status is not passed")
    require(result["physics_modified"] is False, "physics was marked modified")
    require(
        result["production_defaults_modified"] is False,
        "production defaults were marked modified",
    )
    require(
        result["acceptance_thresholds_modified"] is False,
        "acceptance thresholds were marked modified",
    )
    require(
        result["fixed_state_coupling_mode"] == "postprocess_only",
        "fixed-state coupling mode changed",
    )
    require(tuple(result["biases_V"]) == EXPECTED_BIASES, "bias lattice changed")

    require(len(source) == 20, "expected 20 fixed-source substitution rows")
    require(len(newton) == 32, "expected 32 first-update rows")
    for bias in EXPECTED_BIASES:
        require(
            {row["variant"] for row in source if float(row["bias_V"]) == bias}
            == EXPECTED_SOURCE_VARIANTS,
            f"wrong source variants at {bias:g} V",
        )
        require(
            {row["variant"] for row in newton if float(row["bias_V"]) == bias}
            == EXPECTED_NEWTON_VARIANTS,
            f"wrong Newton variants at {bias:g} V",
        )

    by_source = {key(row): row for row in source}
    for bias in EXPECTED_BIASES:
        baseline = by_source[(bias, "vela_baseline")]
        density = by_source[(bias, "sent_density_only")]
        full = by_source[(bias, "sent_all")]
        require(
            float(density["source_A_per_um"])
            == float(baseline["source_A_per_um"]),
            f"density changed the fixed SG/Laux source at {bias:g} V",
        )
        require(
            math.isclose(
                float(full["source_to_sentaurus_ratio"]),
                1.0,
                rel_tol=0.005,
                abs_tol=0.0,
            ),
            f"full Sentaurus state does not close within 0.5% at {bias:g} V",
        )

    qfp_wins = 0
    for bias in EXPECTED_BIASES:
        candidates = {
            variant: float(
                by_source[(bias, variant)]["fraction_of_all_sent_error_removal"]
            )
            for variant in (
                "sent_psi_only",
                "sent_qfp_only",
                "sent_density_only",
            )
        }
        qfp_wins += max(candidates, key=candidates.get) == "sent_qfp_only"
    minus20_qfp = float(
        by_source[(-20.0, "sent_qfp_only")][
            "fraction_of_all_sent_error_removal"
        ]
    )
    require(qfp_wins >= 3, "QFP did not win at least three source comparisons")
    require(minus20_qfp >= 0.60, "QFP did not recover 60% at -20 V")

    by_newton = {key(row): row for row in newton}
    density_projection = float(
        by_newton[(-20.0, "feedback_density_only")][
            "qfp_target_projection_fraction"
        ]
    )
    require(
        density_projection < 0.0,
        "density feedback did not move QFP away from the Sentaurus target",
    )
    for bias in EXPECTED_BIASES:
        independent = by_newton[(bias, "vela_baseline")]
        feedback = by_newton[(bias, "feedback_baseline")]
        for field in (
            "initial_combined_residual",
            "trial_combined_residual",
            "step_norm",
        ):
            require(
                math.isclose(
                    float(independent[field]),
                    float(feedback[field]),
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-18,
                ),
                f"baseline probes disagree for {field} at {bias:g} V",
            )

    verdict = result["verdict"]
    require(verdict["source_outcome"] == "qfp_dominant", "wrong source outcome")
    require(
        verdict["first_update_outcome"]
        == "density_feedback_moves_qfp_away_from_sentaurus",
        "wrong first-update outcome",
    )
    require(
        verdict["typed_outcome"]
        == "qfp_dominant__density_feedback_moves_qfp_away_from_sentaurus",
        "wrong typed outcome",
    )

    require(len(determinism) == 120, "expected 120 determinism records")
    require(
        all(row["repeat_count"] == "2" for row in determinism),
        "an artifact is missing a repeat",
    )
    require(
        all(row["unique_hash_count"] == "1" for row in determinism),
        "a repeated artifact has multiple hashes",
    )
    require(
        all(row["byte_identical"] == "1" for row in determinism),
        "not all repeated artifacts are byte-identical",
    )
    require(
        result["determinism"]["all_byte_identical"] is True,
        "result determinism flag is false",
    )

    print(
        json.dumps(
            {
                "status": "pass",
                "typed_outcome": verdict["typed_outcome"],
                "minus20_qfp_recovery_fraction": minus20_qfp,
                "minus20_density_feedback_qfp_projection_fraction": (
                    density_projection
                ),
                "deterministic_artifacts": len(determinism),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
