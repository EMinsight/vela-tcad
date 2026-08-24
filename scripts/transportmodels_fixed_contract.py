#!/usr/bin/env python3
"""Apply and audit the frozen Sentaurus 2022 TransportModels DD/DG contract."""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = (
    REPO / "configs/regression/transportmodels_dd_dg_sentaurus2022_v1.json"
)
DEFAULT_DD_CONTACT_BASIN_CONTRACT = (
    REPO / "configs/regression/transportmodels_dd_contact_basin_v1.json"
)

# Solver keys that enable or parameterize physical models covered by the fixed
# comparison contract.  A contract is closed over these keys: an unlisted model
# must not survive apply_contract(), and adding one afterwards is a validation
# error instead of an unnoticed change to the simulated device.
CONTRACT_PHYSICS_KEYS = frozenset({
    "mobility",
    "recombination",
    "srh_doping_dependence",
    "srh_density_coupling",
    "bandgap_narrowing",
    "carrier_statistics",
    "electron_quantum_potential",
    "auger_cn_m6_per_s",
    "auger_cp_m6_per_s",
    "band_to_band",
    "impact_ionization",
})

ARTIFACT_ROOT_ENV = "VELA_TRANSPORTMODELS_ARTIFACT_ROOT"


def load_contract(path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("schema") != "vela.transportmodels.dd_dg.fixed_contract.v1":
        raise ValueError(f"Unsupported TransportModels contract: {path}")
    return contract


def load_dd_contact_basin_contract(
    path: Path = DEFAULT_DD_CONTACT_BASIN_CONTRACT,
) -> dict[str, Any]:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("schema") != "vela.transportmodels.dd.numerical_contract.v1":
        raise ValueError(f"Unsupported TransportModels DD numerical contract: {path}")
    base_path = (path.parent / contract["base_contract"]).resolve()
    if load_contract(base_path)["contract_id"] != "transportmodels-dd-dg-sentaurus2022-v1":
        raise ValueError(f"Unexpected TransportModels base contract: {base_path}")
    return contract


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_artifact_root(
    explicit: Path | None,
    environ: dict[str, str] | None = None,
) -> Path:
    """Resolve the external TransportModels artifact bundle fail-closed.

    Large generated meshes and restart states are deliberately not committed to
    Git.  Regression entry points must therefore receive a mounted artifact
    bundle explicitly, either through their CLI or through one documented
    environment variable; they must not silently bind to a developer's ignored
    ``build-release`` directory.
    """
    environment = os.environ if environ is None else environ
    raw = explicit if explicit is not None else environment.get(ARTIFACT_ROOT_ENV)
    if raw is None or not str(raw).strip():
        raise ValueError(
            "TransportModels artifact root is required; pass --artifact-root "
            f"or set {ARTIFACT_ROOT_ENV} to a mounted, versioned artifact bundle"
        )
    root = Path(raw).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"TransportModels artifact root does not exist: {root}")
    return root


def strict_json_payload(value: Any) -> tuple[Any, list[str]]:
    """Return an RFC-8259-safe copy and paths replaced with ``null``."""
    nonfinite: list[str] = []

    def convert(item: Any, path: str) -> Any:
        if isinstance(item, float) and not math.isfinite(item):
            nonfinite.append(path)
            return None
        if isinstance(item, dict):
            return {
                str(key): convert(child, f"{path}.{key}")
                for key, child in item.items()
            }
        if isinstance(item, list):
            return [
                convert(child, f"{path}[{index}]")
                for index, child in enumerate(item)
            ]
        if isinstance(item, tuple):
            return [
                convert(child, f"{path}[{index}]")
                for index, child in enumerate(item)
            ]
        return item

    return convert(value, "$"), nonfinite


def strict_json_text(value: Any, *, indent: int = 2) -> str:
    """Serialize JSON without JavaScript-only NaN/Infinity extensions."""
    payload, _ = strict_json_payload(value)
    return json.dumps(payload, indent=indent, allow_nan=False)


def deep_merge(target: dict[str, Any], source: dict[str, Any]) -> None:
    for key, value in source.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            deep_merge(target[key], value)
        else:
            target[key] = json.loads(json.dumps(value))


def materials_path(contract_path: Path = DEFAULT_CONTRACT) -> Path:
    contract = load_contract(contract_path)
    return (contract_path.parent / contract["materials_file"]).resolve()


def apply_contract(
    config: dict[str, Any], branch: str,
    contract_path: Path = DEFAULT_CONTRACT,
) -> dict[str, Any]:
    """Return a copy of *config* with all fixed physical fields overwritten."""
    if branch not in {"dd", "dg"}:
        raise ValueError(f"Unknown TransportModels branch: {branch}")
    contract = load_contract(contract_path)
    result = json.loads(json.dumps(config))
    result["materials_file"] = str(materials_path(contract_path))
    solver = result.setdefault("solver", {})
    allowed_physics = set(contract["common_solver_physics"])
    if branch == "dg":
        allowed_physics.add("electron_quantum_potential")
    for key in CONTRACT_PHYSICS_KEYS - allowed_physics:
        solver.pop(key, None)
    deep_merge(solver, contract["common_solver_physics"])
    deep_merge(solver, contract["workflow_numerics"])
    cold_start = result.get("sweep", {}).get("initialization", {}).get(
        "mode"
    ) == "poisson_block"
    stage_kind = "cold_start" if cold_start else "restarted"
    deep_merge(solver, contract["stage_numerics"][stage_kind])
    if branch == "dg":
        solver["electron_quantum_potential"] = json.loads(json.dumps(
            contract["dg_quantum_contract"]
        ))
    else:
        solver.pop("electron_quantum_potential", None)
    result["fixed_regression_contract"] = {
        "id": contract["contract_id"],
        "path": str(contract_path.resolve()),
        "sha256": sha256(contract_path),
        "materials_sha256": sha256(materials_path(contract_path)),
        "branch": branch,
    }
    return result


def apply_dd_contact_basin_contract(
    config: dict[str, Any],
    contract_path: Path = DEFAULT_DD_CONTACT_BASIN_CONTRACT,
) -> dict[str, Any]:
    """Apply the frozen DD physics plus the contact-basin numerical overlay."""
    contract = load_dd_contact_basin_contract(contract_path)
    base_path = (contract_path.parent / contract["base_contract"]).resolve()
    result = apply_contract(config, "dd", base_path)
    deep_merge(result.setdefault("solver", {}), contract["solver_numerics"])
    result["fixed_dd_numerical_contract"] = {
        "id": contract["contract_id"],
        "path": str(contract_path.resolve()),
        "sha256": sha256(contract_path),
        "base_contract_sha256": sha256(base_path),
    }
    return result


def controlled_solver_delta(
    dd: dict[str, Any], dg: dict[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    dd_solver = json.loads(json.dumps(dd["solver"]))
    dg_solver = json.loads(json.dumps(dg["solver"]))
    dg_solver.pop("electron_quantum_potential", None)
    return dd_solver, dg_solver


def validate_config(
    config: dict[str, Any], branch: str,
    contract_path: Path = DEFAULT_CONTRACT,
) -> list[str]:
    """Return human-readable violations; an empty list means exact compliance."""
    expected = apply_contract(config, branch, contract_path)
    problems: list[str] = []
    if Path(config.get("materials_file", "")).resolve() != materials_path(contract_path):
        problems.append("materials_file does not select the frozen material payload")
    contract = load_contract(contract_path)
    solver = config.get("solver", {})
    allowed_physics = set(contract["common_solver_physics"])
    if branch == "dg":
        allowed_physics.add("electron_quantum_potential")
    for key in sorted((set(solver) & CONTRACT_PHYSICS_KEYS) - allowed_physics):
        problems.append(f"solver.{key} is not allowed by the fixed physics contract")
    for section_name in ("common_solver_physics", "workflow_numerics"):
        for key, value in contract[section_name].items():
            if solver.get(key) != value:
                problems.append(f"solver.{key} differs from {section_name}")
    cold_start = config.get("sweep", {}).get("initialization", {}).get(
        "mode"
    ) == "poisson_block"
    stage_kind = "cold_start" if cold_start else "restarted"
    for key, value in contract["stage_numerics"][stage_kind].items():
        if solver.get(key) != value:
            problems.append(f"solver.{key} differs from stage_numerics.{stage_kind}")
    quantum = solver.get("electron_quantum_potential")
    if branch == "dg" and quantum != contract["dg_quantum_contract"]:
        problems.append("solver.electron_quantum_potential differs from DG contract")
    if branch == "dd" and quantum is not None:
        problems.append("DD must not contain solver.electron_quantum_potential")
    marker = config.get("fixed_regression_contract", {})
    for key, value in expected["fixed_regression_contract"].items():
        if marker.get(key) != value:
            problems.append(f"fixed_regression_contract.{key} is stale or missing")
    return problems


def validate_dd_contact_basin_config(
    config: dict[str, Any],
    contract_path: Path = DEFAULT_DD_CONTACT_BASIN_CONTRACT,
) -> list[str]:
    """Return violations of the physical DD and contact-basin contracts."""
    contract = load_dd_contact_basin_contract(contract_path)
    base_path = (contract_path.parent / contract["base_contract"]).resolve()
    problems = validate_config(config, "dd", base_path)
    solver = config.get("solver", {})
    for key, value in contract["solver_numerics"].items():
        if solver.get(key) != value:
            problems.append(f"solver.{key} differs from DD numerical contract")
    expected = apply_dd_contact_basin_contract(config, contract_path)[
        "fixed_dd_numerical_contract"
    ]
    marker = config.get("fixed_dd_numerical_contract", {})
    for key, value in expected.items():
        if marker.get(key) != value:
            problems.append(f"fixed_dd_numerical_contract.{key} is stale or missing")
    return problems
