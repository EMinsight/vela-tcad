from __future__ import annotations
import hashlib
import math
import re


def parse_van_overstraeten_de_man(path):
    """Read the carrier-paired vanOverstraetendeMan parameters and source hash."""
    from pathlib import Path
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    match = re.search(r"(?ms)^vanOverstraetendeMan\s+\*.*?^}\s*$", text)
    if match is None:
        raise ValueError("models.par lacks a vanOverstraetendeMan block")
    block = match.group(0)
    pairs = {}
    for key, values in re.findall(r"(?m)^\s*(a\(low\)|a\(high\)|b\(low\)|b\(high\)|E0|hbarOmega)\s*=\s*([^#\n]+)", block):
        tokens = [item.strip() for item in values.split(",")]
        if len(tokens) != 2:
            raise ValueError(f"{key} must contain electron and hole values")
        try:
            pairs[key] = (float(tokens[0]), float(tokens[1]))
        except ValueError as exc:
            raise ValueError(f"{key} contains a non-numeric carrier value") from exc
    required = {"a(low)", "a(high)", "b(low)", "b(high)", "E0", "hbarOmega"}
    missing = sorted(required - set(pairs))
    if missing:
        raise ValueError(f"vanOverstraetendeMan block misses {missing}")
    if pairs["E0"][0] != pairs["E0"][1]:
        raise ValueError("carrier-specific E0 is unsupported by the diagnostic contract")
    return {
        "source": str(source),
        "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "electron": {
            "a_low_cm_inv": pairs["a(low)"][0], "a_high_cm_inv": pairs["a(high)"][0],
            "b_low_v_per_cm": pairs["b(low)"][0], "b_high_v_per_cm": pairs["b(high)"][0],
            "phonon_energy_eV": pairs["hbarOmega"][0],
        },
        "hole": {
            "a_low_cm_inv": pairs["a(low)"][1], "a_high_cm_inv": pairs["a(high)"][1],
            "b_low_v_per_cm": pairs["b(low)"][1], "b_high_v_per_cm": pairs["b(high)"][1],
            "phonon_energy_eV": pairs["hbarOmega"][1],
        },
        "switch_field_v_per_cm": pairs["E0"][0],
    }


def parse_vela_van_overstraeten_defaults(path):
    """Read Vela's production defaults from its C++ configuration header.

    Vela stores inverse lengths and electric fields in SI units.  The
    diagnostic contract converts them to the cm-based units used by the
    tracked Sentaurus parameter file before comparison.
    """
    from pathlib import Path
    source = Path(path)
    text = source.read_text(encoding="utf-8")
    names = {
        "electronALow": ("electron", "a_low_cm_inv"),
        "electronAHigh": ("electron", "a_high_cm_inv"),
        "electronBLow": ("electron", "b_low_v_per_cm"),
        "electronBHigh": ("electron", "b_high_v_per_cm"),
        "holeALow": ("hole", "a_low_cm_inv"),
        "holeAHigh": ("hole", "a_high_cm_inv"),
        "holeBLow": ("hole", "b_low_v_per_cm"),
        "holeBHigh": ("hole", "b_high_v_per_cm"),
    }
    parsed = {"electron": {}, "hole": {}}
    for cpp_name, (carrier, key) in names.items():
        match = re.search(
            rf"\bReal\s+{cpp_name}\s*=\s*([0-9.eE+-]+)\s*;", text
        )
        if match is None:
            raise ValueError(f"Vela production header lacks numeric {cpp_name}")
        parsed[carrier][key] = float(match.group(1)) / 100.0
    switch = re.search(r"\bReal\s+switchField\s*=\s*([0-9.eE+-]+)\s*;", text)
    if switch is None:
        raise ValueError("Vela production header lacks numeric switchField")
    parsed["switch_field_v_per_cm"] = float(switch.group(1)) / 100.0
    phonon = re.search(r"\bReal\s+phononEnergy\s*=\s*([0-9.eE+-]+)\s*;", text)
    if phonon is None:
        raise ValueError("Vela production header lacks numeric phononEnergy")
    phonon_energy_eV = float(phonon.group(1))
    for carrier in ("electron", "hole"):
        parsed[carrier]["phonon_energy_eV"] = phonon_energy_eV
    parsed["source"] = str(source)
    parsed["sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    return parsed


def compare_van_overstraeten_parameters(parsed, production_parameters, *, rel_tol=1.e-12):
    """Compare parsed Sentaurus DeMan coefficients with an explicit Vela table.

    The audit configuration may select a model without serializing numerical
    coefficients.  That is deliberately ``unavailable`` rather than assumed
    equal to Sentaurus defaults.

    Vela's ``aScale`` and ``bScale`` are diagnostic runtime multipliers, not
    coefficients in the Sentaurus material block.  Vela's reference and active
    temperatures are runtime thermal inputs; the Sentaurus block states the
    temperature-factor formula but supplies no numeric ``T0`` or ``T`` value.
    Those controls therefore have no direct default-value comparison here.
    """
    if not isinstance(production_parameters, dict):
        return {"status": "unavailable", "reason": "Vela production configuration does not serialize DeMan coefficients", "comparisons": []}
    keys = (
        "a_low_cm_inv", "a_high_cm_inv", "b_low_v_per_cm",
        "b_high_v_per_cm", "phonon_energy_eV",
    )
    comparisons = []
    try:
        for carrier in ("electron", "hole"):
            for key in keys:
                sentaurus = float(parsed[carrier][key])
                vela = float(production_parameters[carrier][key])
                comparisons.append({"carrier": carrier, "parameter": key, "sentaurus": sentaurus, "vela": vela, "matches": math.isclose(sentaurus, vela, rel_tol=rel_tol, abs_tol=0.0)})
        sentaurus_switch = float(parsed["switch_field_v_per_cm"])
        vela_switch = float(production_parameters["switch_field_v_per_cm"])
        comparisons.append({"carrier": "shared", "parameter": "switch_field_v_per_cm", "sentaurus": sentaurus_switch, "vela": vela_switch, "matches": math.isclose(sentaurus_switch, vela_switch, rel_tol=rel_tol, abs_tol=0.0)})
    except (KeyError, TypeError, ValueError):
        return {"status": "unavailable", "reason": "Vela production configuration lacks a complete numeric DeMan coefficient table", "comparisons": comparisons}
    return {"status": "available" if all(item["matches"] for item in comparisons) else "mismatch", "comparisons": comparisons}

def van_overstraeten_alpha(field_v_per_cm, low_a_cm_inv, low_b_v_per_cm, high_a_cm_inv, high_b_v_per_cm, switch_v_per_cm):
    field = abs(float(field_v_per_cm))
    if field == 0.0: return 0.0
    a,b = (low_a_cm_inv, low_b_v_per_cm) if field < switch_v_per_cm else (high_a_cm_inv, high_b_v_per_cm)
    if a <= 0.0 or b <= 0.0: raise ValueError("Van Overstraeten parameters must be positive")
    return float(a) * math.exp(max(-700.0, -float(b) / field))

def invert_alpha(alpha_cm_inv, a_cm_inv, b_v_per_cm):
    if alpha_cm_inv <= 0.0 or a_cm_inv <= 0.0 or b_v_per_cm <= 0.0 or alpha_cm_inv >= a_cm_inv: raise ValueError("alpha is outside the invertible branch domain")
    return -float(b_v_per_cm) / math.log(float(alpha_cm_inv) / float(a_cm_inv))
def invert_piecewise_alpha(alpha_cm_inv, *, low_a_cm_inv, low_b_v_per_cm, high_a_cm_inv, high_b_v_per_cm, switch_v_per_cm):
    """Return every branch-consistent field candidate for a piecewise alpha law."""
    switch = float(switch_v_per_cm)
    if switch <= 0.0:
        raise ValueError("switch field must be positive")
    candidates = []
    for branch, a, b, is_low in (
        ("low", low_a_cm_inv, low_b_v_per_cm, True),
        ("high", high_a_cm_inv, high_b_v_per_cm, False),
    ):
        try:
            field = invert_alpha(alpha_cm_inv, a, b)
        except ValueError:
            continue
        if (is_low and field < switch) or (not is_low and field >= switch):
            candidates.append({"branch": branch, "field_v_per_cm": field})
    return candidates
def infer_ni_eff(*, psi_V, phin_V, phip_V, n_cm3, p_cm3, thermal_voltage_V):
    """Infer ``ni_eff`` independently from electron and hole relations.

    ``psi_V`` is electrostatic potential, while ``phin_V`` and ``phip_V``
    are the electron and hole quasi-Fermi potentials under the sign convention
    ``n = ni_eff exp((phi_n - psi) / V_T)`` and
    ``p = ni_eff exp((psi - phi_p) / V_T)``.  Therefore the returned estimates
    are ``n exp((psi - phi_n) / V_T)`` and
    ``p exp((phi_p - psi) / V_T)``.  Inconsistent estimates remain separate.
    """
    if thermal_voltage_V <= 0.0 or n_cm3 <= 0.0 or p_cm3 <= 0.0: raise ValueError("densities and thermal voltage must be positive")
    electron = float(n_cm3) * math.exp((float(psi_V) - float(phin_V)) / float(thermal_voltage_V))
    hole = float(p_cm3) * math.exp((float(phip_V) - float(psi_V)) / float(thermal_voltage_V))
    scale = max(abs(electron), abs(hole), 1.0)
    return {"electron_cm3": electron, "hole_cm3": hole, "relative_residual": abs(electron-hole)/scale}