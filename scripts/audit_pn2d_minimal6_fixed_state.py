#!/usr/bin/env python3
"""Independent report generator for Task3 PN2D minimal6 fixed-state roots."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
from pathlib import Path
import subprocess
import tempfile
from types import SimpleNamespace

import numpy as np

SCHEMA = "vela.pn2d_minimal6_fixed_state_audit.v1"
TASK3_SCHEMA = "vela.pn2d_minimal6_states.v1"
BIASES = (0.0, -12.0, -19.0)
TOPOLOGIES = ("sketch", "mirror")
Q_C = 1.602176634e-19
KB_J_PER_K = 1.380649e-23
TEMPERATURE_K = 300.0
THERMAL_VOLTAGE_V = KB_J_PER_K * TEMPERATURE_K / Q_C
INTRINSIC_DENSITY_M3 = 1.0e16
FORMULA_LIMIT = 5.0e-12
STATE_LIMIT = 1.0e-12
COORDINATE_LIMIT_UM = 1.0e-12
GEOMETRIC_ZERO_AREA_M2 = 1.0e-27
GEOMETRIC_ZERO_SOURCE_PER_M_S = 1.0e-285
TASK4_PRODUCER = "build-release/pn2d_minimal6_operator_audit.exe"
TASK4_SOURCE_COMMIT = "37a95459dc5f360bb24b9afa00439301935e98de"
REPO_ROOT = Path(__file__).resolve().parents[1]
SILICON_ELECTRON_MOBILITY_M2_PER_V_S = 0.135
SILICON_HOLE_MOBILITY_M2_PER_V_S = 0.048
VAN_OVERSTRAETEN = {
    "electron": {"a_low": 7.03e7, "a_high": 7.03e7, "b_low": 1.231e8, "b_high": 1.231e8},
    "hole": {"a_low": 1.582e8, "a_high": 6.71e7, "b_low": 2.036e8, "b_high": 1.693e8},
}
VAN_OVERSTRAETEN_SWITCH_FIELD_V_PER_M = 4.0e7
NODES = {1: (0.0, 0.5), 2: (1.0, 0.5), 3: (2.0, 0.5), 4: (2.0, 0.0), 5: (0.0, 0.0), 6: (1.0, 0.0)}
TRIS = {
    "sketch": ((1, 5, 2), (5, 6, 2), (2, 6, 4), (2, 4, 3)),
    "mirror": ((1, 5, 6), (1, 6, 2), (2, 6, 3), (6, 4, 3)),
}
CONTACTS = {"Anode": (1, 5), "Cathode": (3, 4)}
DONORS_CM3 = {1: 0.0, 2: 1e17, 3: 1e17, 4: 1e17, 5: 0.0, 6: 1e17}
ACCEPTORS_CM3 = {1: 1e17, 2: 1e17, 3: 0.0, 4: 0.0, 5: 1e17, 6: 1e17}
FIELD_CONTRACT = {
    "ElectrostaticPotential": (1, "V"),
    "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"),
    "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"),
    "DonorConcentration": (1, "cm^-3"),
    "AcceptorConcentration": (1, "cm^-3"),
    "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"),
    "hCurrentDensity": (2, "A*cm^-2"),
    "eMobility": (1, "cm^2*V^-1*s^-1"),
    "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eAlphaAvalanche": (1, "cm^-1"),
    "hAlphaAvalanche": (1, "cm^-1"),
}
NODE_HEADER = ("node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3")
EDGE_HEADER = (
    "edge_id", "node0", "node1", "length_m",
    "electron_raw_signed_flux_per_m2_s", "hole_raw_signed_flux_per_m2_s",
    "electron_midpoint_density_m3", "hole_midpoint_density_m3",
    "electron_impact_field_V_per_m", "hole_impact_field_V_per_m",
    "electron_alpha_per_m", "hole_alpha_per_m", "edge_area_m2",
)
LOCAL_SUFFIXES = (
    "edge_id", "node0", "node1", "truncated_partial_volume_m2",
    "electron_cell_qf_field_V_per_m", "hole_cell_qf_field_V_per_m",
    "electron_midpoint_density_m3", "hole_midpoint_density_m3",
    "electron_mobility_m2_per_V_s", "hole_mobility_m2_per_V_s",
    "electron_alpha_per_m", "hole_alpha_per_m",
    "electron_flux_proxy_per_m2_s", "hole_flux_proxy_per_m2_s",
    "electron_source_integral_per_m_s", "hole_source_integral_per_m_s",
)
TRIANGLE_HEADER = (
    "cell_id", "node0", "node1", "node2", "signed_double_area_m2",
    "grad_psi_x_V_per_m", "grad_psi_y_V_per_m",
    "grad_phin_x_V_per_m", "grad_phin_y_V_per_m",
    "grad_phip_x_V_per_m", "grad_phip_y_V_per_m",
) + tuple(f"local_edge{i}_{suffix}" for i in range(3) for suffix in LOCAL_SUFFIXES)


class ContractError(RuntimeError):
    pass


def finite(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ContractError(f"invalid numeric value for {label}") from error
    if not math.isfinite(result):
        raise ContractError(f"non-finite value for {label}")
    return result


def hybrid_error(actual, expected, abs_floor=1e-300):
    return abs(actual - expected) / max(abs(actual), abs(expected), abs_floor)


def absolute_error(actual, expected):
    return abs(actual - expected)


def gate(actual, expected, limit, label, abs_floor=1e-300):
    error = hybrid_error(actual, expected, abs_floor)
    if error >= limit:
        raise ContractError(f"{label} {error:.17g} >= {limit:.17g}")
    return error


def zero_normalized_gate(actual, expected, threshold, limit, label):
    if abs(actual) <= threshold and abs(expected) <= threshold:
        return 0.0
    return gate(actual, expected, limit, label)


def geometric_source_gate(actual, expected, actual_partial, expected_partial, limit, label):
    if (abs(actual_partial) <= GEOMETRIC_ZERO_AREA_M2
            and abs(expected_partial) <= GEOMETRIC_ZERO_AREA_M2):
        return zero_normalized_gate(
            actual, expected, GEOMETRIC_ZERO_SOURCE_PER_M_S, limit, label)
    return gate(actual, expected, limit, label)


def van_overstraeten_alpha(field_V_per_m, carrier, temperature_K=TEMPERATURE_K,
                            minimum_field_V_per_m=0.0):
    field = abs(finite(field_V_per_m, "Van Overstraeten field"))
    if field < minimum_field_V_per_m or field <= 0.0:
        return 0.0
    if carrier not in VAN_OVERSTRAETEN:
        raise ContractError("Van Overstraeten carrier must be electron or hole")
    kb_eV_per_K = KB_J_PER_K / Q_C
    phonon_energy_eV = 0.063
    reference_temperature_K = 300.0
    gamma = (math.tanh(phonon_energy_eV / (2.0 * kb_eV_per_K * reference_temperature_K))
             / math.tanh(phonon_energy_eV / (2.0 * kb_eV_per_K * temperature_K)))
    params = VAN_OVERSTRAETEN[carrier]
    suffix = "low" if field < VAN_OVERSTRAETEN_SWITCH_FIELD_V_PER_M else "high"
    prefactor = params[f"a_{suffix}"]
    critical = params[f"b_{suffix}"]
    return gamma * prefactor * math.exp(max(-700.0, min(0.0, -critical * gamma / field)))


def classify_orientation_pair(sketch, mirror):
    result = {"mirror_over_sketch": None, "signed_difference": mirror - sketch,
              "absolute_log10_ratio": None}
    if sketch == 0 and mirror == 0:
        result["zero_classification"] = "both_zero"
    elif sketch == 0:
        result["zero_classification"] = "sketch_zero"
    elif mirror == 0:
        result.update(mirror_over_sketch=0.0, zero_classification="mirror_zero")
    else:
        ratio = mirror / sketch
        result.update(mirror_over_sketch=ratio,
                      absolute_log10_ratio=abs(math.log10(abs(ratio))),
                      zero_classification="neither_zero")
    return result


def bernoulli(value):
    if not math.isfinite(value):
        raise ContractError("non-finite Bernoulli argument")
    if abs(value) < 1e-8:
        return 1.0 - value / 2.0 + value * value / 12.0 - value**4 / 720.0
    if value > 50.0:
        exp_negative = math.exp(-value) if value < 745.0 else 0.0
        return value * exp_negative / (1.0 - exp_negative) if exp_negative else 0.0
    if value < -50.0:
        exp_value = math.exp(value)
        return -value / (1.0 - exp_value)
    return value / math.expm1(value)


def limited_exp(value):
    return math.exp(max(-500.0, min(500.0, value)))


def sg_electron_variable_ni_flux(ni0, ni1, psi0, psi1, phin0, phin1, vt, coef,
                                 include_ni_gradient_drift=True):
    if phin0 == phin1:
        return 0.0
    eta = (psi1 - psi0) / vt
    if ni0 > 0 and ni1 > 0 and include_ni_gradient_drift:
        eta += math.log(ni1 / ni0)
    n0 = ni0 * limited_exp((psi0 - phin0) / vt)
    n1 = ni1 * limited_exp((psi1 - phin1) / vt)
    return coef * (bernoulli(-eta) * n0 - bernoulli(eta) * n1)


def sg_hole_variable_ni_flux(ni0, ni1, psi0, psi1, phip0, phip1, vt, coef,
                             include_ni_gradient_drift=True):
    if phip0 == phip1:
        return 0.0
    eta = (psi1 - psi0) / vt
    if ni0 > 0 and ni1 > 0 and include_ni_gradient_drift:
        eta += math.log(ni0 / ni1)
    p0 = ni0 * limited_exp((phip0 - psi0) / vt)
    p1 = ni1 * limited_exp((phip1 - psi1) / vt)
    return coef * (bernoulli(eta) * p0 - bernoulli(-eta) * p1)


def sg_electron_flux(n0, n1, dpsi, vt, mobility, length):
    return -mobility * vt / length * (bernoulli(-dpsi / vt) * n0 - bernoulli(dpsi / vt) * n1)


def sg_hole_flux(p0, p1, dpsi, vt, mobility, length):
    return -mobility * vt / length * (bernoulli(dpsi / vt) * p0 - bernoulli(-dpsi / vt) * p1)

def aux2(value):
    if value >= 0:
        exp_negative = math.exp(-value) if value < 745 else 0.0
        return exp_negative / (1.0 + exp_negative)
    exp_value = math.exp(value)
    return 1.0 / (1.0 + exp_value)


def gss_logistic_midpoint(d0, d1, v0, v1, vt, carrier):
    if carrier == "electron":
        value = (v1 - v0) / (2.0 * vt)
    elif carrier == "hole":
        value = (v0 - v1) / (2.0 * vt)
    else:
        raise ContractError("carrier must be electron or hole")
    return d0 * aux2(value) + d1 * aux2(-value)


def triangle_gradient(points, values):
    matrix = np.array([[1.0, float(x), float(y)] for x, y in points])
    if abs(np.linalg.det(matrix)) <= 1e-300:
        raise ContractError("degenerate triangle")
    coefficients = np.linalg.inv(matrix) @ np.array(values)
    return float(coefficients[1]), float(coefficients[2])


def area2(points):
    return ((points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1]))


def canonical_projection(vector, point0, point1):
    dx = point1[0] - point0[0]; dy = point1[1] - point0[1]
    length = math.hypot(dx, dy)
    if length <= 1e-300:
        raise ContractError("zero-length projection edge")
    return (vector[0] * dx + vector[1] * dy) / length


def genius_truncated_partial_volume(points, local_edge):
    points = [np.array(point, dtype=float) for point in points]
    determinant = area2(points)
    if abs(determinant) <= 1e-300:
        return 0.0
    squared = [float(point @ point) for point in points]
    center = np.array([
        (squared[0]*(points[1][1]-points[2][1]) + squared[1]*(points[2][1]-points[0][1]) + squared[2]*(points[0][1]-points[1][1]))/(2*determinant),
        (squared[0]*(points[2][0]-points[1][0]) + squared[1]*(points[0][0]-points[2][0]) + squared[2]*(points[1][0]-points[0][0]))/(2*determinant),
    ])
    sides = ((0,1),(1,2),(2,0)); lengths=[]; dual=[]; obtuse=-1
    for index,(i,j) in enumerate(sides):
        opposite=(2+index)%3; lengths.append(float(np.linalg.norm(points[i]-points[j])))
        distance=float(np.linalg.norm((points[i]+points[j])/2-center))
        if float((points[i]-points[opposite]) @ (points[j]-points[opposite])) < 0:
            dual.append(-distance); obtuse=index
        else:
            dual.append(distance)
    if obtuse >= 0:
        i,j=sides[obtuse]; opposite=(2+obtuse)%3; a,b,q=points[i],points[j],points[opposite]
        def angle(x,y):
            denominator=np.linalg.norm(x)*np.linalg.norm(y)
            return 0.0 if denominator <= 1e-300 else math.acos(max(-1.0,min(1.0,float(x@y)/denominator)))
        cosine1=math.cos(angle(b-a,q-a)); cosine2=math.cos(angle(a-b,q-b)); dual[obtuse]=0.0
        if abs(cosine1)>1e-300:
            midpoint=(a+q)/2; foot=a+(b-a)/np.linalg.norm(b-a)*(np.linalg.norm(midpoint-a)/cosine1); dual[(obtuse+2)%3]=float(np.linalg.norm(midpoint-foot))
        if abs(cosine2)>1e-300:
            midpoint=(b+q)/2; foot=b+(a-b)/np.linalg.norm(a-b)*(np.linalg.norm(midpoint-b)/cosine2); dual[(obtuse+1)%3]=float(np.linalg.norm(midpoint-foot))
    return 0.5 * lengths[local_edge] * max(0.0, dual[local_edge])


def require_unique_keys(rows, keys):
    seen=set()
    for row in rows:
        key=tuple(row[name] for name in keys)
        if key in seen:
            raise ContractError(f"duplicate key {key}")
        seen.add(key)


def read_csv(path, expected_header=None):
    path=Path(path)
    with path.open(encoding="utf-8",newline="") as handle:
        reader=csv.DictReader(handle); header=tuple(reader.fieldnames or ()); rows=list(reader)
    if expected_header is not None and header != tuple(expected_header):
        raise ContractError(f"wrong CSV schema in {path.name}")
    for row in rows:
        for name,value in row.items():
            if name not in ("name","node_ids"):
                finite(value,f"{path.name}:{name}")
    return rows


def resolve_input(root, value):
    path=Path(value)
    return path if path.is_absolute() else root/path


def canonical_edges(triangles):
    return sorted({tuple(sorted((triangle[i],triangle[(i+1)%3]))) for triangle in triangles for i in range(3)})


def expected_edge_ids(triangles):
    result={}; next_id=0
    for triangle in triangles:
        for local in range(3):
            key=tuple(sorted((triangle[local],triangle[(local+1)%3])))
            if key not in result:
                result[key]=next_id; next_id+=1
    return result

def canonical_node_map(export_dir):
    rows=read_csv(export_dir/"nodes.csv",("id","x_um","y_um")); mapping={}; reverse={}
    for row in rows:
        source=int(row["id"]); xy=(finite(row["x_um"],"x_um"),finite(row["y_um"],"y_um"))
        candidates=[node for node,expected in NODES.items() if abs(xy[0]-expected[0])<COORDINATE_LIMIT_UM and abs(xy[1]-expected[1])<COORDINATE_LIMIT_UM]
        if len(candidates)!=1 or candidates[0] in reverse or source in mapping:
            raise ContractError("topology coordinate mapping is missing or duplicate")
        mapping[source]=candidates[0]; reverse[candidates[0]]=source
    if set(reverse)!=set(NODES):
        raise ContractError("topology coordinate mapping is incomplete")
    return mapping


def validate_topology(export_dir, state, topology_id):
    contract=state.get("topology_contract")
    if not isinstance(contract,dict) or (contract.get("nodes"),contract.get("triangles"),contract.get("edges"))!=(6,4,9):
        raise ContractError("wrong inline topology contract counts")
    contacts=contract.get("contact_edges")
    if contacts!={name:list(nodes) for name,nodes in CONTACTS.items()}:
        raise ContractError("wrong topology contact ownership")
    expected=TRIS[topology_id]
    contract_triangles=tuple(tuple(int(x) for x in triangle) for triangle in contract.get("triangle_connectivity",[]))
    if contract_triangles != expected:
        raise ContractError("inline topology must match exact approved CCW topology tuples/order")
    source_to_canonical=canonical_node_map(export_dir)
    elements=read_csv(export_dir/"elements.csv",("id","node0","node1","node2"))
    if len(elements) != 4 or {int(row["id"]) for row in elements} != set(range(4)):
        raise ContractError("wrong topology connectivity or counts")
    ordered=[]
    for row in sorted(elements, key=lambda value: int(value["id"])):
        triangle=tuple(source_to_canonical[int(row[f"node{i}"])] for i in range(3))
        if area2([NODES[node] for node in triangle])<=0:
            raise ContractError("reversed triangle orientation")
        ordered.append(triangle)
    if len(set(ordered)) != len(expected) or set(ordered) != set(expected):
        raise ContractError("elements must contain the exact approved CCW topology tuple set")
    contact_rows=read_csv(export_dir/"contacts.csv",("id","name","node_ids")); actual={}
    for row in contact_rows:
        actual[row["name"]]=tuple(sorted(source_to_canonical[int(x)] for x in row["node_ids"].split(";") if x!=""))
    if actual!={name:tuple(sorted(nodes)) for name,nodes in CONTACTS.items()}:
        raise ContractError("wrong topology contact ownership")
    doping_rows=read_csv(export_dir/"doping.csv",("node_id","donors_cm3","acceptors_cm3")); doping={}
    for row in doping_rows:
        node=source_to_canonical[int(row["node_id"])]
        if node in doping:
            raise ContractError("duplicate doping node")
        doping[node]=(finite(row["donors_cm3"],"donors_cm3"),finite(row["acceptors_cm3"],"acceptors_cm3"))
    if set(doping)!=set(NODES):
        raise ContractError("doping semantics are incomplete")
    for node in NODES:
        if doping[node]!=(DONORS_CM3[node],ACCEPTORS_CM3[node]):
            raise ContractError(f"wrong doping semantics at canonical node {node}")
    return source_to_canonical,doping


def select_fields(manifest, bias):
    if "bias_V" in manifest and finite(manifest["bias_V"],"field manifest bias")!=bias:
        raise ContractError("field manifest bias does not match requested/actual bias")
    fields=manifest.get("fields")
    if not isinstance(fields,list):
        raise ContractError("field manifest fields must be a list")
    selected={}
    for name,(components,unit) in FIELD_CONTRACT.items():
        matches=[field for field in fields if field.get("name")==name and field.get("region")==0 and field.get("components")==components]
        valid=[field for field in matches if field.get("unit")==unit and field.get("mapping_status")=="complete" and field.get("global_node_mapping")=="global_vertex_order"]
        if len(valid)!=1:
            if not matches:
                raise ContractError(f"missing required field {name} components={components}")
            raise ContractError(f"wrong unit or incomplete mapping for required field {name}")
        selected[name]=valid[0]
    return selected


def read_raw_fields(export_dir, source_to_canonical, field_manifest, bias):
    select_fields(field_manifest,bias); values={node:{} for node in NODES}
    for name,(components,_) in FIELD_CONTRACT.items():
        header=("node_id","component0") if components==1 else ("node_id","component0","component1")
        rows=read_csv(export_dir/"fields"/f"{name}_region0.csv",header)
        if len(rows)!=6:
            raise ContractError(f"partial required field {name}")
        seen=set()
        for row in rows:
            source=int(row["node_id"])
            if source not in source_to_canonical:
                raise ContractError(f"field {name} has unknown source node")
            node=source_to_canonical[source]
            if node in seen:
                raise ContractError(f"field {name} has duplicate node")
            seen.add(node)
            values[node][name]=tuple(finite(row[f"component{i}"],f"{name}:component{i}") for i in range(components))
    return values


def read_state(path):
    rows=read_csv(path,NODE_HEADER); result={}
    for row in rows:
        node=int(row["node_id"])+1
        if node not in NODES or node in result:
            raise ContractError("partial state or duplicate canonical node")
        result[node]={name:finite(row[name],name) for name in NODE_HEADER[1:]}
    if set(result)!=set(NODES):
        raise ContractError("partial state matrix")
    return result


def task4_paths(root,state,export_dir):
    aliases={
        "node":("vela_node_csv","vela_node_state.csv"),
        "edge":("vela_edge_csv","vela_edge_audit.csv"),
        "triangle":("vela_triangle_csv","vela_triangle_audit.csv"),
    }
    result={}
    for key,(state_key,fallback) in aliases.items():
        result[key]=resolve_input(root,state[state_key]) if state_key in state else export_dir/fallback
        if not result[key].is_file():
            raise ContractError(f"missing Task4 {key} CSV")
    return result


def load_audit_model(export_dir):
    config_path = export_dir / "audit.json"
    mesh_path = export_dir / "mesh.json"
    try:
        config = json.loads(config_path.read_text(encoding="utf-8"))
        mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContractError("missing or invalid immutable Task4 audit config/mesh") from error
    if finite(config.get("temperature_K"), "audit temperature") != TEMPERATURE_K:
        raise ContractError("audit config must use exactly 300 K")
    mobility = config.get("mobility")
    if not isinstance(mobility, dict) or mobility.get("model") != "constant":
        raise ContractError("audit mobility model must be constant")
    impact = config.get("impact_ionization")
    required_impact = {
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "current_approximation": "cell_reconstructed",
        "current_magnitude_mode": "edge_scalar_abs",
        "cell_reconstructed_midpoint_density": "gss_logistic",
        "quasi_fermi_gradient_discretization": "cell_gradient",
        "source_mapping_mode": "triangle_gss_gradqf_truncated",
        "source_volume_policy": "genius_truncated",
        "edge_source_partition": "symmetric",
    }
    if not isinstance(impact, dict) or any(impact.get(key) != value for key, value in required_impact.items()):
        raise ContractError("audit impact config does not match immutable Task4 model")
    for key, value in {
        "source_volume_factor": 0.0,
        "source_geometry_scale": 1.0,
        "quasi_fermi_carrier_truncation": 0.0,
        "minimum_field_V_m": 0.0,
    }.items():
        if finite(impact.get(key), f"audit impact {key}") != value:
            raise ContractError("audit impact config does not match immutable Task4 model")
    regions = mesh.get("regions")
    if not isinstance(regions, list) or len(regions) != 1 or regions[0].get("material") != "Si":
        raise ContractError("audit mesh must use the immutable default Silicon material")
    return SimpleNamespace(
        config_path=config_path,
        mesh_path=mesh_path,
        electron_mobility=SILICON_ELECTRON_MOBILITY_M2_PER_V_S,
        hole_mobility=SILICON_HOLE_MOBILITY_M2_PER_V_S,
        minimum_field=0.0,
    )


def validate_state_root(root):
    root=Path(root); manifest=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema")!=TASK3_SCHEMA or not manifest.get("outputs_complete"):
        raise ContractError("state root is not a complete vela.pn2d_minimal6_states.v1 artifact")
    states=manifest.get("states")
    if not isinstance(states,list) or len(states)!=6:
        raise ContractError("state root must contain six flat states")
    matrix={}; loaded=[]
    for state in states:
        topology_id=state.get("topology_id"); requested=finite(state.get("requested_bias_V"),"requested bias"); actual=finite(state.get("actual_bias_V"),"actual bias")
        if requested!=actual:
            raise ContractError("requested and actual bias must match exactly")
        if topology_id not in TOPOLOGIES or requested not in BIASES or state.get("status")!="passed":
            raise ContractError("state root contains unexpected topology, bias, or status")
        key=(topology_id,requested)
        if key in matrix:
            raise ContractError("duplicate topology/bias state")
        matrix[key]=state
        export_dir=resolve_input(root,state["export_dir"]); source_map,doping=validate_topology(export_dir,state,topology_id)
        field_path=resolve_input(root,state["field_manifest"]); field_manifest=json.loads(field_path.read_text(encoding="utf-8"))
        raw=read_raw_fields(export_dir,source_map,field_manifest,requested)
        state_csv=read_state(resolve_input(root,state["state_csv"])); paths=task4_paths(root,state,export_dir)
        model=load_audit_model(export_dir)
        loaded.append(SimpleNamespace(topology_id=topology_id,bias=requested,state=state,export_dir=export_dir,source_map=source_map,doping=doping,raw=raw,state_csv=state_csv,task4=paths,model=model))
    expected={(topology,bias) for topology in TOPOLOGIES for bias in BIASES}
    if set(matrix)!=expected:
        raise ContractError("states must contain exact biases 0, -12, -19 V for both topologies")
    return manifest,loaded


def raw_si_node(raw):
    return {
        "psi_V":raw["ElectrostaticPotential"][0],
        "phin_V":raw["eQuasiFermiPotential"][0],
        "phip_V":raw["hQuasiFermiPotential"][0],
        "n_m3":raw["eDensity"][0]*1e6,
        "p_m3":raw["hDensity"][0]*1e6,
        "donors_m3":raw["DonorConcentration"][0]*1e6,
        "acceptors_m3":raw["AcceptorConcentration"][0]*1e6,
        "electron_mobility_m2_per_V_s":raw["eMobility"][0]*1e-4,
        "hole_mobility_m2_per_V_s":raw["hMobility"][0]*1e-4,
        "electric_field_x_V_per_m":raw["ElectricField"][0]*100.0,
        "electric_field_y_V_per_m":raw["ElectricField"][1]*100.0,
        "electron_current_x_A_per_m2":raw["eCurrentDensity"][0]*1e4,
        "electron_current_y_A_per_m2":raw["eCurrentDensity"][1]*1e4,
        "hole_current_x_A_per_m2":raw["hCurrentDensity"][0]*1e4,
        "hole_current_y_A_per_m2":raw["hCurrentDensity"][1]*1e4,
        "electron_alpha_per_m":raw["eAlphaAvalanche"][0]*100.0,
        "hole_alpha_per_m":raw["hAlphaAvalanche"][0]*100.0,
    }

def build_report(state_root):
    root=Path(state_root); manifest,states=validate_state_root(root)
    node_rows=[]; edge_rows=[]; triangle_rows=[]; formula_errors=[]; state_errors=[]; external_errors=[]
    for item in states:
        topology=item.topology_id; bias=item.bias; triangles=TRIS[topology]; edge_ids=expected_edge_ids(triangles)
        sent={node:{**item.state_csv[node],**raw_si_node(item.raw[node])} for node in NODES}
        for node in NODES:
            donor_cm3, acceptor_cm3 = item.doping[node]
            gate(item.raw[node]["DonorConcentration"][0], donor_cm3, STATE_LIMIT,
                 "raw Task3 donor/doping.csv parity error")
            gate(item.raw[node]["AcceptorConcentration"][0], acceptor_cm3, STATE_LIMIT,
                 "raw Task3 acceptor/doping.csv parity error")
            for name in ("psi_V","phin_V","phip_V","n_m3","p_m3"):
                error=gate(sent[node][name],item.state_csv[node][name],STATE_LIMIT,"Task3 raw/state parity error"); state_errors.append(error)
        vela_nodes=read_csv(item.task4["node"],NODE_HEADER); vela_node_map={}
        for raw_node in vela_nodes:
            node=int(raw_node["node_id"])+1
            if node in vela_node_map: raise ContractError("duplicate Vela node key")
            vela_node_map[node]=raw_node
        if set(vela_node_map)!=set(NODES): raise ContractError("partial Vela node matrix")
        for node in NODES:
            raw=item.raw[node]; si=sent[node]; vela=vela_node_map[node]
            row={"topology_id":topology,"bias_V":bias,"node_id":node,"x_um":NODES[node][0],"y_um":NODES[node][1],
                 "raw_ElectrostaticPotential_V":raw["ElectrostaticPotential"][0],"raw_eQuasiFermiPotential_V":raw["eQuasiFermiPotential"][0],"raw_hQuasiFermiPotential_V":raw["hQuasiFermiPotential"][0],
                 "raw_eDensity_cm3":raw["eDensity"][0],"raw_hDensity_cm3":raw["hDensity"][0],"raw_DonorConcentration_cm3":raw["DonorConcentration"][0],"raw_AcceptorConcentration_cm3":raw["AcceptorConcentration"][0],
                 "raw_ElectricField_x_V_per_cm":raw["ElectricField"][0],"raw_ElectricField_y_V_per_cm":raw["ElectricField"][1],
                 "raw_eCurrentDensity_x_A_per_cm2":raw["eCurrentDensity"][0],"raw_eCurrentDensity_y_A_per_cm2":raw["eCurrentDensity"][1],
                 "raw_hCurrentDensity_x_A_per_cm2":raw["hCurrentDensity"][0],"raw_hCurrentDensity_y_A_per_cm2":raw["hCurrentDensity"][1],
                 "raw_eMobility_cm2_per_V_s":raw["eMobility"][0],"raw_hMobility_cm2_per_V_s":raw["hMobility"][0],
                 "raw_eAlphaAvalanche_per_cm":raw["eAlphaAvalanche"][0],"raw_hAlphaAvalanche_per_cm":raw["hAlphaAvalanche"][0]}
            for key,value in si.items(): row[f"sentaurus_{key}"]=value
            for name in NODE_HEADER[1:]:
                actual=finite(vela[name],f"Vela {name}"); expected=item.state_csv[node][name]
                error=gate(actual,expected,STATE_LIMIT,"state parity error"); state_errors.append(error)
                row[f"vela_{name}"]=actual; row[f"abs_error_{name}"]=absolute_error(actual,expected); row[f"hybrid_error_{name}"]=error
            node_rows.append(row)
        vela_triangles=read_csv(item.task4["triangle"],TRIANGLE_HEADER)
        if len(vela_triangles)!=4: raise ContractError("triangle row count mismatch")
        triangle_by_cell={int(row["cell_id"]):row for row in vela_triangles}
        if set(triangle_by_cell)!=set(range(4)): raise ContractError("duplicate or missing Vela triangle cell")
        edge_mobility={}
        for cell_id,triangle in enumerate(triangles):
            raw_triangle=triangle_by_cell[cell_id]
            raw_nodes=tuple(int(raw_triangle[f"node{i}"])+1 for i in range(3))
            if raw_nodes!=triangle: raise ContractError("Vela triangle identity/order mismatch")
            for local in range(3):
                a=triangle[local]; b=triangle[(local+1)%3]; key=tuple(sorted((a,b))); prefix=f"local_edge{local}_"
                if int(raw_triangle[prefix+"node0"])+1!=a or int(raw_triangle[prefix+"node1"])+1!=b or int(raw_triangle[prefix+"edge_id"])!=edge_ids[key]:
                    raise ContractError("Vela triangle local edge identity/order mismatch")
                mobility=(finite(raw_triangle[prefix+"electron_mobility_m2_per_V_s"],"electron mobility"),finite(raw_triangle[prefix+"hole_mobility_m2_per_V_s"],"hole mobility"))
                gate(mobility[0], item.model.electron_mobility, FORMULA_LIMIT,
                     "independent electron mobility error")
                gate(mobility[1], item.model.hole_mobility, FORMULA_LIMIT,
                     "independent hole mobility error")
                if key in edge_mobility and any(hybrid_error(x,y)>=FORMULA_LIMIT for x,y in zip(edge_mobility[key],mobility)):
                    raise ContractError("inconsistent production edge mobility")
                edge_mobility[key]=(item.model.electron_mobility, item.model.hole_mobility)
        vela_edges=read_csv(item.task4["edge"],EDGE_HEADER)
        if len(vela_edges)!=9: raise ContractError("edge row count mismatch")
        edge_map={}
        for raw_edge in vela_edges:
            a=int(raw_edge["node0"])+1; b=int(raw_edge["node1"])+1; key=(a,b)
            if a>=b or key in edge_map or key not in edge_ids or int(raw_edge["edge_id"])!=edge_ids[key]:
                raise ContractError("Vela edge_id or canonical endpoint identity/order mismatch")
            edge_map[key]=raw_edge
        if set(edge_map)!=set(canonical_edges(triangles)): raise ContractError("wrong topology edge matrix")
        state_edges={}
        for a,b in canonical_edges(triangles):
            s0=sent[a]; s1=sent[b]; point0=tuple(np.array(NODES[a])*1e-6); point1=tuple(np.array(NODES[b])*1e-6)
            dx=point1[0]-point0[0]; dy=point1[1]-point0[1]; length=math.hypot(dx,dy); raw_edge=edge_map[(a,b)]; mun,mup=edge_mobility[(a,b)]
            coefficient_e=mun*THERMAL_VOLTAGE_V/length; coefficient_h=mup*THERMAL_VOLTAGE_V/length
            electron_flux=sg_electron_variable_ni_flux(INTRINSIC_DENSITY_M3,INTRINSIC_DENSITY_M3,s0["psi_V"],s1["psi_V"],s0["phin_V"],s1["phin_V"],THERMAL_VOLTAGE_V,coefficient_e,True)
            hole_flux=sg_hole_variable_ni_flux(INTRINSIC_DENSITY_M3,INTRINSIC_DENSITY_M3,s0["psi_V"],s1["psi_V"],s0["phip_V"],s1["phip_V"],THERMAL_VOLTAGE_V,coefficient_h,True)
            vela_e=finite(raw_edge["electron_raw_signed_flux_per_m2_s"],"electron flux"); vela_h=finite(raw_edge["hole_raw_signed_flux_per_m2_s"],"hole flux")
            midpoint_e=gss_logistic_midpoint(s0["n_m3"],s1["n_m3"],s0["psi_V"],s1["psi_V"],THERMAL_VOLTAGE_V,"hole")
            midpoint_h=gss_logistic_midpoint(s0["p_m3"],s1["p_m3"],s0["psi_V"],s1["psi_V"],THERMAL_VOLTAGE_V,"electron")
            errors=[gate(finite(raw_edge["length_m"],"length"),length,FORMULA_LIMIT,"formula error"),gate(vela_e,electron_flux,FORMULA_LIMIT,"formula error"),gate(vela_h,hole_flux,FORMULA_LIMIT,"formula error"),gate(finite(raw_edge["electron_midpoint_density_m3"],"midpoint"),midpoint_e,FORMULA_LIMIT,"formula error"),gate(finite(raw_edge["hole_midpoint_density_m3"],"midpoint"),midpoint_h,FORMULA_LIMIT,"formula error")]
            formula_errors.extend(errors)
            sent_e_vector=((s0["electron_current_x_A_per_m2"]+s1["electron_current_x_A_per_m2"])/2,(s0["electron_current_y_A_per_m2"]+s1["electron_current_y_A_per_m2"])/2)
            sent_h_vector=((s0["hole_current_x_A_per_m2"]+s1["hole_current_x_A_per_m2"])/2,(s0["hole_current_y_A_per_m2"]+s1["hole_current_y_A_per_m2"])/2)
            sent_e_projection=canonical_projection(sent_e_vector,point0,point1); sent_h_projection=canonical_projection(sent_h_vector,point0,point1)
            adjacency=sum(a in triangle and b in triangle for triangle in triangles)
            edge_class="contact" if (a,b) in CONTACTS.values() else ("interior" if adjacency==2 else "boundary")
            eta_e=(s1["psi_V"]-s0["psi_V"])/THERMAL_VOLTAGE_V; eta_h=eta_e
            electron_formula_error=hybrid_error(vela_e,electron_flux); hole_formula_error=hybrid_error(vela_h,hole_flux)
            current_diag_e=hybrid_error(sent_e_projection,Q_C*vela_e); current_diag_h=hybrid_error(sent_h_projection,Q_C*vela_h); external_errors.extend((current_diag_e,current_diag_h))
            vela_edge_field_e=finite(raw_edge["electron_impact_field_V_per_m"],"impact field")
            vela_edge_field_h=finite(raw_edge["hole_impact_field_V_per_m"],"impact field")
            vela_edge_alpha_e=finite(raw_edge["electron_alpha_per_m"],"alpha")
            vela_edge_alpha_h=finite(raw_edge["hole_alpha_per_m"],"alpha")
            formula_errors.extend((
                gate(vela_edge_alpha_e, van_overstraeten_alpha(vela_edge_field_e,"electron"), FORMULA_LIMIT, "independent electron alpha error"),
                gate(vela_edge_alpha_h, van_overstraeten_alpha(vela_edge_field_h,"hole"), FORMULA_LIMIT, "independent hole alpha error"),
            ))
            row={"topology_id":topology,"bias_V":bias,"edge_id":edge_ids[(a,b)],"node0":a,"node1":b,"canonical_sign":"+1 node0-to-node1","edge_class":edge_class,"dx_m":dx,"dy_m":dy,"length_m":length,
                 "node0_psi_V":s0["psi_V"],"node1_psi_V":s1["psi_V"],"node0_phin_V":s0["phin_V"],"node1_phin_V":s1["phin_V"],"node0_phip_V":s0["phip_V"],"node1_phip_V":s1["phip_V"],"node0_n_m3":s0["n_m3"],"node1_n_m3":s1["n_m3"],"node0_p_m3":s0["p_m3"],"node1_p_m3":s1["p_m3"],
                 "delta_phin_over_h_V_per_m":(s1["phin_V"]-s0["phin_V"])/length,"delta_phip_over_h_V_per_m":(s1["phip_V"]-s0["phip_V"])/length,
                 "electron_eta":eta_e,"hole_eta":eta_h,"bernoulli_eta":bernoulli(eta_e),"bernoulli_minus_eta":bernoulli(-eta_e),"electron_mobility_m2_per_V_s":mun,"hole_mobility_m2_per_V_s":mup,
                 "electron_midpoint_density_m3":midpoint_e,"hole_midpoint_density_m3":midpoint_h,"python_electron_flux_per_m2_s":electron_flux,"python_hole_flux_per_m2_s":hole_flux,"pdf_electron_grad_qf_flux_per_m2_s":mun*midpoint_e*(s1["phin_V"]-s0["phin_V"])/length,"pdf_hole_grad_qf_flux_per_m2_s":-mup*midpoint_h*(s1["phip_V"]-s0["phip_V"])/length,
                 "vela_electron_flux_per_m2_s":vela_e,"vela_hole_flux_per_m2_s":vela_h,"vela_electron_current_A_per_m2":Q_C*vela_e,"vela_hole_current_A_per_m2":Q_C*vela_h,"electron_formula_abs_error":abs(vela_e-electron_flux),"hole_formula_abs_error":abs(vela_h-hole_flux),"electron_formula_hybrid_error":electron_formula_error,"hole_formula_hybrid_error":hole_formula_error,
                 "sentaurus_electron_projection_A_per_m2":sent_e_projection,"sentaurus_hole_projection_A_per_m2":sent_h_projection,"sentaurus_electron_magnitude_A_per_m2":math.hypot(*sent_e_vector),"sentaurus_hole_magnitude_A_per_m2":math.hypot(*sent_h_vector),"sentaurus_vs_vela_electron_current_diagnostic":current_diag_e,"sentaurus_vs_vela_hole_current_diagnostic":current_diag_h,
                 "vela_electron_impact_field_V_per_m":vela_edge_field_e,"vela_hole_impact_field_V_per_m":vela_edge_field_h,"vela_electron_alpha_per_m":vela_edge_alpha_e,"vela_hole_alpha_per_m":vela_edge_alpha_h,"edge_area_m2":finite(raw_edge["edge_area_m2"],"edge area")}
            row["vela_electron_edge_source_per_s"]=row["vela_electron_alpha_per_m"]*abs(vela_e)*row["edge_area_m2"]
            row["vela_hole_edge_source_per_s"]=row["vela_hole_alpha_per_m"]*abs(vela_h)*row["edge_area_m2"]
            edge_rows.append(row); state_edges[(a,b)]=row
        for cell_id,triangle in enumerate(triangles):
            raw_triangle=triangle_by_cell[cell_id]; points=[tuple(np.array(NODES[node])*1e-6) for node in triangle]; signed_area2=area2(points)
            gradients={name:triangle_gradient(points,[sent[node][name] for node in triangle]) for name in ("psi_V","phin_V","phip_V")}
            shape_gradients=[triangle_gradient(points,[1.0 if i==basis else 0.0 for i in range(3)]) for basis in range(3)]
            errors=[gate(finite(raw_triangle["signed_double_area_m2"],"area"),signed_area2,FORMULA_LIMIT,"formula error")]
            for name,prefix in (("psi_V","grad_psi"),("phin_V","grad_phin"),("phip_V","grad_phip")):
                errors.append(gate(finite(raw_triangle[f"{prefix}_x_V_per_m"],prefix),gradients[name][0],FORMULA_LIMIT,"formula error"))
                errors.append(gate(finite(raw_triangle[f"{prefix}_y_V_per_m"],prefix),gradients[name][1],FORMULA_LIMIT,"formula error"))
            matrices=[]; electron_projections=[]; hole_projections=[]
            for local in range(3):
                a=triangle[local]; b=triangle[(local+1)%3]; key=tuple(sorted((a,b))); edge=state_edges[key]
                point_a=tuple(np.array(NODES[key[0]])*1e-6); point_b=tuple(np.array(NODES[key[1]])*1e-6)
                matrices.append([(point_b[0]-point_a[0])/edge["length_m"],(point_b[1]-point_a[1])/edge["length_m"]])
                electron_projections.append(edge["vela_electron_flux_per_m2_s"]); hole_projections.append(edge["vela_hole_flux_per_m2_s"])
            electron_vector=np.linalg.lstsq(np.array(matrices),np.array(electron_projections),rcond=None)[0]
            hole_vector=np.linalg.lstsq(np.array(matrices),np.array(hole_projections),rcond=None)[0]
            vela_electron_source=0.0; vela_hole_source=0.0; python_electron_source=0.0; python_hole_source=0.0; sentaurus_source=0.0
            vela_electron_partition={node:0.0 for node in NODES}; vela_hole_partition={node:0.0 for node in NODES}
            python_electron_partition={node:0.0 for node in NODES}; python_hole_partition={node:0.0 for node in NODES}
            row={"topology_id":topology,"bias_V":bias,"cell_id":cell_id,"node0":triangle[0],"node1":triangle[1],"node2":triangle[2],"orientation":"CCW","signed_double_area_m2":signed_area2,"area_m2":signed_area2/2,
                 "shape_grad_N0_x_per_m":shape_gradients[0][0],"shape_grad_N0_y_per_m":shape_gradients[0][1],"shape_grad_N1_x_per_m":shape_gradients[1][0],"shape_grad_N1_y_per_m":shape_gradients[1][1],"shape_grad_N2_x_per_m":shape_gradients[2][0],"shape_grad_N2_y_per_m":shape_gradients[2][1],
                 "python_grad_psi_x_V_per_m":gradients["psi_V"][0],"python_grad_psi_y_V_per_m":gradients["psi_V"][1],"python_grad_phin_x_V_per_m":gradients["phin_V"][0],"python_grad_phin_y_V_per_m":gradients["phin_V"][1],"python_grad_phip_x_V_per_m":gradients["phip_V"][0],"python_grad_phip_y_V_per_m":gradients["phip_V"][1],
                 "vela_grad_psi_x_V_per_m":finite(raw_triangle["grad_psi_x_V_per_m"],"grad psi"),"vela_grad_psi_y_V_per_m":finite(raw_triangle["grad_psi_y_V_per_m"],"grad psi"),"vela_grad_phin_x_V_per_m":finite(raw_triangle["grad_phin_x_V_per_m"],"grad phin"),"vela_grad_phin_y_V_per_m":finite(raw_triangle["grad_phin_y_V_per_m"],"grad phin"),"vela_grad_phip_x_V_per_m":finite(raw_triangle["grad_phip_x_V_per_m"],"grad phip"),"vela_grad_phip_y_V_per_m":finite(raw_triangle["grad_phip_y_V_per_m"],"grad phip"),
                 "reconstructed_electron_particle_flux_x_per_m2_s":float(electron_vector[0]),"reconstructed_electron_particle_flux_y_per_m2_s":float(electron_vector[1]),"reconstructed_hole_particle_flux_x_per_m2_s":float(hole_vector[0]),"reconstructed_hole_particle_flux_y_per_m2_s":float(hole_vector[1]),"reconstructed_electron_current_x_A_per_m2":Q_C*float(electron_vector[0]),"reconstructed_electron_current_y_A_per_m2":Q_C*float(electron_vector[1]),"reconstructed_hole_current_x_A_per_m2":Q_C*float(hole_vector[0]),"reconstructed_hole_current_y_A_per_m2":Q_C*float(hole_vector[1])}
            for local in range(3):
                prefix=f"local_edge{local}_"; a=triangle[local]; b=triangle[(local+1)%3]; s0=sent[a]; s1=sent[b]
                partial=genius_truncated_partial_volume(points,local); electron_mid=gss_logistic_midpoint(s0["n_m3"],s1["n_m3"],s0["psi_V"],s1["psi_V"],THERMAL_VOLTAGE_V,"electron"); hole_mid=gss_logistic_midpoint(s0["p_m3"],s1["p_m3"],s0["psi_V"],s1["psi_V"],THERMAL_VOLTAGE_V,"hole")
                electron_field=math.hypot(*gradients["phin_V"]); hole_field=math.hypot(*gradients["phip_V"])
                cpp_mun=finite(raw_triangle[prefix+"electron_mobility_m2_per_V_s"],"mobility"); cpp_mup=finite(raw_triangle[prefix+"hole_mobility_m2_per_V_s"],"mobility")
                mun=item.model.electron_mobility; mup=item.model.hole_mobility
                edge_length=math.dist(tuple(np.array(NODES[a])*1e-6),tuple(np.array(NODES[b])*1e-6)); electron_edge_field=abs(s1["phin_V"]-s0["phin_V"])/edge_length; hole_edge_field=abs(s1["phip_V"]-s0["phip_V"])/edge_length; electron_proxy=mun*electron_mid*electron_edge_field; hole_proxy=mup*hole_mid*hole_edge_field
                electron_alpha=van_overstraeten_alpha(electron_field,"electron",minimum_field_V_per_m=item.model.minimum_field); hole_alpha=van_overstraeten_alpha(hole_field,"hole",minimum_field_V_per_m=item.model.minimum_field)
                cpp_electron_alpha=finite(raw_triangle[prefix+"electron_alpha_per_m"],"alpha"); cpp_hole_alpha=finite(raw_triangle[prefix+"hole_alpha_per_m"],"alpha")
                electron_local=electron_alpha*electron_proxy*partial; hole_local=hole_alpha*hole_proxy*partial
                cpp_partial=finite(raw_triangle[prefix+"truncated_partial_volume_m2"],"partial volume")
                cpp_electron_local=finite(raw_triangle[prefix+"electron_source_integral_per_m_s"],"source"); cpp_hole_local=finite(raw_triangle[prefix+"hole_source_integral_per_m_s"],"source")
                local_errors=[zero_normalized_gate(cpp_partial,partial,GEOMETRIC_ZERO_AREA_M2,FORMULA_LIMIT,"formula error"),gate(finite(raw_triangle[prefix+"electron_cell_qf_field_V_per_m"],"cell field"),electron_field,FORMULA_LIMIT,"formula error"),gate(finite(raw_triangle[prefix+"hole_cell_qf_field_V_per_m"],"cell field"),hole_field,FORMULA_LIMIT,"formula error"),gate(finite(raw_triangle[prefix+"electron_midpoint_density_m3"],"midpoint"),electron_mid,FORMULA_LIMIT,"formula error"),gate(finite(raw_triangle[prefix+"hole_midpoint_density_m3"],"midpoint"),hole_mid,FORMULA_LIMIT,"formula error"),gate(cpp_mun,mun,FORMULA_LIMIT,"independent electron mobility error"),gate(cpp_mup,mup,FORMULA_LIMIT,"independent hole mobility error"),gate(cpp_electron_alpha,electron_alpha,FORMULA_LIMIT,"independent electron alpha error"),gate(cpp_hole_alpha,hole_alpha,FORMULA_LIMIT,"independent hole alpha error"),gate(finite(raw_triangle[prefix+"electron_flux_proxy_per_m2_s"],"flux proxy"),electron_proxy,FORMULA_LIMIT,"formula error"),gate(finite(raw_triangle[prefix+"hole_flux_proxy_per_m2_s"],"flux proxy"),hole_proxy,FORMULA_LIMIT,"formula error"),geometric_source_gate(cpp_electron_local,electron_local,cpp_partial,partial,FORMULA_LIMIT,"formula error"),geometric_source_gate(cpp_hole_local,hole_local,cpp_partial,partial,FORMULA_LIMIT,"formula error")]
                errors.extend(local_errors)
                vela_electron_source+=cpp_electron_local; vela_hole_source+=cpp_hole_local; python_electron_source+=electron_local; python_hole_source+=hole_local
                vela_electron_partition[a]+=cpp_electron_local/2; vela_electron_partition[b]+=cpp_electron_local/2; vela_hole_partition[a]+=cpp_hole_local/2; vela_hole_partition[b]+=cpp_hole_local/2
                python_electron_partition[a]+=electron_local/2; python_electron_partition[b]+=electron_local/2; python_hole_partition[a]+=hole_local/2; python_hole_partition[b]+=hole_local/2
                sent_e=((s0["electron_current_x_A_per_m2"]+s1["electron_current_x_A_per_m2"])/2,(s0["electron_current_y_A_per_m2"]+s1["electron_current_y_A_per_m2"])/2); sent_h=((s0["hole_current_x_A_per_m2"]+s1["hole_current_x_A_per_m2"])/2,(s0["hole_current_y_A_per_m2"]+s1["hole_current_y_A_per_m2"])/2)
                sent_alpha_e=(s0["electron_alpha_per_m"]+s1["electron_alpha_per_m"])/2; sent_alpha_h=(s0["hole_alpha_per_m"]+s1["hole_alpha_per_m"])/2
                sentaurus_source+=(sent_alpha_e*math.hypot(*sent_e)/Q_C+sent_alpha_h*math.hypot(*sent_h)/Q_C)*partial
                for suffix in LOCAL_SUFFIXES: row[f"vela_{prefix}{suffix}"]=finite(raw_triangle[prefix+suffix],prefix+suffix)
                row[f"python_{prefix}truncated_partial_volume_m2"]=partial; row[f"python_{prefix}electron_source_integral_per_m_s"]=electron_local; row[f"python_{prefix}hole_source_integral_per_m_s"]=hole_local
            row["normalized_geometric_zero"] = any(abs(finite(raw_triangle[f"local_edge{i}_truncated_partial_volume_m2"],"partial volume")) <= GEOMETRIC_ZERO_AREA_M2 and abs(genius_truncated_partial_volume(points,i)) <= GEOMETRIC_ZERO_AREA_M2 for i in range(3))
            aggregate_error=gate(vela_electron_source+vela_hole_source,python_electron_source+python_hole_source,FORMULA_LIMIT,"raw C++ versus independent Python aggregate source error")
            errors.append(aggregate_error); formula_errors.extend(errors); row["max_formula_hybrid_error"]=max(errors)
            row["vela_electron_source_integral_per_m_s"]=vela_electron_source; row["vela_hole_source_integral_per_m_s"]=vela_hole_source; row["vela_total_source_integral_per_m_s"]=vela_electron_source+vela_hole_source
            row["python_electron_source_integral_per_m_s"]=python_electron_source; row["python_hole_source_integral_per_m_s"]=python_hole_source; row["python_total_source_integral_per_m_s"]=python_electron_source+python_hole_source
            row["vela_vs_python_total_source_hybrid_error"]=aggregate_error; row["sentaurus_total_source_integral_per_m_s"]=sentaurus_source
            diagnostic=hybrid_error(sentaurus_source,vela_electron_source+vela_hole_source); row["sentaurus_vs_vela_total_source_diagnostic"]=diagnostic; external_errors.append(diagnostic)
            for node in NODES:
                row[f"vela_electron_node{node}_source_partition_per_m_s"]=vela_electron_partition[node]; row[f"vela_hole_node{node}_source_partition_per_m_s"]=vela_hole_partition[node]
                row[f"python_electron_node{node}_source_partition_per_m_s"]=python_electron_partition[node]; row[f"python_hole_node{node}_source_partition_per_m_s"]=python_hole_partition[node]
            triangle_rows.append(row)
    require_unique_keys(node_rows,("topology_id","bias_V","node_id")); require_unique_keys(edge_rows,("topology_id","bias_V","node0","node1")); require_unique_keys(triangle_rows,("topology_id","bias_V","cell_id"))
    if (len(node_rows),len(edge_rows),len(triangle_rows))!=(36,54,24): raise ContractError("output row-count mismatch")
    orientation=build_orientation_summary(node_rows,edge_rows,triangle_rows)
    summary={"schema":SCHEMA,"status":"UNVERIFIED","scope":"fixed-state operator audit, not a BV curve","row_counts":{"node_state":36,"edge_audit":54,"triangle_audit":24},"gates":{"passed":False,"provenance_replay_validated":False,"state_parity_limit":STATE_LIMIT,"formula_limit":FORMULA_LIMIT,"sentaurus_vs_vela_current_source_threshold":None,"max_state_parity_hybrid_error":max(state_errors,default=0.0),"max_cpp_python_formula_hybrid_error":max(formula_errors,default=0.0),"max_external_diagnostic_hybrid_error":max(external_errors,default=0.0)},"orientation_sensitivity":orientation,"figure_count":14,"qa_notes":["Three biases are separate fixed-state diagnostic samples, not a BV curve.","External Sentaurus/Vela comparisons are diagnostic-only.","Zero pairs have no log ratio."]}
    return SimpleNamespace(node_rows=node_rows,edge_rows=edge_rows,triangle_rows=triangle_rows,summary=summary,manifest=manifest,fixture_root=root)


def numeric_columns(rows, excluded):
    return [name for name,value in rows[0].items() if name not in excluded and isinstance(value,(int,float)) and not isinstance(value,bool)]


def build_orientation_summary(node_rows,edge_rows,triangle_rows):
    result=[]
    def add(bias,quantity,sketch,mirror): result.append({"bias_V":bias,"quantity":quantity,**classify_orientation_pair(sketch,mirror)})
    for bias in BIASES:
        for node in NODES:
            pair={row["topology_id"]:row for row in node_rows if row["bias_V"]==bias and row["node_id"]==node}
            for column in numeric_columns(list(pair.values()),{"bias_V","node_id"}): add(bias,f"node_{node}_{column}",pair["sketch"][column],pair["mirror"][column])
        shared=set((row["node0"],row["node1"]) for row in edge_rows if row["topology_id"]=="sketch") & set((row["node0"],row["node1"]) for row in edge_rows if row["topology_id"]=="mirror")
        for edge in sorted(shared):
            pair={row["topology_id"]:row for row in edge_rows if row["bias_V"]==bias and (row["node0"],row["node1"])==edge}
            for column in numeric_columns(list(pair.values()),{"bias_V","edge_id","node0","node1"}): add(bias,f"edge_{edge[0]}_{edge[1]}_{column}",pair["sketch"][column],pair["mirror"][column])
        for kind,rows in (("edge",edge_rows),("triangle",triangle_rows)):
            sample=[row for row in rows if row["bias_V"]==bias]
            excluded={"bias_V","edge_id","node0","node1","cell_id","node2"}
            for column in numeric_columns(sample,excluded):
                sums={topology:sum(row[column] for row in sample if row["topology_id"]==topology) for topology in TOPOLOGIES}
                add(bias,f"integrated_{kind}_{column}",sums["sketch"],sums["mirror"])
    return result


def orientation_quantity_names(report):
    return {row["quantity"] for row in report.summary["orientation_sensitivity"]}

def write_csv(path,rows):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def command_status_from_provenance(manifest, replay_validated=False):
    replays=manifest.get("task4_provenance",{}).get("replays",[])
    replay_status="PASS" if replay_validated else "FAIL"
    return {"task4_replay":replay_status,"task4_replay_exit_codes":[replay.get("exit_code") for replay in replays],"report_generation":"PASS" if replay_validated else "BLOCKED"}

def plot_reports(report,out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure_dir=Path(out)/"figures"; figure_dir.mkdir(parents=True,exist_ok=True); blue="#2166ac"; orange="#e08214"
    figure,axes=plt.subplots(1,2,figsize=(13,4.8),constrained_layout=True)
    for axis,(topology,triangles) in zip(axes,TRIS.items()):
        for triangle in triangles:
            polygon=np.array([NODES[node] for node in triangle+(triangle[0],)]); axis.plot(polygon[:,0],polygon[:,1],color="#555555")
            center=polygon[:-1].mean(0); axis.text(center[0],center[1],"T"+"-".join(map(str,triangle)),fontsize=8,ha="center",bbox={"facecolor":"white","alpha":.8,"edgecolor":"none"})
        for a,b in canonical_edges(triangles):
            midpoint=(np.array(NODES[a])+np.array(NODES[b]))/2; axis.text(midpoint[0],midpoint[1]+.025,f"E{a}-{b}",fontsize=7,ha="center")
        for node,(x,y) in NODES.items(): axis.scatter(x,y,s=65,color=blue,edgecolor="black",zorder=3); axis.text(x,y-.055,f"N{node}",ha="center",fontsize=9,fontweight="bold")
        for node in (2,6): x,y=NODES[node]; axis.scatter(x,y,s=155,facecolor="none",edgecolor=orange,linewidth=2,zorder=4)
        axis.plot([0,0],[0,.5],color=orange,linewidth=4,label="Anode contact"); axis.plot([2,2],[0,.5],color=blue,linewidth=4,linestyle="--",label="Cathode contact")
        axis.set(title=f"{topology}: nodes, edges, triangles, contacts\nrings: compensated nodes 2 and 6",xlabel="x (um)",ylabel="y (um)",xlim=(-.12,2.12),ylim=(-.12,.64),aspect="equal"); axis.legend(fontsize=8,ncol=2)
    figure.suptitle("PN2D minimal6 canonical topologies — fixed-state operator audit")
    for extension in ("png","pdf"): figure.savefig(figure_dir/f"minimal6-topologies.{extension}",dpi=180)
    plt.close(figure)
    for bias,slug in ((0.0,"0v"),(-12.0,"minus12v"),(-19.0,"minus19v")):
        rows=[row for row in report.edge_rows if row["bias_V"]==bias]; x=np.arange(len(rows)); labels=[f"{row['topology_id'][0].upper()} {row['node0']}-{row['node1']}" for row in rows]
        figure,axis=plt.subplots(figsize=(14,6),constrained_layout=True)
        axis.scatter(x,[row["vela_electron_current_A_per_m2"] for row in rows],color=blue,marker="o",label="Vela electron current")
        axis.scatter(x,[row["vela_hole_current_A_per_m2"] for row in rows],color=orange,marker="s",label="Vela hole current")
        axis.axhline(0,color="#555",linewidth=.8); axis.set_yscale("symlog",linthresh=1e-12); axis.set_xticks(x,labels,rotation=55,ha="right",fontsize=8)
        axis.set(xlabel="topology and canonical edge (unconnected samples)",ylabel="signed reconstructed current (A m$^{-2}$)",title=f"PN2D minimal6 edge current at {bias:g} V\nfixed-state operator audit, not a BV curve")
        axis.legend(ncol=2); axis.grid(alpha=.2)
        for extension in ("png","pdf"): figure.savefig(figure_dir/f"minimal6-edge-current-audit-{slug}.{extension}",dpi=180)
        plt.close(figure)
        rows=[row for row in report.triangle_rows if row["bias_V"]==bias]; x=np.arange(len(rows)); labels=[f"{row['topology_id'][0].upper()} T{row['node0']}-{row['node1']}-{row['node2']}" for row in rows]
        figure,axis=plt.subplots(figsize=(12,5.5),constrained_layout=True)
        axis.scatter(x,[row["vela_electron_source_integral_per_m_s"] for row in rows],color=blue,marker="o",label="Vela electron source")
        axis.scatter(x,[row["vela_hole_source_integral_per_m_s"] for row in rows],color=orange,marker="s",label="Vela hole source")
        axis.axhline(0,color="#555",linewidth=.8); axis.ticklabel_format(axis="y",style="sci",scilimits=(0,0)); axis.set_xticks(x,labels,rotation=35,ha="right")
        axis.set(xlabel="topology and canonical triangle (unconnected samples)",ylabel="source integral (m$^{-1}$ s$^{-1}$)",title=f"PN2D minimal6 triangle source at {bias:g} V\nfixed-state operator audit, not a BV curve")
        axis.legend(); axis.grid(alpha=.2)
        for extension in ("png","pdf"): figure.savefig(figure_dir/f"minimal6-triangle-source-audit-{slug}.{extension}",dpi=180)
        plt.close(figure)


def write_report(report,out):
    failures=verify_task4_replay(report.fixture_root, REPO_ROOT/TASK4_PRODUCER)
    if failures:
        raise ContractError("Task4 provenance/replay validation failed: " + "; ".join(failures))
    report.summary["status"]="PASS"
    report.summary["gates"]["passed"]=True
    report.summary["gates"]["provenance_replay_validated"]=True
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    write_csv(out/"node_state.csv",report.node_rows); write_csv(out/"edge_audit.csv",report.edge_rows); write_csv(out/"triangle_audit.csv",report.triangle_rows)
    (out/"summary.json").write_text(json.dumps(report.summary,indent=2)+"\n",encoding="utf-8")
    all_files=sorted(path for path in report.fixture_root.rglob("*") if path.is_file())
    hashes={path.relative_to(report.fixture_root).as_posix():hashlib.sha256(path.read_bytes()).hexdigest() for path in all_files}
    status=command_status_from_provenance(report.manifest, replay_validated=True)
    output_manifest={"schema":SCHEMA,"source_state_schema":TASK3_SCHEMA,"topology_definitions":report.manifest["states"],"bias_states_V":list(BIASES),"model_configuration":{"workflow":"immutable_fixed_state","temperature_K":TEMPERATURE_K,"thermal_voltage_V":THERMAL_VOLTAGE_V,"includeNiGradientDrift":True,"solver_runs":False},"task4_provenance":report.manifest.get("task4_provenance"),"tool_versions":{"python":platform.python_version(),"numpy":np.__version__},"command_status":status,"input_sha256":hashes,"row_counts":report.summary["row_counts"],"gate_status":"PASS" if status["task4_replay"]=="PASS" and status["report_generation"]=="PASS" else "FAIL"}
    (out/"manifest.json").write_text(json.dumps(output_manifest,indent=2)+"\n",encoding="utf-8")
    gates=report.summary["gates"]
    markdown=f'''# PN2D Minimal6 Fixed-State Operator Audit

**Answer first:** the six-state C++ replay and all strict topology, state-parity, independent-formula, completeness, finiteness, and uniqueness gates pass. This is a fixed-state operator audit, not a BV curve, and supports no physical breakdown-voltage conclusion.

## Numeric results

- Rows: {len(report.node_rows)} node, {len(report.edge_rows)} edge, {len(report.triangle_rows)} triangle.
- Maximum state parity hybrid error: {gates["max_state_parity_hybrid_error"]:.17g} (gate `< {STATE_LIMIT:.1e}`).
- Maximum C++/Python formula hybrid error: {gates["max_cpp_python_formula_hybrid_error"]:.17g} (gate `< {FORMULA_LIMIT:.1e}`).
- Maximum diagnostic-only Sentaurus/Vela hybrid difference: {gates["max_external_diagnostic_hybrid_error"]:.17g} (no automatic threshold).
- Task 4 replay provenance status: {status["task4_replay"]} for six executable invocations.
- Orientation rows: {len(report.summary["orientation_sensitivity"])}.

## Scope and method

The input is an actual `vela.pn2d_minimal6_states.v1` root with flat sketch/mirror states at exactly 0 V, -12 V, and -19 V. Raw Task 3 units are preserved and converted to SI. Python independently recomputes variable-intrinsic-density SG continuity fluxes with `includeNiGradientDrift=true`, inverse-matrix gradients, GSS midpoint densities, and Genius-truncated partial volumes; columns labeled Vela come only from the recorded Task 4 C++ executable replay.

## Limitations

The committed state root is synthetic and exercises the real interfaces; it does not replace a live Sentaurus physics result. The three biases are discrete fixed-state samples, not a sweep or trend. Sentaurus-versus-Vela current/source comparisons are diagnostic-only.
'''
    (out/"summary.md").write_text(markdown,encoding="utf-8"); plot_reports(report,out)


def _manifest_path(value):
    return str(value).replace("\\", "/")


def _joined_manifest_path(base, name):
    base=_manifest_path(base).rstrip("/")
    return f"{base}/{name}"


def _expected_replay_arguments(state):
    export=state["export_dir"]
    node_out=_manifest_path(state.get("vela_node_csv", _joined_manifest_path(export,"vela_node_state.csv")))
    edge_out=_manifest_path(state.get("vela_edge_csv", _joined_manifest_path(export,"vela_edge_audit.csv")))
    triangle_out=_manifest_path(state.get("vela_triangle_csv", _joined_manifest_path(export,"vela_triangle_audit.csv")))
    return [
        "--mesh", _joined_manifest_path(export,"mesh.json"),
        "--doping", _joined_manifest_path(export,"doping.csv"),
        "--state", _manifest_path(state["state_csv"]),
        "--config", _joined_manifest_path(export,"audit.json"),
        "--node-out", node_out,
        "--edge-out", edge_out,
        "--triangle-out", triangle_out,
    ]


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_task4_replay(root,executable):
    root=Path(root); executable=Path(executable); failures=[]
    try:
        manifest=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    except (OSError,json.JSONDecodeError) as error:
        return [f"invalid provenance manifest: {error}"]
    provenance=manifest.get("task4_provenance")
    if not isinstance(provenance,dict):
        return ["missing Task4 provenance"]
    if provenance.get("producer")!=TASK4_PRODUCER:
        failures.append("unexpected producer identity")
    if provenance.get("task4_source_commit")!=TASK4_SOURCE_COMMIT:
        failures.append("unexpected Task4 source commit")
    if not executable.is_file():
        failures.append("producer executable is missing")
    elif _sha256(executable)!=provenance.get("producer_sha256"):
        failures.append("producer hash mismatch")

    states=manifest.get("states")
    replays=provenance.get("replays")
    expected_matrix={(topology,bias) for topology in TOPOLOGIES for bias in BIASES}
    if not isinstance(states,list) or len(states)!=6:
        failures.append("state manifest does not contain six exact replay identities")
        states=[]
    state_map={}
    for state in states:
        try:
            identity=(state.get("topology_id"),finite(state.get("actual_bias_V"),"provenance state bias"))
        except ContractError:
            failures.append("state manifest has invalid replay identity")
            continue
        if identity in state_map:
            failures.append("state manifest has duplicate replay identity")
        state_map[identity]=state
    if set(state_map)!=expected_matrix:
        failures.append("state manifest does not contain six exact replay identities")
    if not isinstance(replays,list) or len(replays)!=6:
        failures.append("provenance must contain six exact replay identities")
        replays=[]
    replay_map={}
    for replay in replays:
        if not isinstance(replay,dict):
            failures.append("invalid replay record")
            continue
        try:
            identity=(replay.get("topology_id"),finite(replay.get("bias_V"),"replay bias"))
        except ContractError:
            failures.append("invalid replay identity")
            continue
        if identity in replay_map:
            failures.append("provenance must contain six exact replay identities")
        replay_map[identity]=replay
    if set(replay_map)!=expected_matrix:
        failures.append("provenance must contain six exact replay identities")
    if failures:
        return sorted(set(failures))

    options=("--mesh","--doping","--state","--config","--node-out","--edge-out","--triangle-out")
    input_options=options[:4]; output_options=options[4:]
    environment=os.environ.copy(); tool_dirs=[r"D:\msys64\ucrt64\bin",r"D:\msys64\usr\bin"]
    environment["PATH"]=os.pathsep.join(tool_dirs+[environment.get("PATH","")])
    for identity in sorted(expected_matrix):
        state=state_map[identity]; replay=replay_map[identity]; expected_arguments=_expected_replay_arguments(state)
        arguments=replay.get("arguments")
        label=f"{identity[0]} {identity[1]:g}V"
        if replay.get("producer")!=TASK4_PRODUCER:
            failures.append(f"replay producer identity mismatch: {label}")
        if arguments!=expected_arguments:
            failures.append(f"replay arguments mismatch: {label}")
            continue
        expected_command=" ".join([TASK4_PRODUCER,*expected_arguments])
        if replay.get("command")!=expected_command:
            failures.append(f"replay command mismatch: {label}")
        if replay.get("exit_code")!=0:
            failures.append(f"recorded replay exit code is not zero: {label}")
        values=dict(zip(arguments[0::2],arguments[1::2]))
        if tuple(arguments[0::2])!=options or len(values)!=len(options):
            failures.append(f"replay option identity/order mismatch: {label}")
            continue
        expected_inputs={values[option] for option in input_options}
        expected_outputs={values[option] for option in output_options}
        input_hashes=replay.get("input_sha256")
        output_hashes=replay.get("output_sha256")
        if not isinstance(input_hashes,dict) or set(input_hashes)!=expected_inputs:
            failures.append(f"recorded input hashes are incomplete or extra: {label}")
            continue
        if not isinstance(output_hashes,dict) or set(output_hashes)!=expected_outputs:
            failures.append(f"committed output hashes are incomplete or extra: {label}")
            continue
        for relative,expected_hash in input_hashes.items():
            path=Path(relative); path=path if path.is_absolute() else root/path
            if not path.is_file() or _sha256(path)!=expected_hash:
                failures.append(f"input hash mismatch: {relative}")
        for relative,expected_hash in output_hashes.items():
            path=Path(relative); path=path if path.is_absolute() else root/path
            if not path.is_file() or _sha256(path)!=expected_hash:
                failures.append(f"committed output hash mismatch: {relative}")
        with tempfile.TemporaryDirectory() as directory:
            directory=Path(directory); replay_arguments=list(arguments); fresh={
                "--node-out":directory/"node.csv",
                "--edge-out":directory/"edge.csv",
                "--triangle-out":directory/"triangle.csv",
            }
            for index in range(0,len(replay_arguments),2):
                option=replay_arguments[index]; value=Path(replay_arguments[index+1])
                if option in fresh:
                    replay_arguments[index+1]=str(fresh[option])
                else:
                    replay_arguments[index+1]=str(value if value.is_absolute() else root/value)
            completed=subprocess.run([str(executable),*replay_arguments],cwd=REPO_ROOT,capture_output=True,text=True,env=environment)
            if completed.returncode!=0:
                failures.append(f"fresh replay exit code is not zero: {label}")
                continue
            for option,path in fresh.items():
                relative=values[option]
                if not path.is_file() or _sha256(path)!=output_hashes[relative]:
                    failures.append(f"fresh replay output hash mismatch: {relative}")
    return failures


def main():
    parser=argparse.ArgumentParser(); group=parser.add_mutually_exclusive_group(required=True); group.add_argument("--state-root",type=Path); group.add_argument("--fixture",type=Path,help="alias for --state-root"); parser.add_argument("--out-dir",type=Path,required=True); args=parser.parse_args()
    report=build_report(args.state_root or args.fixture); write_report(report,args.out_dir)
    print(f"PASS {SCHEMA}: node={len(report.node_rows)} edge={len(report.edge_rows)} triangle={len(report.triangle_rows)} figures=14")


if __name__=="__main__":
    main()