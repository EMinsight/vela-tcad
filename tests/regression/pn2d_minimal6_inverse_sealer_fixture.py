import csv
import hashlib
import json
from pathlib import Path


REQUIRED = {
    "ElectrostaticPotential": (1, "V"), "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"), "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"), "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"), "hCurrentDensity": (2, "A*cm^-2"),
    "eAlphaAvalanche": (1, "cm^-1"), "hAlphaAvalanche": (1, "cm^-1"),
    "ImpactIonization": (1, "cm^-3*s^-1"),
}
SUPPLEMENTAL = {
    "eMobility": (1, "cm^2*V^-1*s^-1"),
    "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eVelocity": (1, "cm*s^-1"), "hVelocity": (1, "cm*s^-1"),
}
COORDINATES = {
    0: (0.0, 0.5), 1: (1.0, 0.5), 2: (2.0, 0.5),
    3: (2.0, 0.0), 4: (0.0, 0.0), 5: (1.0, 0.0),
}
TRIANGLES = {
    "sketch": ((0, 4, 1), (4, 5, 1), (1, 5, 3), (1, 3, 2)),
    "mirror": ((0, 4, 5), (0, 5, 1), (1, 5, 2), (5, 3, 2)),
}
SOURCE_IDS = {0: 5, 1: 2, 2: 4, 3: 1, 4: 3, 5: 0}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def mutate_json(path: Path, mutate) -> None:
    value = json.loads(path.read_text(encoding="utf-8"))
    mutate(value)
    write_json(path, value)


def mesh(topology: str) -> dict:
    return {
        "coordinate_unit": "um",
        "nodes": [{"id": node, "x": xy[0], "y": xy[1]}
                  for node, xy in COORDINATES.items()],
        "triangles": [{"id": index, "node_ids": list(nodes), "region_id": 0}
                      for index, nodes in enumerate(TRIANGLES[topology])],
        "regions": [{"id": 0, "name": "R.Si", "material": "Si",
                     "cell_ids": [0, 1, 2, 3]}],
        "contacts": [
            {"id": 0, "name": "Anode", "node_ids": [0, 4], "region_id": 0},
            {"id": 1, "name": "Cathode", "node_ids": [2, 3], "region_id": 0},
        ],
    }


def raw_tdr(topology: str, bias: float, fields: dict[str, tuple[int, str]]) -> dict:
    values = {}
    for field_index, (name, (components, unit)) in enumerate(fields.items()):
        by_source = {}
        for canonical, source in SOURCE_IDS.items():
            by_source[str(source)] = [
                float(1000 * (field_index + 1) + 10 * canonical + component + abs(bias))
                for component in range(components)
            ]
        values[name] = {"components": components, "unit": unit, "values": by_source,
                        "mapping_status": "complete"}
    return {
        "topology": topology, "bias_V": bias,
        "nodes": [{"id": source, "x_um": COORDINATES[canonical][0],
                   "y_um": COORDINATES[canonical][1]}
                  for canonical, source in SOURCE_IDS.items()],
        "triangles": [[SOURCE_IDS[node] for node in nodes] for nodes in TRIANGLES[topology]],
        "contacts": {"Anode": [SOURCE_IDS[0], SOURCE_IDS[4]],
                     "Cathode": [SOURCE_IDS[2], SOURCE_IDS[3]]},
        "fields": values,
    }


def fake_importer(tdr: Path, output: Path) -> None:
    raw = json.loads(tdr.read_text(encoding="utf-8"))
    output.mkdir(parents=True, exist_ok=False)
    with (output / "nodes.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("id", "x_um", "y_um"))
        for row in raw["nodes"]:
            writer.writerow((row["id"], row["x_um"], row["y_um"]))
    with (output / "elements.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("id", "node0", "node1", "node2", "region", "material"))
        for index, nodes in enumerate(raw["triangles"]):
            writer.writerow((index, *nodes, "R.Si", "Si"))
    with (output / "contacts.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("name", "node_ids", "region"))
        for name in ("Cathode", "Anode"):
            writer.writerow((name, ";".join(str(node) for node in raw["contacts"][name]), "R.Si"))
    fields_dir = output / "fields"
    fields_dir.mkdir()
    manifest_fields = []
    for name, field in raw["fields"].items():
        manifest_fields.append({
            "name": name, "components": field["components"], "unit": field["unit"],
            "region": 0, "mapping_status": field.get("mapping_status", "complete"),
            "values": 6, "global_vertex_count": 6, "global_node_mapping": "global_vertex_order",
        })
        with (fields_dir / f"{name}_region0.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(("node_id", *(f"component{i}" for i in range(field["components"]))))
            for source in sorted(int(value) for value in field["values"]):
                writer.writerow((source, *field["values"][str(source)]))
    for region, region_name, value in (
            (1, "Cathode", 0.0), (2, "Anode", raw["bias_V"])):
        with (fields_dir / f"ContactExternalVoltage_region{region}.csv").open(
                "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(("node_id", "component0"))
            writer.writerow((raw["contacts"][region_name][0], value))
        manifest_fields.append({"name": "ContactExternalVoltage", "components": 1,
                                "unit": "V", "region": region, "mapping_status": "scalar",
                                "values": 1, "global_vertex_count": 6,
                                "global_node_mapping": "contact_scalar",
                                "region_name": region_name})
    write_json(output / "field_manifest.json", {"fields": manifest_fields})


def make_fixture(base: Path) -> tuple[Path, Path, Path, Path, Path]:
    vela = base / "vela-sweep"
    sentaurus = base / "sentaurus-sweep"
    supplemental = base / "supplemental"
    runner, importer = base / "vela_example_runner.exe", base / "sentaurus_import.exe"
    runner.write_bytes(b"runner")
    importer.write_bytes(b"importer")
    vela_segments, vela_accepted = [], []
    sent_segments, sent_accepted = [], []
    topology_hashes = {}
    for topology in ("sketch", "mirror"):
        mesh_path = vela / "inputs" / topology / "mesh.json"
        doping_path = vela / "inputs" / topology / "doping.csv"
        write_json(mesh_path, mesh(topology))
        doping_path.parent.mkdir(parents=True, exist_ok=True)
        doping_path.write_text("node_id,donors_cm3,acceptors_cm3\n" +
                               "\n".join(f"{node},1e15,0" for node in range(6)) + "\n",
                               encoding="utf-8")
        topology_hashes[topology] = {"mesh.json": sha256(mesh_path),
                                     "doping.csv": sha256(doping_path)}
        for magnitude in range(1, 21):
            bias, index = -float(magnitude), magnitude - 1
            token = f"m{magnitude}p000000"
            deck = vela / "vela" / topology / "decks" / f"segment_{index:02d}.json"
            write_json(deck, {"topology": topology, "stop": bias})
            state = vela / "vela" / topology / "states" / f"segment_{index:02d}_bias_{token}.csv"
            state.parent.mkdir(parents=True, exist_ok=True)
            with state.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(("node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"))
                for node in range(6):
                    writer.writerow((node, bias + node / 10, bias + node / 20,
                                     bias - node / 20, 1.0e18 + node, 2.0e18 + node))
            vela_segments.append({"solver": "vela", "topology": topology,
                                  "start_bias_V": bias + 1.0, "target_bias_V": bias,
                                  "status": "accepted", "deck": deck.relative_to(vela).as_posix(),
                                  "deck_sha256": sha256(deck)})
            vela_accepted.append({"solver": "vela", "topology": topology,
                                  "target_bias_V": bias, "actual_bias_V": bias,
                                  "status": "accepted", "state_path": state.relative_to(vela).as_posix(),
                                  "state_sha256": sha256(state)})
        for magnitude in range(21):
            bias = -float(magnitude)
            token = "0p000000" if magnitude == 0 else f"m{magnitude}p000000"
            deck = sentaurus / "sentaurus" / topology / "decks" / f"{topology}_{token}.cmd"
            deck.parent.mkdir(parents=True, exist_ok=True)
            deck.write_text(f"bias={bias:.17g}\n", encoding="utf-8")
            checkpoint = sentaurus / "sentaurus" / topology / "checkpoints" / f"{topology}_{token}.tdr"
            write_json(checkpoint, raw_tdr(topology, bias, REQUIRED))
            sent_segments.append({"solver": "sentaurus", "topology": topology,
                                  "start_bias_V": bias + 1.0, "target_bias_V": bias,
                                  "status": "accepted", "deck": deck.relative_to(sentaurus).as_posix(),
                                  "deck_sha256": sha256(deck)})
            old_export = sentaurus / "sentaurus" / topology / "exports" / f"{topology}_{token}"
            write_json(old_export / "field_manifest.json", {"old": True})
            sent_accepted.append({"solver": "sentaurus", "topology": topology,
                                  "target_bias_V": bias, "actual_bias_V": bias,
                                  "status": "accepted", "state_path": checkpoint.relative_to(sentaurus).as_posix(),
                                  "state_sha256": sha256(checkpoint),
                                  "export_dir": str(old_export.resolve()),
                                  "export_field_manifest_sha256": sha256(old_export / "field_manifest.json")})
    common = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1",
              "targets_V": [float(-value) for value in range(21)],
              "interpolation": "forbidden", "failed_transition": None,
              "failed_transitions": []}
    write_json(vela / "sweep_manifest.json",
               {**common, "topology_input_sha256": topology_hashes,
                "segments": vela_segments, "sentaurus_segments": [],
                "accepted_checkpoints": vela_accepted})
    write_json(sentaurus / "sweep_manifest.json",
               {**common, "segments": [], "sentaurus_segments": sent_segments,
                "accepted_checkpoints": sent_accepted})
    supplemental_states, expected_matrix = [], []
    for topology in ("sketch", "mirror"):
        for magnitude in range(1, 21):
            bias, tag = -float(magnitude), f"m{magnitude}V"
            state_root = supplemental / "states" / topology / tag
            tdr = state_root / "artifacts" / f"pn2d_minimal6_state_{tag}.tdr"
            write_json(tdr, raw_tdr(topology, bias, SUPPLEMENTAL))
            export = state_root / "export"
            fake_importer(tdr, export)
            ledger = {path.relative_to(export).as_posix(): sha256(path)
                      for path in export.rglob("*") if path.is_file()}
            supplemental_states.append({
                "topology_id": topology, "requested_bias_V": bias, "actual_bias_V": bias,
                "bias_tag": tag, "status": "passed", "sentaurus_version": "O-2018.06-SP2",
                "final_tdr_name": tdr.name, "member_sha256": ledger,
                "export_dir": str(export.resolve()), "artifacts_dir": str(tdr.parent.resolve()),
            })
            expected_matrix.append([topology, bias])
    write_json(supplemental / "manifest.json", {
        "schema": "vela.pn2d_minimal6_states.v1", "run_id": "fixture",
        "bias_tolerance_V": 1.0e-12, "expected_matrix": expected_matrix,
        "sentaurus_version": "O-2018.06-SP2", "outputs_complete": True,
        "importer": str(importer.resolve()),
        "states": supplemental_states,
        "manifest_path": str((supplemental / "manifest.json").resolve()),
    })
    return vela, sentaurus, supplemental, importer, runner
