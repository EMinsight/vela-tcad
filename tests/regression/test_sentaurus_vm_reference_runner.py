#!/usr/bin/env python3
"""Regression coverage for the opt-in Sentaurus VM reference runner."""

from __future__ import annotations

import json
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from scripts import run_bvmethods_nmos_low_bias_sentaurus_vm as low_bias_vm
from scripts import run_bvmethods_nmos_mean_ionization_controls_vm as mean_controls_vm
from scripts import run_bvmethods_nmos_multibias_sentaurus_vm as multibias_vm


REPO = Path(__file__).resolve().parents[2]
RUNNER = REPO / "scripts" / "run_sentaurus_vm_reference.py"
HOST_OVERRIDE_MODULES = (low_bias_vm, multibias_vm, mean_controls_vm)


class SentaurusVmHostConfigurationTest(unittest.TestCase):
    def test_specialized_runners_default_to_ssh_config(self) -> None:
        for module in HOST_OVERRIDE_MODULES:
            with self.subTest(module=module.__name__):
                self.assertEqual(module.ssh_host_options(None), [])
                self.assertEqual(module.ssh_host_options(""), [])

    def test_specialized_runners_keep_explicit_hostname_override(self) -> None:
        for module in HOST_OVERRIDE_MODULES:
            with self.subTest(module=module.__name__):
                self.assertEqual(
                    module.ssh_host_options("sentaurus.example"),
                    ["-o", "HostName=sentaurus.example"],
                )

    def test_specialized_runners_do_not_embed_ipv4_addresses(self) -> None:
        ipv4 = re.compile(r"(?<![0-9])(?:[0-9]{1,3}[.]){3}[0-9]{1,3}(?![0-9])")
        for module in HOST_OVERRIDE_MODULES:
            path = Path(module.__file__)
            with self.subTest(path=path.name):
                self.assertEqual(ipv4.findall(path.read_text(encoding="utf-8")), [])


class SentaurusVmReferenceRunnerTest(unittest.TestCase):
    def test_dry_run_writes_manifest_without_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vm_dry_") as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            for name in [
                "pn2d_sde.cmd",
                "pn2d_bv_sdevice.cmd",
                "models.par",
            ]:
                (source / name).write_text(f"{name}\n")
            out = root / "runs"

            subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--ssh-target", "sentaurus",
                "--source-dir", str(source),
                "--local-output-dir", str(out),
                "--remote-root", "~/sentaurus_runs/vela_oracle",
                "--run-id", "pn2d_bv_vm_dry_run",
                "--stages", "bv",
                "--dry-run",
            ], check=True)

            manifest = json.loads(
                (out / "pn2d_bv_vm_dry_run" / "sentaurus_vm_run_manifest.json").read_text()
            )
            self.assertEqual(manifest["ssh_target"], "sentaurus")
            self.assertIsNone(manifest["sentaurus_version"])
            self.assertEqual(
                manifest["remote_source_dir"],
                "~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source",
            )
            self.assertEqual(manifest["stages"], ["bv"])
            self.assertEqual(manifest["commands"], [
                "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source && sde -e -l pn2d_sde.cmd",
                "cd ~/sentaurus_runs/vela_oracle/pn2d_bv_vm_dry_run/source && sdevice pn2d_bv_sdevice.cmd > run_pn2d_bv.out 2>&1",
            ])

    def test_2022_0v_smoke_plan_records_release_and_uses_versioned_output(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus2022_0v_") as tmp:
            root = Path(tmp)
            source = root / "pn2d_sentaurus2018" / "source"
            source.mkdir(parents=True)
            for name in [
                "pn2d_sde.cmd",
                "pn2d_0v_sdevice.cmd",
                "models.par",
            ]:
                (source / name).write_text(f"{name}\n")
            out = root / "pn2d_sentaurus2022" / "sentaurus_vm_runs"

            subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--ssh-target", "sentaurus",
                "--sentaurus-version", "T-2022.03-SP2",
                "--source-dir", str(source),
                "--local-output-dir", str(out),
                "--remote-root", "~/sentaurus_runs/vela_oracle_2022",
                "--run-id", "sentaurus2022_license_smoke",
                "--stages", "0v",
                "--dry-run",
            ], check=True)

            manifest = json.loads(
                (out / "sentaurus2022_license_smoke" / "sentaurus_vm_run_manifest.json").read_text()
            )
            self.assertEqual(manifest["sentaurus_version"], "T-2022.03-SP2")
            self.assertEqual(manifest["stages"], ["0v"])
            self.assertEqual(manifest["required_files"], [
                "pn2d_sde.cmd",
                "models.par",
                "pn2d_0v_sdevice.cmd",
            ])
            self.assertEqual(manifest["commands"], [
                "cd ~/sentaurus_runs/vela_oracle_2022/sentaurus2022_license_smoke/source && sde -e -l pn2d_sde.cmd",
                "cd ~/sentaurus_runs/vela_oracle_2022/sentaurus2022_license_smoke/source && sdevice pn2d_0v_sdevice.cmd > run_pn2d_0v.out 2>&1",
            ])

    def test_missing_required_deck_fails_before_ssh(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_sentaurus_vm_missing_") as tmp:
            root = Path(tmp)
            source = root / "source"
            source.mkdir()
            (source / "pn2d_sde.cmd").write_text("mesh\n")
            completed = subprocess.run([
                sys.executable,
                str(RUNNER),
                "pn2d",
                "--source-dir", str(source),
                "--local-output-dir", str(root / "runs"),
                "--run-id", "missing_bv",
                "--stages", "bv",
                "--dry-run",
            ], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            self.assertNotEqual(completed.returncode, 0)
            self.assertIn("missing required source file", completed.stderr)


if __name__ == "__main__":
    unittest.main()
