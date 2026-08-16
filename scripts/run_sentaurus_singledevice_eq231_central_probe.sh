#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 SUPPORT_DIR PROBE_ARCHIVE PROBE_DIR" >&2
  exit 2
fi

support_dir=$1
archive=$2
probe_dir=$3

mkdir -p "$probe_dir"
tar -xzf "$archive" -C "$probe_dir"

for perturbed in "$probe_dir"/lin_eq231_perturbed_n*_des.tdr; do
  label=$(basename "$perturbed" _des.tdr)
  cp "$perturbed" "$support_dir/lin_eq231_perturbed_des.tdr"
  (
    cd "$support_dir"
    sdevice sentaurus_singledevice_eq231_perturbation_probe.cmd \
      > "$probe_dir/${label}.stdout" 2>&1
  )
  mv "$support_dir/eq231_perturbation_newton_0.tdr" \
    "$probe_dir/${label}_newton0.tdr"
  mv "$support_dir/eq231_perturbation_newton_1.tdr" \
    "$probe_dir/${label}_newton1.tdr"
  echo "completed $label"
done
