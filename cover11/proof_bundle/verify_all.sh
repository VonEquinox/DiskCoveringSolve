#!/usr/bin/env bash
set -euo pipefail
MODE="${1:-quick}"
HERE="$(cd "$(dirname "$0")" && pwd)"
# Existing scripts use /mnt/data.  When launched elsewhere, expose package
# files there via symlinks.
if [[ "$HERE" != "/mnt/data" ]]; then
  mkdir -p /mnt/data
  for f in "$HERE"/*; do ln -sfn "$f" "/mnt/data/$(basename "$f")"; done
fi
cd /mnt/data
run(){ echo "=== $*"; "$@"; }
run python cover11_fixed_trig.py
run python cover11_reduction_constants_verify.py
run python cover11_enumeration_verify.py
run python cover11_metric_Tplus.py verify cover11_metric_Tplus.json
run python cover11_topology_linkage_verify.py
run python cover11_full_kkt_certificate.py verify
run python cover11_local_exact_audit2.py
run python cover11_local_theorem_verify_v2.py
run python cover11_candidate_structure_verify.py
run python cover11_exact_upper_candidate_v2.py
run python cover11_wheel_row_validity_verify.py
if [[ "$MODE" == "full" ]]; then
  ranges=("0 21228" "21228 42456" "42456 63684" "63684 84912" "84912 106136")
  for i in 0 1 2 3 4; do
    read -r a b <<<"${ranges[$i]}"
    run python cover11_candidate_base_full_verify.py --start "$a" --end "$b" --out "candbaseverify_${i}.json"
  done
  ranges=("0 1060" "1060 2120" "2120 3180" "3180 4240" "4240 5302")
  for i in 0 1 2 3 4; do
    read -r a b <<<"${ranges[$i]}"
    run python cover11_candidate_refined_full_verify.py --start "$a" --end "$b" --out "candrefverify_${i}.json"
  done
  ranges=("0 607" "607 1214" "1214 1821" "1821 2428" "2428 3035" "3035 3642" "3642 4249" "4249 4853")
  for i in 0 1 2 3 4 5 6 7; do
    read -r a b <<<"${ranges[$i]}"
    run python cover11_wheel_verify_chunk.py --start "$a" --end "$b" --out "wheelchunk_${i}.json"
  done
  ranges=("0 11" "11 22" "22 33" "33 43" "43 52" "52 53")
  for i in 0 1 2 3 4 5; do
    read -r a b <<<"${ranges[$i]}"
    run python cover11_residual_all_verify.py --start "$a" --end "$b" --out "resverify_${i}.json"
  done
fi
run python cover11_aggregate_results.py
printf '\nALL REQUESTED COVER11 VERIFICATIONS PASSED\n'
