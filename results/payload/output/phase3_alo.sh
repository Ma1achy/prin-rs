#!/bin/zsh
# Phase 3, the area floor and merging: alpha_lo x whiteness on the caches, one pinned binary.
# usage: phase3_alo.sh <pinned payload_metric binary> <root>
BIN=$1; ROOT=$2; S=/private/tmp/claude-501/-Users-malachy-src-principia-rs-test/260ec0e2-d203-4f2d-85ff-f64a85497918/scratchpad
mkdir -p $ROOT/output
for tgt in near-field preset_shape_h1 config_stability preset_shape deep_interior preset_prho; do
  f=$(ls $S/payload/payload/${tgt}_t13_L6.fcache 2>/dev/null | head -1)
  [ -z "$f" ] && f=$(ls results/payload/${tgt}_t13_L6.fcache 2>/dev/null | head -1)
  [ -z "$f" ] && { echo "no fcache for $tgt"; continue; }
  for stat in 0 1; do
    for alo in 0.2 0; do
      echo "=== live $tgt stationary=$stat alpha_lo=$alo  $(date +%H:%M:%S)"
      $BIN live "$f" tolerance 0.01 0.25 $ROOT $stat 0.01 $alo 2>&1 | grep -E "descent:|indicator|resolvable|panicked|error:" | cut -c1-220
    done
  done
  echo "=== march $tgt default (alpha_lo 0.2, merge on)  $(date +%H:%M:%S)"
  $BIN march "$f" 0.01 0.25 $ROOT 0 4 2>&1 | grep -E "live descent:|memory:|indicator|resolvable|panicked|error:" | cut -c1-220
done
echo "phase3_alo done $(date +%H:%M:%S)"
