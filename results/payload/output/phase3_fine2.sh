#!/bin/zsh
# Completing the case for alpha_lo = 0.005 as the default: preset_shape's fine rungs, the three
# no-sea charts at 0.005, and the three marches at 0.005 (the merge shares the no-gain test).
BIN=$1; ROOT=$2; S=/private/tmp/claude-501/-Users-malachy-src-principia-rs-test/260ec0e2-d203-4f2d-85ff-f64a85497918/scratchpad
mkdir -p $ROOT/output
f=$(ls $S/payload/payload/preset_shape_t13_L6.fcache | head -1)
for alo in 0.02 0.005 0.001; do
  echo "=== live preset_shape stationary=0 alpha_lo=$alo  $(date +%H:%M:%S)"
  $BIN live "$f" tolerance 0.01 0.25 $ROOT 0 0.01 $alo 2>&1 | grep -E "descent:|indicator|resolvable|panicked|error:" | cut -c1-220
done
for tgt in near-field deep_interior preset_prho; do
  f=$(ls $S/payload/payload/${tgt}_t13_L6.fcache | head -1)
  echo "=== live $tgt stationary=0 alpha_lo=0.005  $(date +%H:%M:%S)"
  $BIN live "$f" tolerance 0.01 0.25 $ROOT 0 0.01 0.005 2>&1 | grep -E "descent:|indicator|resolvable|panicked|error:" | cut -c1-220
done
for tgt in preset_shape_h1 config_stability preset_shape; do
  f=$(ls $S/payload/payload/${tgt}_t13_L6.fcache | head -1)
  echo "=== march $tgt alpha_lo=0.005  $(date +%H:%M:%S)"
  $BIN march "$f" 0.01 0.25 $ROOT 0 4 0.005 1 1 1 2>&1 | grep -E "live descent:|memory:|indicator|resolvable|panicked|error:" | cut -c1-220
done
echo "phase3_fine2 done $(date +%H:%M:%S)"
