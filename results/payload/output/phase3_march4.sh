#!/bin/zsh
# The three marches at the shipped default (alpha_lo 0.005) with the live_nonfinite fix, so the
# live view no longer inherits the run's non-finite verdict. Same arguments as march3, which is
# the row this is diffed against.
BIN=$1; ROOT=$2; S=/private/tmp/claude-501/-Users-malachy-src-principia-rs-test/260ec0e2-d203-4f2d-85ff-f64a85497918/scratchpad
mkdir -p $ROOT/output
for tgt in preset_shape_h1 config_stability preset_shape; do
  f=$(ls $S/payload/payload/${tgt}_t13_L6.fcache | head -1)
  echo "=== march $tgt alpha_lo=0.005 livenf  $(date +%H:%M:%S)"
  $BIN march "$f" 0.01 0.25 $ROOT 0 4 0.005 1 1 1 2>&1 | grep -E "live descent:|memory:|indicator|resolvable|panicked|error:" | cut -c1-220
done
echo "phase3_march4 done $(date +%H:%M:%S)"
