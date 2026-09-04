#!/bin/zsh
# The live march again under the structured-weight expiry, on the two charts that matter.
BIN=$1; ROOT=$2; S=/private/tmp/claude-501/-Users-malachy-src-principia-rs-test/260ec0e2-d203-4f2d-85ff-f64a85497918/scratchpad
mkdir -p $ROOT/output
for tgt in preset_shape_h1 near-field config_stability; do
  f=$(ls $S/payload/payload/${tgt}_t13_L6.fcache 2>/dev/null | head -1)
  [ -z "$f" ] && { echo "no fcache for $tgt"; continue; }
  echo "=== march $tgt expiry-on-structure  $(date +%H:%M:%S)"
  $BIN march "$f" 0.01 0.25 $ROOT 0 4 2>&1 | grep -E "live descent:|memory:|indicator|resolvable|panicked|error:" | cut -c1-220
done
echo "phase3_march2 done $(date +%H:%M:%S)"
