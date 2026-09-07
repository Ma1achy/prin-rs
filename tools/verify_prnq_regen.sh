#!/bin/sh
# Verify a .prnq regeneration over the WHOLE corpus, not a sample.
#
# Usage: tools/verify_prnq_regen.sh [dir=results/charts]
#
# RESULTS §13 already reported "reproduces bitwise" from eleven dumps while nineteen had
# moved, so this walks every .prnq in results/charts and diffs it against HEAD. The one
# fields allowed to move are `wall_seconds` and the `[kernel: ...]` stamp; anything else --
# and especially the `decision` column, which is what moves when a parameter is wrong and
# the tree is not -- is a failure. A same-size binary difference is exactly what a moved
# decision column looks like, so file size is not the test.
#
# The kernel stamp landed 2026-09-06, so every dump committed before it lacks the suffix
# entirely and a first re-run differs on every file for a FORMAT reason. It is stripped for
# the comparison and reported separately, per file, on its own line: discarding it silently
# would put the one field that says which integrator ran into the same blind spot the stamp
# was added to close. A kernel that genuinely changed also moves the record block, which is
# what the comparison reads -- the stamp line is the label, not the evidence.
cd "$(git rev-parse --show-toplevel)" || exit 1
ok=0; wall=0; bad=0; stamped=0
for f in "${1:-results/charts}"/*.prnq; do
  git show "HEAD:$f" > /tmp/vg_old.prnq 2>/dev/null || { echo "  NEW (not in HEAD): $f"; continue; }
  if cmp -s /tmp/vg_old.prnq "$f"; then ok=$((ok+1)); continue; fi
  strings /tmp/vg_old.prnq | sed -e 's/wall_seconds=[0-9.]*//' -e 's/ *\[kernel:[^]]*\]//' > /tmp/vg_a.txt
  strings "$f"            | sed -e 's/wall_seconds=[0-9.]*//' -e 's/ *\[kernel:[^]]*\]//' > /tmp/vg_b.txt
  ka=$(strings /tmp/vg_old.prnq | sed -n 's/.*\[kernel: \([^]]*\)\].*/\1/p' | head -1)
  kb=$(strings "$f"             | sed -n 's/.*\[kernel: \([^]]*\)\].*/\1/p' | head -1)
  if [ "$ka" != "$kb" ]; then
    stamped=$((stamped+1))
    echo "  kernel stamp: $f"
    echo "      was: ${ka:-<absent, predates 2026-09-06>}"
    echo "      now: ${kb:-<absent>}"
  fi
  if cmp -s /tmp/vg_a.txt /tmp/vg_b.txt; then
    wall=$((wall+1))
  else
    bad=$((bad+1))
    echo "  MOVED beyond wall_seconds: $f"
    diff /tmp/vg_a.txt /tmp/vg_b.txt | head -6
  fi
done
echo
echo "identical=$ok  differ only in wall_seconds/kernel stamp=$wall  MOVED=$bad  (total $((ok+wall+bad)))"
echo "kernel stamp differs on $stamped file(s)"
[ "$bad" -eq 0 ] || echo "*** the decision column or the record block moved -- do not commit until this is understood"
