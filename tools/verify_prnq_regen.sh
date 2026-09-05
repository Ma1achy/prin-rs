#!/bin/sh
# Verify a .prnq regeneration over the WHOLE corpus, not a sample.
#
# Usage: tools/verify_prnq_regen.sh [dir=results/charts]
#
# RESULTS §13 already reported "reproduces bitwise" from eleven dumps while nineteen had
# moved, so this walks every .prnq in results/charts and diffs it against HEAD. The one
# field allowed to move is `wall_seconds`; anything else -- and especially the `decision`
# column, which is what moves when a parameter is wrong and the tree is not -- is a
# failure. A same-size binary difference is exactly what a moved decision column looks
# like, so file size is not the test.
cd "$(git rev-parse --show-toplevel)" || exit 1
ok=0; wall=0; bad=0
for f in "${1:-results/charts}"/*.prnq; do
  git show "HEAD:$f" > /tmp/vg_old.prnq 2>/dev/null || { echo "  NEW (not in HEAD): $f"; continue; }
  if cmp -s /tmp/vg_old.prnq "$f"; then ok=$((ok+1)); continue; fi
  strings /tmp/vg_old.prnq | sed 's/wall_seconds=[0-9.]*//' > /tmp/vg_a.txt
  strings "$f"            | sed 's/wall_seconds=[0-9.]*//' > /tmp/vg_b.txt
  if cmp -s /tmp/vg_a.txt /tmp/vg_b.txt; then
    wall=$((wall+1))
  else
    bad=$((bad+1))
    echo "  MOVED beyond wall_seconds: $f"
    diff /tmp/vg_a.txt /tmp/vg_b.txt | head -6
  fi
done
echo
echo "identical=$ok  differ only in wall_seconds=$wall  MOVED=$bad  (total $((ok+wall+bad)))"
[ "$bad" -eq 0 ] || echo "*** the decision column or the record block moved -- do not commit until this is understood"
