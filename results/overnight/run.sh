#!/bin/zsh
# Overnight staged run. Cheapest-and-most-decisive first; each stage writes on completion so a
# later crash loses nothing earlier. Nothing here writes into a committed render directory.
cd /Users/malachy/src/principia-rs-test
export PATH="$HOME/.cargo/bin:$PATH"
OUT=results/overnight
mkdir -p $OUT

log() { echo "[$(date -u +%H:%M:%S)] $*" >> $OUT/RUNLOG.txt; }

log "START. Predictions were written to PREDICTIONS.md before any stage ran."

# --- STAGE 0: far at f32 vs f64, AZ vs Heggie. The conditioning question. -----------------
log "STAGE 0 begin: far 128^2, az/heggie x f64/f32"
if ./target/release/examples/stage0_far_precision 128 $OUT 400000 > $OUT/stage0_far_precision.txt 2>&1; then
  log "STAGE 0 done -> stage0_far_precision.txt"
else
  log "STAGE 0 FAILED (exit $?) -- see stage0_far_precision.txt; later stages continue"
fi

# --- STAGE 1: reversibility, per case. ------------------------------------------------------
for c in near-field far deep_interior config_stability; do
  log "STAGE 1 begin: $c at 128^2"
  if ./target/release/examples/stage1_reversibility 128 "$c" 400000 > $OUT/stage1_reversibility_$c.txt 2>&1; then
    log "STAGE 1 $c done -> stage1_reversibility_$c.txt"
  else
    log "STAGE 1 $c FAILED (exit $?) -- continuing to next case"
  fi
done

# --- STAGES 2-4: NOT BUILT. Recorded rather than silently skipped. --------------------------
cat > $OUT/STAGES_2_4_NOT_BUILT.md <<'MD'
# Stages 2, 3 and 4 are NOT BUILT, and this file says so rather than leaving a gap

Stage 0 and Stage 1 ran. Chain coordinates, TTL and IAS15 did not, because each is a new
integrator or stepper and none could be written and *validated* in the time available. An
unvalidated integrator run overnight produces numbers that have to be withdrawn, which is worse
than no numbers -- this project has withdrawn three results that way already.

**Stage 2 (chain coordinates) is gated on Stage 0 by the brief's own logic**, so it was never
going to run in the same session regardless: read `stage0_far_precision.txt` first, and in
particular its SATURATION GUARD. If the f32 row reads DEAD or SATURATED, Stage 0 does not
adjudicate the conditioning question and Stage 2's target is unproven.

What each stage needs, so the next session starts from a spec and not from scratch:

- **Chain coordinates** (Mikkola & Aarseth 1993). Two inter-particle vectors for three bodies.
  A round-off fix, so it must be run at BOTH precisions or the result is uninterpretable, and
  `far` is the predicted target.
- **TTL** (Mikkola & Aarseth 2002). Needs a deliberate mass-ratio sweep through `z8`/`z9`; the
  default charts are equal-mass and would exercise nothing. The equal-mass tie is the control.
- **IAS15** (Rein & Spiegel 2015). A reference arm, not a production candidate -- its per-lane
  variable work is the property already measured as fatal on GPU (`warps hit 1.0000`).
MD
log "STAGES 2-4 not built -- see STAGES_2_4_NOT_BUILT.md"
log "END"
echo OVERNIGHT_COMPLETE >> $OUT/RUNLOG.txt
