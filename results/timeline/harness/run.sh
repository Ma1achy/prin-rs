#!/bin/bash
# Bisect runner. One worktree, checked out at each commit in turn; the harness is REGENERATED
# from one template per commit so the parameters come from the template and never from the
# code's own defaults -- which are exactly what changes between these commits.
set -u
export PATH="$HOME/.cargo/bin:$PATH"
S="$(cd "$(dirname "$0")" && pwd)"
REPO=/Users/malachy/src/principia-rs-test
WT="$S/wt"
OUT="${OUT:-$S/out}"
RES="${RES:-1024}"
mkdir -p "$OUT"

COMMITS="${COMMITS:-2596830 483b630 077b092 e53223d 71de13f 5cc8dec 220d928 f7d2a31}"

if [ ! -d "$WT" ]; then
  git -C "$REPO" worktree add --detach "$WT" f7d2a31 >/dev/null 2>&1 || exit 1
fi

# name -> literal. Emitted only where the field exists at that commit; an absent field is
# logged as absent, which is the bisect signal.
emit_cfg() {
  local px="$1"
  local -a NAMES=(n_extra jitter_frac jitter_scheme seed t_max n_sync escape_every escape_confirm \
                  r_esc_frac escape_all_bodies escape_rule closure_k stop_on_escape dtau_mode \
                  clamp_final_step eta max_steps ref_policy lc_stable r_coll_frac stop_on_event \
                  refine_flagged refine_threshold refine_eta_factor refine_max_passes \
                  keep_copy_outcomes keep_copy_shapes keep_boundary_shapes keep_drift_hist ftle ftle_dt decode_path)
  local -a VALS=("7" "0.5" "prin_rs::ensemble::jitter::Scheme::Halton" "0" "50.0" "32" "0" "true" \
                 "12.0" "true" "prin_rs::outcome::EscapeRule::Closure(prin_rs::outcome::CLOSURE_TAU)" "1" "false" \
                 "prin_rs::integrate::az::DtauMode::PerStepInterval" "true" "0.01" "30_000" \
                 "prin_rs::integrate::az::RefPolicy::PerCopy" "true" "0.005" "true" \
                 "false" "10.0" "0.25" "3" "false" "false" "false" "false" "None" "1e-4" \
                 "prin_rs::decode::Path::DirectF64")
  local i present=() absent=()
  : > "$S/cfg.rs"
  for i in "${!NAMES[@]}"; do
    if grep -qE "^    pub ${NAMES[$i]}: " "$px"; then
      printf '        %s: %s,\n' "${NAMES[$i]}" "${VALS[$i]}" >> "$S/cfg.rs"
      present+=("${NAMES[$i]}")
    else
      absent+=("${NAMES[$i]}")
    fi
  done
  echo "  fields set: ${#present[@]}   ABSENT AT THIS COMMIT: ${absent[*]:-none}"
}

for c in $COMMITS; do
  echo "=================================================================="
  echo "== $c  $(git -C "$REPO" log -1 --format=%s "$c")"
  git -C "$WT" checkout --detach "$c" >/dev/null 2>&1 || { echo "  CHECKOUT FAILED"; continue; }
  emit_cfg "$WT/src/ensemble/pixel.rs"
  python3 - "$S/bisect_template.rs" "$S/cfg.rs" "$WT/examples/bisect_slice.rs" <<'PY'
import sys
tpl, cfg, out = sys.argv[1], sys.argv[2], sys.argv[3]
t = open(tpl).read().replace("// @@CFG@@", open(cfg).read().rstrip("\n"))
open(out, "w").write(t)
PY
  cp "$WT/examples/bisect_slice.rs" "$OUT/harness_$c.rs"
  ( cd "$WT" && cargo build --release --example bisect_slice 2>&1 | tail -3 )
  if [ ! -x "$WT/target/release/examples/bisect_slice" ]; then
    echo "  BUILD FAILED -- skipping $c"; continue
  fi
  ( cd "$WT" && ./target/release/examples/bisect_slice "$c" "$RES" "$OUT" )
done
echo "ALL DONE"
