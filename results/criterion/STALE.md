# What in this directory is stale, and why

**Everything here was measured under a colouring that no longer ships.** The `error(B)`
curves, the `err_sum` in every `.qcache`, and every rendered image were produced with
`bivariate::rgb`, whose lightness ramp was linear over a window an order of magnitude too
wide and whose hue map was 2-to-1 in `n0`. The standing rule is *choose a criterion under
the colouring that will ship*, so these numbers score the criteria against the wrong
target and are kept only as the record of that.

They are not deleted. They are the measurement of what a criterion looks like under an
instrument that could not see its own signal, which is itself a finding.

## Additionally inconsistent right now

`curve_*_t13.png` / `.svg` were regenerated at **`levels = 4` (128^2)** as a validation
run for the new figure writer, while `curve_*_t20.*` are still the original
**`levels = 6` (512^2)**. The two are not comparable and neither is a result. Both are
replaced by the Stage 5 re-measure.

`.qcache` files here cannot be replayed under the new colouring: `PRQC` v1 stores
per-quad reductions and a baked `err_sum`, not per-footprint state. That is what the
footprint cache exists to fix, so this is the last time a colouring change costs an
integration.

## What in this directory is NOT stale

Nothing that is a picture of a tree rather than of a colouring: the `_wire` twins still
say where the tree cut, and the tree structure did not change. But they are drawn over
the old colouring, so read the wire, not the fill.
