# Prediction for the horizon x eps cross product, written before the 14 cells land

Standing at 00:40 on 2026-09-08, from the 12 cells already measured.

1. **eps partly buys the horizon back on Burrau, and `t` and `eps` are not separable.**
   `near-field` runs 37.51x (t13) -> 2.89x (t23) -> 1.00x (t50) at eps=1e-2, and the two cells
   already in hand say the ladder shifts rather than falls: t23 eps=1e-1 reads **13.72x** against
   t23 eps=1e-2's 2.89x. So I expect t50 eps=1e-1 to be well above its 1.00x at eps=1e-2 -- call
   it 3-15x -- and t23/t50 at eps=1e-3 to be at or below 1x, degenerate the way t13 eps=1e-3
   already is (0.30x). Mechanism: the sea grows with `t` (things terminate and stop agreeing) and
   shrinks with `eps`; only the *ratio* of the two matters, so the grid should be closer to
   diagonal than to either axis.

2. **The sea charts stay flat near 1x in every cell.** `config_stability` and `tilt_plambda` sit
   in 0.82-1.16x across all six cells measured. I do not expect any (t, eps) to lift them above
   ~1.5x, because their `dp/u` ceiling is 1.08-1.58x -- there is nothing for a policy to find. If
   one of them exceeds 2x I have the mechanism wrong.

3. **`Undetermined` stays at exactly 0.0000 in all 16.** It is 0.0000 in every cell so far,
   including t50, where the budget has the most room to run out. The user named this as the thing
   most likely to break the saving; on this evidence it does not fire at all at these settings.

4. **The floored-child `alpha_area` p50 stays at or below 0.0000 in the degenerate cells.** The
   empty-mask defect returns exactly 0.0000, indistinguishable from a space-filling sea, so I
   expect every degenerate cell to read +0.0000 with a small `n` -- which is the *defect*
   signature and not a measurement of dimension.

## The thing most likely to prove me wrong

That the grid is diagonal at all. If `eps=1e-1` at `t=50` reads ~1x like its eps=1e-2 sibling,
then `t` dominates and loosening the tolerance does not buy the horizon back -- which would mean
the saving is a property of the horizon alone and no tolerance setting recovers it.
