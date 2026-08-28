"""
All candidate fields, integrated with Aarseth-Zare global regularisation.

The Benettin shadow is carried as extra rows of the SAME batch (main trajectories first, then
shadows), so each gets its own reference-body assignment from `choose_reference` -- correct,
since AZ is per-trajectory. Renormalisation, diffusion sampling and event detection all happen
at the sync boundaries, where the state is Cartesian and every trajectory shares a playhead.
"""
import numpy as np
import tb, tb_az, tb_ftle


def integrate_all_az(r0, v0, t_max, n_sync=32, eta=0.01, max_steps=30000,
                     d0=1e-8, seed=0, M=None, dtau_mode=tb_az.DTAU_MODE_DEFAULT,
                     clamp_final=tb_az.CLAMP_FINAL_DEFAULT):
    Mv = tb.M if M is None else M
    n = r0.shape[0]

    rng = np.random.default_rng(seed)
    d1 = rng.normal(size=(1,) + r0.shape[1:]); d1 /= np.linalg.norm(d1)
    Rg0 = np.sqrt(np.maximum(tb_ftle.inertia(r0), 1e-300) / Mv.sum())
    rs0 = r0 + (d0 * Rg0)[:, None, None] * np.repeat(d1, n, axis=0)

    R = np.concatenate([r0, rs0], axis=0)          # [main | shadow]
    V = np.concatenate([v0, v0.copy()], axis=0)
    E0 = tb.energy(r0, v0, 0.0)

    S = np.zeros(n); nre = np.zeros(n, dtype=np.int64)
    cnt = 0.0; St = 0.0; Stt = 0.0
    Sy = np.zeros(n); Sty = np.zeros(n)
    dmin = np.full(n, np.inf)
    t_end = np.full(n, np.nan)

    for kk in range(n_sync):
        step = t_max / n_sync
        out = tb_az.integrate_az(R, V, t_max=step, n_sync=1, eta=eta, max_steps=max_steps, M=Mv,
                                dtau_mode=dtau_mode, clamp_final=clamp_final)
        R, V = out['r'], out['v']
        dmin = np.minimum(dmin, out['dmin'][:n])

        r = R[:n]; v = V[:n]
        rs = R[n:]; vs = V[n:]

        # ---- Benettin: dimensionless separation, renormalise ----
        dr = rs - r; dv = vs - v
        Rg = np.sqrt(np.maximum(tb_ftle.inertia(r), 1e-300) / Mv.sum())
        d = np.sqrt(np.einsum('nki,nki->n', dr, dr) / Rg**2
                    + np.einsum('nki,nki->n', dv, dv) * Rg / Mv.sum())
        ok = np.isfinite(d) & (d > 1e-300)
        S[ok] += np.log(d[ok] / d0); nre[ok] += 1
        sc = np.where(ok, d0 / np.maximum(d, 1e-300), 1.0)[:, None, None]
        R[n:] = r + dr * sc
        V[n:] = v + dv * sc

        # ---- diffusion regression ----
        t_now = (kk + 1) * step
        y = np.log(np.maximum(tb_ftle.inertia(r), 1e-300))
        cnt += 1.0; St += t_now; Stt += t_now**2; Sy += y; Sty += t_now * y

        # ---- events ----
        pd = tb.pair_dists(r)
        third = np.array([2, 1, 0])[np.argmin(pd, axis=1)]
        live = np.isnan(t_end)
        if live.any():
            for b in range(3):
                sel = live & (third == b)
                if not sel.any():
                    continue
                o = [q for q in range(3) if q != b]
                mb = Mv[o].sum()
                rc = (Mv[o[0]]*r[sel][:, o[0], :] + Mv[o[1]]*r[sel][:, o[1], :]) / mb
                vc = (Mv[o[0]]*v[sel][:, o[0], :] + Mv[o[1]]*v[sel][:, o[1], :]) / mb
                dr2 = r[sel][:, b, :] - rc; dv2 = v[sel][:, b, :] - vc
                dist = np.sqrt(np.einsum('ij,ij->i', dr2, dr2))
                spec = 0.5*np.einsum('ij,ij->i', dv2, dv2) - mb/np.maximum(dist, 1e-12)
                esc = (spec > 0) & (np.einsum('ij,ij->i', dr2, dv2) > 0)
                t_end[np.nonzero(sel)[0][esc]] = t_now

    r = R[:n]; v = V[:n]
    ftle = np.where(nre > 0, S / max(t_max, 1e-12), np.nan)
    den = cnt*Stt - St*St
    diffusion = np.where(abs(den) > 1e-12, (cnt*Sty - St*Sy) / den, np.nan)
    censored = np.isnan(t_end)
    drift = np.abs((tb.energy(r, v, 0.0) - E0) / np.maximum(np.abs(E0), 1e-30))
    return dict(r=r, v=v, ftle=ftle, diffusion=diffusion, dmin=dmin,
                t_end=np.where(censored, t_max, t_end), censored=censored,
                binary_id=np.argmin(tb.pair_dists(r), axis=1), T=t_max, drift=drift)
