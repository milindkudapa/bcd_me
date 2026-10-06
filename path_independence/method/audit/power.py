from common import *
import time
MOD = 'CanESM5-1'
rng = np.random.default_rng(0)
for k in range(2):
    t = time.time(); kt.block_perm_map(rng.standard_normal((3000, 120)), rng.standard_normal((3000, 120)))
    print(f'block_perm_map 3000 px: {time.time()-t:.1f} s' + (' (incl. JIT)' if k == 0 else ''), flush=True)
rows = []
for var in ['tas', 'pr']:
    fns = kt.member_files(RAMIP, MOD, 'ssp370', var); runs = sorted(fns)
    W = [kt.window_vals(fns[r], var, 2050)[0] for r in runs[:7]]
    _, lat, lon, _ = kt.window_vals(fns[runs[0]], var, 2050); lat, lon = lat.values, lon.values
    land = (kt.land_mask(lat, lon) & (lat > -60)[:, None]).ravel()
    trop = np.repeat(np.abs(lat) <= 23.5, len(lon))
    clim12 = W[6].reshape(-1, 10, 12).mean(1)             # common climatology from an independent 7th member
    clim = np.tile(clim12, 10); amp = clim12.std(1)        # seasonal amplitude per pixel
    sig = (W[6] - clim).std(1)                             # monthly-anomaly sd per pixel
    if var == 'tas':
        shifts = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]; mk = lambda B, s: B + np.float32(s)
    else:
        shifts = [0, 0.05, 0.1, 0.2, 0.3, 0.5]; mk = lambda B, s: B * np.float32(1 + s)
    for mode in ['raw', 'anom']:
        for i, j in [(0, 1), (2, 3), (4, 5)]:
            A, B = W[i][land], W[j][land]
            for s in shifts:
                Bs = mk(B, s)
                A2, B2 = (A - clim[land], Bs - clim[land]) if mode == 'anom' else (A, Bs)
                rej = kt.block_perm_map(A2, B2) < 0.05
                a, g = amp[land], sig[land]
                rows.append(dict(var=var, mode=mode, pair=f'{i}{j}', shift=s, all=rej.mean(),
                                 tropics=rej[trop[land]].mean(), extratropics=rej[~trop[land]].mean(),
                                 amp_lo=rej[a < np.median(a)].mean(), amp_hi=rej[a >= np.median(a)].mean()))
            print(var, mode, i, j, 'done', flush=True)
    print(var, 'land px', land.sum(), 'median seasonal amp', np.median(amp[land]), 'median anomaly sd', np.median(sig[land]))
df = pd.DataFrame(rows); df.to_csv(OUT + 'power.csv', index=False)
print(df.groupby(['var', 'mode', 'shift'])[['all', 'tropics', 'extratropics', 'amp_lo', 'amp_hi']].mean().round(3).to_string())
