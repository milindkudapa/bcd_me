from common import *
g = notebook_setup(); ARMS, PARENT, win_df = g['ARMS'], g['PARENT'], g['win_df']
def box_series(da):
    lat, lon = da.lat.values, da.lon.values; land = kt.land_mask(lat, lon); w = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon))
    v = da.values
    return {b: (v * (w * (box_mask(lat, lon, bb) & land))[None]).sum((1, 2)) / (w * (box_mask(lat, lon, bb) & land)).sum() for b, bb in BOXES.items()}
def clim(series):   # tiled 12-month climatology from a list of (120,) series
    return np.tile(np.mean([s.reshape(10, 12) for s in series], axis=(0, 1)), 10)
rows, nulls = [], []
for mod in PARENT:
    for var in ['tas', 'pr']:
        src = kt.member_files(RAMIP, mod, PARENT[mod], var)
        wd = win_df[(win_df.model == mod) & win_df.usable]
        PS = {}
        for _, rw in wd.iterrows():
            if rw.run not in src: continue
            P = kt.window_da(src[rw.run], var, int(rw.anchor_end))
            if P.sizes['time'] == 120: PS[rw.run] = box_series(P)
        runs = sorted(PS)
        if len(runs) < 3: print(mod, var, 'too few parents', len(runs)); continue
        for i, ri in enumerate(runs):
            for rj in runs[i + 1:]:
                others = [r for r in runs if r not in (ri, rj)]
                for b in BOXES:
                    c = clim([PS[r][b] for r in others])
                    nulls.append(dict(model=mod, var=var, box=b, run_p=ri, run_a=rj, p=kt.ks_block_permutation_test(PS[ri][b] - c, PS[rj][b] - c)))
        for arm, exp in list(ARMS.items()) + [('glob', 'ssp370-126aer')]:
            ae = gwl_ends(mod, exp); afn = kt.member_files(RAMIP, mod, exp, var)
            for r in runs:
                if r not in ae or r not in afn: continue
                A = kt.window_da(afn[r], var, ae[r])
                if A.sizes['time'] != 120: continue
                AS = box_series(A); others = [x for x in runs if x != r]
                for b in (BOXES if arm == 'glob' else [arm]):
                    c = clim([PS[x][b] for x in others])
                    d = float((AS[b] - PS[r][b]).mean()); d = 100 * d / float(PS[r][b].mean()) if var == 'pr' else d
                    rows.append(dict(model=mod, var=var, arm=arm, box=b, run=r, p=kt.ks_block_permutation_test(PS[r][b] - c, AS[b] - c), d=d))
            print(mod, var, arm, 'done', flush=True)
df, nu = pd.DataFrame(rows), pd.DataFrame(nulls)
df.to_csv(OUT + 'boxseries.csv', index=False); nu.to_csv(OUT + 'boxseries_null.csv', index=False)
print('--- box-mean series test, GWL-matched arm vs parent: fraction of pairs with p<0.05')
print(df.groupby(['var', 'arm', 'box']).agg(models=('model', 'nunique'), pairs=('run', 'size'), rej=('p', lambda s: (s < .05).mean()), d=('d', 'mean')).round(3).to_string())
print('--- same-forcing parent pairs (null), by box'); print(nu.groupby(['var', 'box']).agg(pairs=('p', 'size'), rej=('p', lambda s: (s < .05).mean())).round(3).to_string())
print('--- arm rejection by model (tas)'); print(df[df['var'] == 'tas'].groupby(['arm', 'model']).p.agg(lambda s: (s < .05).mean()).unstack().round(2).to_string())
