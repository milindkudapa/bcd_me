from common import *
g = notebook_setup(); ARMS, PARENT = g['ARMS'], g['PARENT']
rows = []
def wmean(F, w, m): return float((F * w)[m].sum() / w[m].sum())
for mod in PARENT:
    for var in ['tas', 'pr']:
        src = kt.member_files(RAMIP, mod, PARENT[mod], var); cacheP = {}
        for arm, exp in list(ARMS.items()) + [('glob', 'ssp370-126aer')]:
            afn = kt.member_files(RAMIP, mod, exp, var)
            for r in sorted(set(src) & set(afn)):
                if r not in cacheP:
                    P = kt.window_da(src[r], var, 2050)
                    cacheP[r] = (P.mean('time').values, P.lat.values, P.lon.values) if P.sizes['time'] == 120 else None
                if cacheP[r] is None: continue
                Pm, lat, lon = cacheP[r]
                A = kt.window_da(afn[r], var, 2050)
                if A.sizes['time'] != 120 or Pm.shape != A.shape[1:]: continue
                Am = A.mean('time').values
                land = kt.land_mask(lat, lon); w = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon)); d = Am - Pm
                for b, bb in (BOXES.items() if arm == 'glob' else [(arm, BOXES[arm])]):
                    box = box_mask(lat, lon, bb) & land
                    row = dict(model=mod, arm=arm, box=b, var=var, run=r, d_in=wmean(d, w, box), d_out=wmean(d, w, land & ~box), d_land=wmean(d, w, land))
                    if var == 'pr':
                        for k, m in [('d_in', box), ('d_out', land & ~box), ('d_land', land)]:
                            row[k] = 100 * row[k] / wmean(Pm, w, m)
                    rows.append(row)
            print(mod, var, arm, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT + 'effects_ym.csv', index=False)
print(df.groupby(['var', 'arm', 'box']).agg(models=('model', 'nunique'), pairs=('run', 'size'), d_in=('d_in', 'mean'), d_in_sd=('d_in', 'std'), d_out=('d_out', 'mean'), d_land=('d_land', 'mean')).round(3).to_string())
