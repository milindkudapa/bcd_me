from common import *
g = notebook_setup(); JOBS, ARMS, PARENT, parent_pairs = g['JOBS'], g['ARMS'], g['PARENT'], g['parent_pairs']
rows = []
def wmean(F, w, m): return float((F * w)[m].sum() / w[m].sum())
for mod in PARENT:
    for var in ['tas', 'pr']:
        src = kt.member_files(RAMIP, mod, PARENT[mod], var); cacheP = {}
        for arm, exp in list(ARMS.items()) + [('glob', 'ssp370-126aer')]:
            prs = parent_pairs(mod, exp) if arm == 'glob' else JOBS.get(('370', arm, mod), (None,) * 5)[4]
            if not prs: continue
            afn = kt.member_files(RAMIP, mod, exp, var)
            for r, ep, _, ea in prs:
                if r not in src or r not in afn: continue
                if (r, ep) not in cacheP:
                    P = kt.window_da(src[r], var, ep)
                    cacheP[(r, ep)] = (P.mean('time').values, (P.groupby('time.month') - P.groupby('time.month').mean('time')).std('time').values, P.lat.values, P.lon.values)
                Pm, Psd, lat, lon = cacheP[(r, ep)]
                A = kt.window_da(afn[r], var, ea)
                if A.sizes['time'] != 120 or Pm.shape != A.shape[1:]: continue
                Am = A.mean('time').values
                land = kt.land_mask(lat, lon); w = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon))
                for b, bb in (BOXES.items() if arm == 'glob' else [(arm, BOXES[arm])]):
                    box = box_mask(lat, lon, bb) & land
                    d = Am - Pm
                    row = dict(model=mod, arm=arm, box=b, var=var, run=r, parent_end=ep, arm_end=ea,
                               d_in=wmean(d, w, box), d_out=wmean(d, w, land & ~box), d_land=wmean(d, w, land), sd_in=wmean(Psd, w, box))
                    if var == 'pr':   # percent of the parent's box mean
                        for k, m in [('d_in', box), ('d_out', land & ~box), ('d_land', land), ('sd_in', box)]:
                            row[k] = 100 * row[k] / wmean(Pm, w, m)
                    rows.append(row)
            print(mod, var, arm, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT + 'effects.csv', index=False)
s = df.groupby(['var', 'arm', 'box']).agg(models=('model', 'nunique'), pairs=('run', 'size'), d_in=('d_in', 'mean'), d_in_sd=('d_in', 'std'),
                                          d_out=('d_out', 'mean'), d_land=('d_land', 'mean'), sd_in=('sd_in', 'mean'))
print(s.round(3).to_string())
