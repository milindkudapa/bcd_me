from common import *
g = notebook_setup(); PARENT, parent_pairs = g['PARENT'], g['parent_pairs']
EXP = 'ssp370-EAS126aer'; rows = []
def wmean(F, w, m): return float((F * w)[m].sum() / w[m].sum())
for mod in PARENT:
    for var in ['tas', 'pr']:
        src = kt.member_files(RAMIP, mod, PARENT[mod], var); afn = kt.member_files(RAMIP, mod, EXP, var)
        jobs = [('gwl', r, ep, ea) for r, ep, _, ea in parent_pairs(mod, EXP)] + [('ym', r, 2050, 2050) for r in sorted(set(src) & set(afn))]
        cache = {}
        for kind, r, ep, ea in jobs:
            if r not in src or r not in afn: continue
            if (r, ep) not in cache:
                P = kt.window_da(src[r], var, ep); cache[(r, ep)] = (P.mean('time').values, P.lat.values, P.lon.values) if P.sizes['time'] == 120 else None
            if cache[(r, ep)] is None: continue
            Pm, lat, lon = cache[(r, ep)]; A = kt.window_da(afn[r], var, ea)
            if A.sizes['time'] != 120 or A.shape[1:] != Pm.shape: continue
            Am = A.mean('time').values; land = kt.land_mask(lat, lon); w = np.cos(np.deg2rad(lat))[:, None] * np.ones(len(lon))
            box = box_mask(lat, lon, BOXES['eas']) & land; d = Am - Pm
            row = dict(kind=kind, model=mod, var=var, run=r, parent_end=ep, arm_end=ea, d_in=wmean(d, w, box), d_out=wmean(d, w, land & ~box), d_land=wmean(d, w, land))
            if var == 'pr':
                for k, m in [('d_in', box), ('d_out', land & ~box), ('d_land', land)]: row[k] = 100 * row[k] / wmean(Pm, w, m)
            rows.append(row)
        print(mod, var, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT + 'effects_eas.csv', index=False)
print(df.groupby(['kind', 'var']).agg(models=('model', 'nunique'), pairs=('run', 'size'), d_in=('d_in', 'mean'), d_in_sd=('d_in', 'std'), d_out=('d_out', 'mean')).round(3).to_string())
print(df.groupby(['kind', 'var', 'model']).d_in.mean().unstack('model').round(3).to_string())
