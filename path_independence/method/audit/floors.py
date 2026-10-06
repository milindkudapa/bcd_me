from common import *
PARENT = {m: 'ssp370' for m in ['CanESM5-1', 'CNRM-ESM2-1', 'EC-Earth3-AerChem', 'GISS-E2-1-G', 'MIROC6', 'MRI-ESM2-0', 'NorESM2-LM', 'UKESM1-0-LL']}
PARENT['CESM2'] = 'ssp370-LE'
PAIRS = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)]
rows = []
for mod, exp in PARENT.items():
    fns = kt.member_files(RAMIP, mod, exp, 'tas'); runs = sorted(fns); W = {}
    def win(i):
        if i not in W:
            v, lat, lon, n = kt.window_vals(fns[runs[i]], 'tas', 2050); W[i] = (v if n == 120 else None, lat.values, lon.values)
        return W[i]
    for i, j in PAIRS:
        if j >= len(runs): continue
        (A, lat, lon), (B, _, _) = win(i), win(j)
        if A is None or B is None: continue
        land = kt.land_mask(lat, lon); noant = land & (lat > -60)[:, None]; ocean = ~land & (np.abs(lat) < 60)[:, None]
        rej = kt.block_perm_map(A, B) < 0.05
        rows.append(dict(model=mod, pair=f'{runs[i]}|{runs[j]}', land_all=rej[land.ravel()].mean(), land_noant=rej[noant.ravel()].mean(),
                         ocean=rej[ocean.ravel()].mean(), npx_land=int(land.sum())))
    print(mod, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT + 'floors.csv', index=False)
print(df.groupby('model')[['land_all', 'land_noant', 'ocean']].agg(['mean', 'std']).round(3).to_string())
print(df.groupby('model').npx_land.first().to_string())
