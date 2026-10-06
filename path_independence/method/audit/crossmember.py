from common import *
MOD = 'CanESM5-1'
rows = []
ep = gwl_ends(MOD, 'ssp370'); fp = kt.member_files(RAMIP, MOD, 'ssp370', 'tas')
WP = {r: kt.window_vals(fp[r], 'tas', ep[r])[0] for r in sorted(set(ep) & set(fp))}
_, lat, lon, _ = kt.window_vals(fp[next(iter(WP))], 'tas', 2050); lat, lon = lat.values, lon.values
land = kt.land_mask(lat, lon); landf = land.ravel()
runs = sorted(WP)
for i, ri in enumerate(runs):                      # GWL-matched same-forcing null (parent vs parent)
    for rj in runs[i + 1:]:
        rej = kt.block_perm_map(WP[ri], WP[rj]) < 0.05
        rows.append(dict(arm='null', kind='null', run_p=ri, run_a=rj, land=rej[landf].mean(),
                         **{f'inbox_{a}': rej[(box_mask(lat, lon, b) & land).ravel()].mean() for a, b in BOXES.items()}))
print('null done', flush=True)
for arm, exp in [('eas', 'ssp370-EAS126aer'), ('sas', 'ssp370-SAS126aer')]:
    ea = gwl_ends(MOD, exp); fa = kt.member_files(RAMIP, MOD, exp, 'tas')
    WA = {r: kt.window_vals(fa[r], 'tas', ea[r])[0] for r in runs if r in ea and r in fa}
    boxf = (box_mask(lat, lon, BOXES[arm]) & land).ravel()
    for ri in runs:
        for rj in WA:
            rej = kt.block_perm_map(WP[ri], WA[rj]) < 0.05
            rows.append(dict(arm=arm, kind='same' if ri == rj else 'cross', run_p=ri, run_a=rj, land=rej[landf].mean(),
                             inbox=rej[boxf].mean(), outside=rej[landf & ~boxf].mean()))
    print(arm, 'done', flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT + 'crossmember.csv', index=False)
print(df.groupby(['arm', 'kind'])[['land', 'inbox', 'outside']].agg(['mean', 'std', 'count']).round(3).to_string())
print('null in-box by arm box:'); print(df[df.kind == 'null'][[c for c in df if c.startswith('inbox_')]].mean().round(3).to_string())
