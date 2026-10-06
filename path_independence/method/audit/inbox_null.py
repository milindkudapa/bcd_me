from common import *
KS_DIR = AUX + 'diagnostics/regional_arms/'
def box_stats(ds, boxes):
    p = ds['p_blockperm']; land = kt.land_mask(ds.lat.values, ds.lon.values)
    box = box_mask(ds.lat.values, ds.lon.values, boxes) & land; rej = (p < 0.05).values
    pb = p.where(xr.DataArray(box, coords={'lat': ds.lat, 'lon': ds.lon}))
    fdr = (kt.sig_fdr(pb, 0.2) == 1).sum(('lat', 'lon')) / box.sum()
    return dict(n_box_px=int(box.sum()), in_box=rej[:, box].mean(), fdr_in_box=float(fdr.mean()), fdr_pairs_nonzero=int((fdr > 0).sum()))
rows = []
nat = xr.open_dataset(AUX + 'diagnostics/ks_blockperm_enspairs/CanESM5-1.nc').load()
for arm, b in BOXES.items():
    rows.append(dict(box=arm, var='tas', grid='native 2.8°', **box_stats(nat, b)))
    for var in ['tas', 'pr', 'gdp']:
        rows.append(dict(box=arm, var=var, grid='1°', **box_stats(xr.open_dataset(f'{KS_DIR}null_1deg_{var}.nc').load(), b)))
df = pd.DataFrame(rows); df.to_csv(OUT + 'inbox_null.csv', index=False)
print(df.round(3).to_string(index=False))
bx = pd.read_csv(KS_DIR + 'box_summary.csv')
print('\narms (370 anchor) in-box FDR for reference:'); print(bx[bx.anchor == 370].pivot_table(index='arm', columns='var', values='fdr_in_box').round(3).to_string())
