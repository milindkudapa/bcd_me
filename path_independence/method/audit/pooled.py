from common import *
g = notebook_setup(); PARENT, parent_pairs = g['PARENT'], g['parent_pairs']
rows = []
def stats(p, lat, lon, land, box):
    P2 = p.reshape(len(lat), len(lon)); da = xr.DataArray(P2, coords={'lat': lat, 'lon': lon}, dims=('lat', 'lon'))
    pb = da.where(xr.DataArray(box, coords=da.coords)).expand_dims('pair')
    return dict(land=(P2 < .05)[land].mean(), inbox=(P2 < .05)[box].mean(), outside=(P2 < .05)[land & ~box].mean(),
                fdr_in=float((kt.sig_fdr(pb, 0.2) == 1).sum()) / box.sum())
for mod in ['MRI-ESM2-0', 'MIROC6']:            # matched windows >= 2048: emission pathways fully diverged
    for arm, exp in [('sas', 'ssp370-SAS126aer'), ('eas', 'ssp370-EAS126aer')]:
        fp = kt.member_files(RAMIP, mod, PARENT[mod], 'tas'); fa = kt.member_files(RAMIP, mod, exp, 'tas')
        P, A = [], []
        for r, ep, _, ea in parent_pairs(mod, exp):
            if r in fp and r in fa:
                vp, la, lo, n1 = kt.window_vals(fp[r], 'tas', ep); va, _, _, n2 = kt.window_vals(fa[r], 'tas', ea)
                if n1 == n2 == 120: P.append(vp); A.append(va)
        lat, lon = la.values, lo.values; land = kt.land_mask(lat, lon); box = box_mask(lat, lon, BOXES[arm]) & land
        h = len(P) // 2
        p_arm = kt.block_perm_map(np.concatenate(P, 1), np.concatenate(A, 1))                     # all members pooled, n = 120*k per side
        p_nul = kt.block_perm_map(np.concatenate(P[:h], 1), np.concatenate(P[h:2 * h], 1))       # parents split in two halves
        p_one = kt.block_perm_map(P[0], A[0])                                                      # single pair, for reference
        rows += [dict(model=mod, arm=arm, kind='pooled arm vs parent', members=len(P), **stats(p_arm, lat, lon, land, box)),
                 dict(model=mod, arm=arm, kind='pooled null (parents split)', members=h, **stats(p_nul, lat, lon, land, box)),
                 dict(model=mod, arm=arm, kind='single pair', members=1, **stats(p_one, lat, lon, land, box))]
        print(pd.DataFrame(rows[-3:]).round(3).to_string(index=False), flush=True)
pd.DataFrame(rows).to_csv(OUT + 'pooled.csv', index=False)
