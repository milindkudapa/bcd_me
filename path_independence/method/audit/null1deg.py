from common import *
MOD = 'CanESM5-1'; KS_DIR = AUX + 'diagnostics/regional_arms/'
runs = sorted(kt.member_files(RAMIP, MOD, 'ssp370', 'tas'))
for var in ['tas', 'pr', 'gdp']:
    fns = kt.member_files(RAMIP, MOD, 'ssp370', kt.SRC_VAR.get(var, var))
    bp, ids = [], []
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            out = kt.pair_ks_1deg(fns[runs[i]], fns[runs[j]], var, 2050, 2050)
            bp.append(out[0]); ids.append(f'{runs[i]}|{runs[j]}')
    ds = xr.Dataset({'p_blockperm': (('pair', 'lat', 'lon'), np.stack(bp).astype('float32'))},
                    coords={'pair': ids, 'lat': kt.GRID_OUT['lat'], 'lon': kt.GRID_OUT['lon']})
    ds.attrs.update(variable=var, model=MOD, exp='ssp370', window='2041-2050', note='same-forcing member pairs: every rejection is a false positive',
                    regrid=f"xESMF {kt.METHOD.get(var, 'bilinear')}, 1deg, transform-before-regrid", created_by='validation audit 2026-10-05')
    ds.to_netcdf(KS_DIR + f'null_1deg_{var}.nc')
    print(var, 'written', len(ids), 'pairs; land fraction', round(kt.land_fraction(ds, 'pair'), 4), 'fdr', round(kt.fdr_land_fraction(ds, 'pair'), 4), flush=True)
