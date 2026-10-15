from common import *
MOD = 'CanESM5-1'; KS_DIR = AUX + 'diagnostics/regional_arms/'
runs = sorted(kt.member_files(RAMIP, MOD, 'ssp370', 'tas'))
VARS = re.split('[:,]', os.environ.get('NULL_VARS', 'tas:pr:gdp'))   # ':' because sbatch --export splits on ','
SEASONAL = bool(os.environ.get('NULL_SEASONAL'))   # per-season nulls, 3-month blocks, n = 30 per side
for var in VARS:
    fns = kt.member_files(RAMIP, MOD, 'ssp370', kt.SRC_VAR.get(var, var))
    seasons = list(kt.SEASONS) if SEASONAL else [None]
    bp, ids = {s: [] for s in seasons}, []
    for i in range(len(runs)):
        for j in range(i + 1, len(runs)):
            ids.append(f'{runs[i]}|{runs[j]}')
            for s in seasons:
                out = kt.pair_ks_1deg(fns[runs[i]], fns[runs[j]], var, 2050, 2050, season=s)
                bp[s].append(out[0] if out is not None else np.full((kt.NLAT, kt.NLON), np.nan))
    if SEASONAL:
        ds = xr.Dataset({'p_blockperm': (('season', 'pair', 'lat', 'lon'),
                                         np.stack([np.stack(bp[s]) for s in seasons]).astype('float32'))},
                        coords={'season': seasons, 'pair': ids, 'lat': kt.GRID_OUT['lat'], 'lon': kt.GRID_OUT['lon']})
    else:
        ds = xr.Dataset({'p_blockperm': (('pair', 'lat', 'lon'), np.stack(bp[None]).astype('float32'))},
                        coords={'pair': ids, 'lat': kt.GRID_OUT['lat'], 'lon': kt.GRID_OUT['lon']})
    ds.attrs.update(variable=var, model=MOD, exp='ssp370', window='2041-2050', blocks='3-month, per season' if SEASONAL else '12-month', note='same-forcing member pairs: every rejection is a false positive',
                    regrid=f"xESMF {kt.METHOD.get(var, 'bilinear')}, 1deg, transform-before-regrid", created_by='validation audit 2026-10-05')
    out_fn = (AUX + 'diagnostics/regional_arms_seasonal/' + f'null_1deg_seasonal_{var}.nc') if SEASONAL else KS_DIR + f'null_1deg_{var}.nc'
    ds.to_netcdf(out_fn)
    lf = {str(s): round(kt.land_fraction(ds.sel(season=s) if SEASONAL else ds, 'pair'), 4) for s in seasons}
    print(var, 'written', len(ids), 'pairs ->', os.path.basename(out_fn), 'land fraction', lf, flush=True)
