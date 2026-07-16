#!/usr/bin/env python
"""Pull GWL0.61->GWL2 change in 20-yr mean daily tas per member from BCD-ME.
Saves (member, lat, lon) change + experiment label for Fig-2-style regional boxplots.
ERA5 reanalysis only. Caps members/model to keep I/O bounded."""
import arraylake, xarray as xr, numpy as np, time
import os

TOKEN = os.environ["ARRAYLAKE_TOKEN"]
MODELS = ['CanESM5', 'MIROC6', 'MPI-ESM1-2-LR', 'IPSL-CM6A-LR']  # GFDL-ESM4 not in BCD-ME
MAXMEM = 12  # per model, enough for boxplot spread
OUT = '/burg/glab/users/mck2199/summer_data/aux/bcdme_gwl061to2_change_ERA5.nc'

c = arraylake.Client()
s = c.get_repo("ClimateUncertaintyLab/bcd_me_qdm").readonly_session(branch="main")

parts = []
for model in MODELS:
    ds = xr.open_zarr(s.store, group=model)
    # has_data: (variable, idv, gwl) -> need tas present at both 0.61 and 2.0
    hd = ds['has_data'].sel(variable='tas', gwl=[0.61, 2.0]).compute()
    ok = np.where(hd.all('gwl').values)[0][:MAXMEM]
    if len(ok) == 0:
        print(f"{model}: no members with gwl2 data, skip", flush=True)
        continue
    da = ds['tas'].sel(proj_base='ERA5', gwl=[0.61, 2.0]).isel(idv=ok)
    t = time.time()
    m = da.mean(('year', 'dayofyear')).compute()        # 20-yr daily mean
    change = (m.sel(gwl=2.0) - m.sel(gwl=0.61)).drop_vars('gwl', errors='ignore')
    change = change.assign_coords(
        experiment=('idv', ds['experiment'].isel(idv=ok).compute().values),
        model=('idv', ds['model'].isel(idv=ok).compute().values),
        run=('idv', ds['run'].isel(idv=ok).compute().values))
    parts.append(change)
    print(f"{model}: {len(ok)} members in {time.time()-t:.0f}s", flush=True)

out = xr.concat(parts, dim='idv')
out.name = 'tas_change'
out.to_netcdf(OUT, engine='h5netcdf')
print(f"\nsaved {OUT}  shape {dict(out.sizes)}")
print("experiments:", np.unique(out.experiment.values, return_counts=True))
