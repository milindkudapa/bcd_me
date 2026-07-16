#!/usr/bin/env python
"""Faithful port of BCD-ME authors' figure_ssps_gwls.ipynb (Fig 2).
Replicates exactly: member-balancing across SSPs, all models, ERA5, xagg
aggregation to IPCC-WGI v4 regions, 5-biggest+5-smallest inter-SSP sort,
3 panels (mean / max / std of 20-yr daily T, GWL0.61->GWL2).

Deviation from authors (documented):
  - landmask: authors use Carleton impact-region shapefile (unavailable here);
    we substitute Natural Earth land (regionmask). Affects only coastal cells.

Run with .venv312 (arraylake 1.x + xagg + matplotlib).
"""
import os, re, string, warnings
import numpy as np, pandas as pd, xarray as xr
import geopandas as gpd, xagg as xa, regionmask
import arraylake, zarr
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches

AUX = '/burg-archive/home/mck2199/summer/summer_data/aux/'
FIGDIR = '/burg/glab/users/mck2199/summer_data/figures/'
POLY_FN = AUX + 'geo_data/IPCC-WGI-reference-regions-v4.shp'
CACHE = AUX + 'diagnostics/gwl_test_data_gwl2_faithful.nc'        # post-xagg (regional)
RAWCACHE = AUX + 'diagnostics/gwl_test_rawchange_gwl2_faithful.nc'  # pre-xagg (per-gridcell)
os.makedirs(AUX + 'diagnostics', exist_ok=True)

VAR, GWL = 'tas', 2
MODELS = ['CanESM5', 'MIROC6', 'MPI-ESM1-2-LR', 'IPSL-CM6A-LR']  # cap to 4 (set None for full set)
EXPS = ['ssp245', 'ssp370', 'ssp585']
COLORS = {'ssp245': '#f69320', 'ssp370': '#df0000', 'ssp585': '#980002'}

# ── WGI v4 regions (authors: drop ocean + E/W Antarctica) ────────────────────
gdf = gpd.read_file(POLY_FN)
gdf = gdf.query('Type != "Ocean" and Acronym != "EAN" and Acronym != "WAN"')


def landmask_for(ds_grid):
    """Substitute for authors' Carleton get_landmask: Natural Earth land -> 1/0."""
    land = regionmask.defined_regions.natural_earth_v5_0_0.land_110
    m = land.mask(ds_grid.lon, ds_grid.lat)        # 0 on land, NaN ocean
    return xr.where(~np.isnan(m), 1, 0)


def build():
    al = arraylake.Client()
    session = al.get_repo("ClimateUncertaintyLab/bcd_me_qdm").readonly_session(branch="main")
    modlist = list(zarr.open_group(session.store, mode='r').group_keys())
    if MODELS is not None:
        modlist = [m for m in modlist if m in MODELS]
    print("models:", modlist, flush=True)

    ds_grid = xr.open_zarr(session.store, zarr_format=3, group=modlist[0])[['lat', 'lon']]
    landmask = landmask_for(ds_grid)

    # has_data across all models (authors cell 9)
    hd = [xr.open_zarr(session.store, zarr_format=3, group=m) for m in modlist]
    hd = [h.has_data for h in hd if 'has_data' in h]
    has_data_all = xr.concat(hd, dim='idv', join='outer').load()
    has_data_all = has_data_all.where(~np.isnan(has_data_all), 0)

    # ── member selection: MATCHED runs (same model+run under all target SSPs) ─
    # Authors balance then rely on dropna(how='any') to keep cross-SSP-matched
    # members. We enforce matching up front so path-independence is like-for-like
    # and multiple models survive (their data overlaps enough that first-N works;
    # our smaller pull needs the explicit intersection).
    htmp = has_data_all.sel(gwl=GWL, variable=VAR)
    htmp = htmp.where(htmp, drop=True)
    htmp = htmp.set_index(idv=['model', 'experiment', 'run']).to_dataframe().reset_index()
    htmp = htmp[htmp.experiment.isin(EXPS)]
    # keep (model,run) present in ALL target experiments
    nexp = htmp.groupby(['model', 'run'])['experiment'].transform('nunique')
    subset = htmp[nexp == len(EXPS)][['model', 'experiment', 'run']]
    keys = subset.set_index(['model', 'experiment', 'run']).index
    print("matched runs/model (in all SSPs):\n",
          subset.groupby(['model', 'experiment']).size(), flush=True)

    # ── load balanced members, compute stats, change to GWL ──────────────────
    dss = []
    for mod in keys.levels[0]:
        k = keys[keys.get_loc(mod)]
        if len(k) == 0:
            continue
        ds = (xr.open_zarr(session.store, zarr_format=3, group=mod)
              .set_index(idv=['model', 'experiment', 'run'])
              .sel(idv=k, gwl=[0.61, GWL], proj_base='ERA5'))[[VAR]]
        ds = ds.where(landmask)
        stats = ds.mean(('dayofyear', 'year'), skipna=False).rename({VAR: VAR + 'mean'})
        dss.append(stats.load())
        print(f"  loaded {mod}: {len(k)} members", flush=True)

    dss = xr.concat(dss, dim='idv')
    ddss = dss.sel(gwl=GWL) - dss.sel(gwl=0.61)
    # persist per-gridcell change BEFORE xagg so aggregation tweaks skip the zarr reload
    ddss.reset_index('idv').to_netcdf(RAWCACHE, engine='h5netcdf')
    print("saved rawcache", RAWCACHE, flush=True)
    return aggregate(ddss)


def _add_bnds(d):
    """Precompute lat/lon cell bounds so xagg.get_bnds returns early.
    Avoids its in-place edge mutation that fails on read-only index coords."""
    d = d.copy()
    for var in ['lat', 'lon']:
        c = d[var]
        diff = c.diff(var)
        diff = xr.concat([xr.DataArray([diff[0].values], coords={var: [c[0].values]}, dims=var),
                          diff], dim=var)
        b = xr.concat([c - 0.5 * diff, c + 0.5 * diff], dim='bnds').transpose(var, 'bnds')
        b = (b.where(b <= 180, b - 360).where(b >= -180, b + 360)) if var == 'lon' \
            else (b.where(b <= 90, 90).where(b >= -90, -90))
        d[var + '_bnds'] = b
    return d


def aggregate(ddss):
    """xagg pixel-overlap aggregate to WGI regions (authors cell 10)."""
    ddss = _add_bnds(ddss)
    wm = xa.pixel_overlaps(ddss, gdf)
    agg = xa.aggregate(ddss, wm)        # default impl (== numba result, slower)
    out = agg.to_dataset()
    for v in ['model', 'experiment', 'run']:
        out = out.assign_coords({v: (('idv'), ddss[v].values)})
    out = out.set_index(idv=['model', 'experiment', 'run']).unstack()
    out = out.stack(idv=['model', 'run']).dropna('idv', how='all')
    out.attrs = {'SOURCE': 'figure_bcdme_fig2_faithful.py'}
    out.reset_index('idv').to_netcdf(CACHE, engine='h5netcdf')
    print("saved cache", CACHE, flush=True)
    return out


if os.path.exists(CACHE):
    print("using cache", CACHE)
    ddss_out = xr.open_dataset(CACHE).set_index(idv=['model', 'run'])
elif os.path.exists(RAWCACHE):
    print("using rawcache (skip zarr reload)", RAWCACHE)
    ddss_out = aggregate(xr.open_dataset(RAWCACHE))
else:
    ddss_out = build()
ddss_out = ddss_out.dropna('idv')   # authors' how='any'; no-op now (members pre-matched)

# ── figure (authors cell 13): mean-only panel, 5 biggest + 5 smallest ────────
METRICS = {'tasmean': ('20-year mean of daily T', 'K', [-1, 5])}
BW, SHOW = 0.2, 5
# count member-runs feeding the boxplots, per SSP (members differ by SSP)
valid = ddss_out['tasmean'].sel(experiment=EXPS).notnull().any('poly_idx')  # (experiment, idv)
nruns = max(int(valid.sel(experiment=e).sum()) for e in EXPS)
nmods = len(np.unique([str(m) for m in ddss_out.model.values]))
run_desc = f'{nruns} runs/SSP from\n {nmods} models'

fig = plt.figure(figsize=(8, 9))
axs = [fig.add_subplot(1, 1, 1)]

for vi, var in enumerate(METRICS):
    pdata = ddss_out[var].sel(experiment=EXPS)
    # sort regions by max pairwise |mean_i - mean_j|
    diffs, labels = [], []
    for e1 in pdata.experiment.values:
        for e2 in pdata.experiment.values:
            diffs.append(np.abs(pdata.sel(experiment=e1).mean('idv')
                                - pdata.sel(experiment=e2).mean('idv'))
                         .drop_vars('experiment', errors='ignore'))
            labels.append(f'{e1}-{e2}')
    diffs = xr.concat(diffs, dim=pd.Index(labels, name='diffexps'))
    ds_sorted = diffs.max('diffexps').sortby(diffs.max('diffexps'), ascending=False).dropna('poly_idx')
    poly_locs = np.r_[ds_sorted.poly_idx.values[:SHOW], ds_sorted.poly_idx.values[-SHOW:]][::-1]

    ax = axs[vi]
    for ri, reg in enumerate(poly_locs):
        for ei, exp in enumerate(pdata.experiment.values):
            d = pdata.sel(experiment=exp, poly_idx=reg)
            d = d[~np.isnan(d)]
            if len(d) == 0:
                continue
            pos = ri - 1.5 * BW + BW * ei
            bp = ax.boxplot(d, positions=[pos], patch_artist=True, widths=BW,
                            whis=[5, 95], flierprops={'marker': '.', 'markersize': 2}, vert=False)
            for it in ['boxes', 'whiskers', 'fliers', 'caps']:
                plt.setp(bp[it], color=COLORS[exp])
            plt.setp(bp['medians'], color='white')
            plt.setp(bp['fliers'], markeredgecolor=COLORS[exp])
            ax.plot([d.mean()], [pos], '-x', color='k')
        ax.axhline(ri + 0.5, color='grey', linewidth=0.3)
    ax.axvline(0, color='k', linestyle='--')
    ax.set_yticks(range(len(poly_locs)),
                  labels=[re.sub(r'\&', ' & ', re.sub(r'\.', '. ', re.sub(r'\-', ' ', t)))
                          for t in ddss_out.Name.sel(poly_idx=poly_locs).values])
    ax.tick_params(axis='x', labeltop=True, top=True)
    ax.set_xlim(METRICS[var][2])
    ax.axhline(SHOW - 1 + 0.5, color='k', linewidth=2)
    if vi == 0:
        xmax = ax.get_xlim()[1]
        for text, yt in zip([f'{SHOW} biggest SSP mean diffs.', f'{SHOW} smallest SSP mean diffs.'],
                            [SHOW * 2 - 1 - 0.2, SHOW - 1 - 0.2]):
            ax.text(xmax - 0.1, yt + 0.5, text, ha='right', va='top', fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.3", fc='w', ec='k', alpha=0.8))
        ax.legend(handles=[mpatches.Patch(facecolor=COLORS[e], label=e.upper())
                           for e in pdata.experiment.values][::-1],
                  loc='upper left', title=run_desc + ' of')
    ax.set_title(f'Change in {METRICS[var][0]}', fontweight='bold')
    ax.set_xlabel(f'Change GWL{GWL} - GWL0.61 [{METRICS[var][1]}]')
    ax.text(-0.1, 1.1, string.ascii_lowercase[vi] + '.', transform=ax.transAxes,
            ha='left', va='top', fontsize=18, fontweight='bold')

plt.tight_layout()
OUT = FIGDIR + 'bcdme_fig2_faithful.png'
plt.savefig(OUT, dpi=150, bbox_inches='tight')
print('saved', OUT)
