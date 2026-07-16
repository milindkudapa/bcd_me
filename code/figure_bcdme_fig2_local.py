#!/usr/bin/env python
"""Fig-2-style: change in 20-yr mean daily T, GWL0.61->GWL2, by IPCC AR6 region, split by SSP.
Mirrors BCD-ME Fig 2a (top regions = largest inter-SSP mean diff)."""
import xarray as xr, numpy as np, pandas as pd, regionmask
import matplotlib; matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import patches as mpatches

IN = '/burg/glab/users/mck2199/summer_data/aux/bcdme_gwl061to2_change_ERA5.nc'
OUT = '/burg/glab/users/mck2199/summer_data/figures/bcdme_fig2_local.png'
SSP_COLORS = {'ssp245': '#f69320', 'ssp370': '#df0000', 'ssp585': '#980002'}
SSPS = list(SSP_COLORS)

da = xr.open_dataarray(IN).load()
# normalize experiment labels (zarr stores like 'ssp245' or 'SSP245')
exp = np.array([str(e).lower().replace('-', '') for e in da.experiment.values])

# AR6 land regions -> region index per gridcell, then mean per member per region
regions = regionmask.defined_regions.ar6.land
mask = regions.mask(da.lon, da.lat)            # (lat, lon) region id
rows = []
for rid in np.unique(mask.values[~np.isnan(mask.values)]):
    rid = int(rid)
    name = regions[rid].name
    cell = da.where(mask == rid)
    w = np.cos(np.deg2rad(da.lat))
    rmean = cell.weighted(w.broadcast_like(cell.isel(idv=0)).fillna(0)).mean(('lat', 'lon'))
    for i in range(da.sizes['idv']):
        if exp[i] in SSP_COLORS:
            rows.append({'region': name, 'ssp': exp[i], 'val': float(rmean.isel(idv=i))})
df = pd.DataFrame(rows).dropna()

# rank regions by spread of per-SSP means (inter-SSP difference)
piv = df.groupby(['region', 'ssp'])['val'].mean().unstack()
spread = (piv.max(1) - piv.min(1)).sort_values(ascending=False)
# preprint Fig 2a: 5 biggest + 5 smallest SSP mean diffs
biggest = list(spread.head(5).index)
smallest = list(spread.tail(5).index)
top = biggest + smallest               # top->bottom: biggest group, then smallest

fig, ax = plt.subplots(figsize=(8, 9))
# top[::-1] -> y=0 bottom (smallest5) ... y=9 top (biggest1)
for r_idx, region in enumerate(top[::-1]):
    for s_idx, ssp in enumerate(SSPS):
        d = df.query('region==@region and ssp==@ssp')['val']
        if len(d) == 0:
            continue
        y = r_idx + (s_idx - 1) * 0.25
        bp = ax.boxplot(d, positions=[y], widths=0.22, vert=False,
                        patch_artist=True, whis=[5, 95],
                        flierprops={'marker': '.', 'markersize': 3})
        for it in ['boxes', 'whiskers', 'caps', 'fliers', 'medians']:
            plt.setp(bp[it], color=SSP_COLORS[ssp])
        plt.setp(bp['boxes'], facecolor=SSP_COLORS[ssp], alpha=0.7)
        plt.setp(bp['medians'], color='white')
        ax.plot(d.mean(), y, 'x', color='k', markersize=5, markeredgewidth=1.2)

ax.set_yticks(range(len(top)))
ax.set_yticklabels(top[::-1], fontsize=8)
ax.axvline(0, color='k', lw=0.8, ls='--', alpha=0.6)
# divider between smallest group (y 0-4) and biggest group (y 5-9)
ax.axhline(4.5, color='k', lw=1.2)
ax.text(0.98, 0.97, '5 biggest SSP mean diffs.', transform=ax.transAxes,
        ha='right', va='top', fontsize=8, bbox=dict(fc='white', ec='0.6', lw=0.5))
ax.text(0.98, 0.47, '5 smallest SSP mean diffs.', transform=ax.transAxes,
        ha='right', va='top', fontsize=8, bbox=dict(fc='white', ec='0.6', lw=0.5))
ax.set_xlim(-1, 5)
ax.set_xlabel('Change GWL2 − GWL0.61 [K]')
ax.set_title('Change in 20-year mean of daily T', fontsize=12, loc='center')
ax.legend(handles=[mpatches.Patch(facecolor=SSP_COLORS[s], label=s.upper()) for s in SSPS],
          loc='upper left', fontsize=9)
plt.tight_layout()
plt.savefig(OUT, dpi=150, bbox_inches='tight')
print('saved', OUT)
