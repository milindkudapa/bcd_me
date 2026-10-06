from common import *
import glob, collections
# 1. RAMIP file time spans per model/experiment
print('--- RAMIP data spans (last year per experiment)')
span = collections.defaultdict(dict)
for f in glob.glob(RAMIP + '*/tas_Amon_*.zarr'):
    mod = f.split('/')[-2]; m = re.search(r'_(ssp[^_]+)_r\d+i\d+p\d+f\d+_(.+?)\.zarr$', f)
    if m: span[mod].setdefault(m.group(1), set()).add(m.group(2))
for mod in sorted(span):
    print(f'{mod:18s}', {e: sorted(y) for e, y in sorted(span[mod].items())})
# 2. GMST mismatch between the paired windows
g = notebook_setup(); JOBS, ARMS = g['JOBS'], g['ARMS']
rows = []
for (anc, arm, mod), (aroot, amod, aexp, bmod, prs) in JOBS.items():
    for ra, ea, rb, eb in prs:
        ga = float(GMST.gmst_anom.sel(model=amod, exp=aexp, member=member_idx(amod, aexp, ra), year=slice(ea - 9, ea)).mean())
        gb = float(GMST.gmst_anom.sel(model=bmod, exp=ARMS[arm], member=member_idx(bmod, ARMS[arm], rb), year=slice(eb - 9, eb)).mean())
        rows.append(dict(anchor=anc, arm=arm, model=mod, g_anchor=ga, g_arm=gb, diff=ga - gb))
d = pd.DataFrame(rows)
print('\n--- realised 10-yr GMST anomaly in the paired windows (K): level and anchor-minus-arm mismatch')
print(d.groupby('anchor').agg(level_anchor=('g_anchor', 'mean'), level_arm=('g_arm', 'mean'), diff_mean=('diff', 'mean'), diff_absmean=('diff', lambda s: s.abs().mean()), diff_absmax=('diff', lambda s: s.abs().max()), n=('diff', 'size')).round(3).to_string())
print(d.groupby(['anchor', 'model'])['diff'].agg(['mean', lambda s: s.abs().max()]).round(3).to_string())
# 3. gdp vs tas redundancy on the cached 1deg maps
print('\n--- gdp vs tas p-maps (land): fraction of pixels with |p_gdp - p_tas| > 0.1, and correlation')
out = []
for ft in sorted(glob.glob(AUX + 'diagnostics/regional_arms/*_tas_*.nc')):
    fg = ft.replace('_tas_', '_gdp_')
    if not os.path.exists(fg): continue
    a, b = xr.open_dataset(ft).p_blockperm, xr.open_dataset(fg).p_blockperm
    land = kt.land_mask(a.lat.values, a.lon.values)
    pa, pb = a.values[:, land], b.values[:, land]
    out.append(dict(file=os.path.basename(ft), frac_diff=(np.abs(pa - pb) > 0.1).mean(), corr=np.corrcoef(pa.ravel(), pb.ravel())[0, 1],
                    agree_rej=((pa < .05) == (pb < .05)).mean()))
o = pd.DataFrame(out); print(o.describe().loc[['mean', 'min', 'max']].round(3).to_string())
# 4. FDR resolution: with p >= 1/(NDRAWS+1) the BH step-up needs >= N*pmin/FDR pixels at the floor
N = int(kt.land_mask(kt.GRID_OUT['lat'], kt.GRID_OUT['lon']).sum())
pmin = 1 / (kt.NDRAWS + 1)
print(f'\n--- FDR resolution: N land px (1deg, incl. Antarctica) = {N}; p floor = {pmin}; BH needs >= {N*pmin/0.2:.0f} px at the floor -> min nonzero FDR fraction {pmin/0.2:.4f}')
vals = []
for f in glob.glob(AUX + 'diagnostics/regional_arms/*_tas_*.nc'):
    ds = xr.open_dataset(f); land = kt.land_mask(ds.lat.values, ds.lon.values)
    p = ds.p_blockperm.where(xr.DataArray(land, coords={'lat': ds.lat, 'lon': ds.lon}))
    sig = kt.sig_fdr(p, 0.2); vals += list(((sig == 1).sum(('lat', 'lon')) / land.sum()).values)
vals = np.array(vals); nz = vals[vals > 0]
print(f'per-pair FDR fractions over {len(vals)} tas pairs: nonzero in {len(nz)} pairs; min nonzero = {nz.min():.4f}; values: {np.unique(np.round(nz, 4))[:12]}')
print(f'for the box sizes: eas ~{int((box_mask(kt.GRID_OUT["lat"], kt.GRID_OUT["lon"], BOXES["eas"]) & kt.land_mask(kt.GRID_OUT["lat"], kt.GRID_OUT["lon"])).sum())} land px -> needs {int((box_mask(kt.GRID_OUT["lat"], kt.GRID_OUT["lon"], BOXES["eas"]) & kt.land_mask(kt.GRID_OUT["lat"], kt.GRID_OUT["lon"])).sum()*pmin/0.2)+1} px at floor; native CanESM5-1 sas box 65 px -> 1 px')
