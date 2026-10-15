"""Audit of the Burke GDP numbers (pipeline/ramip_impacts_gdp_country.ipynb and the 1deg `gdp` KS variable).

1. Burke, Hsiang & Miguel's main regression re-estimated from their replication panel (raw/burke2015/): exact
   coefficients and 1000 country-bootstrap draws, as in their GenerateBootstrapData.do.
2. The country dg recomputed four ways, cumulatively: V0 as published (rounded coefficients, each model's own absolute
   temperature); V1 exact coefficients; V2 a common temperature level per country (multi-model mean), which removes the
   models' disagreement about absolute temperature; V3 an observation-anchored level, Burke's own construction (UDel
   1980-2010 population-weighted mean + CMIP5 pattern-scaled warming to 2 K, response capped at 30 C).
3. dg tested against the same-forcing null of dg itself, per country and for the world total.
4. KS on monthly g(T) vs KS on monthly tas on the 1deg global-arm maps: how often the transform changes the test.
Outputs: out/gdp_check_*.csv, figures/audit/gdp_check.png.
"""
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import *

BUR = dir_list['raw'] + 'burke2015/replication/data/input/'
IMP = AUX + 'diagnostics/impacts/'
KS1 = AUX + 'diagnostics/ramip_gwl/ks_1deg/'
FIG = AUX + '../figures/audit/'
GWL_BASE = 0.65                          # GMST of 1980-2010 above 1850-1900, K (observational datasets: ~0.6-0.7)
FIX = {'ZAR': 'COD', 'ROM': 'ROU'}       # Burke's country codes -> ISO3
ARMS_SHOW = ['glob', 'eas', 'sas', 'afr', 'nae']
pd.set_option('display.width', 220)


def g(T, b1, b2, cap=None):
    T = np.minimum(T, cap) if cap is not None else T
    return b1 * T + b2 * T ** 2


# ------------------------------------------------------------------------------------------- 1. the regression
pan = pd.read_csv(BUR + 'GrowthClimateDataset.csv').copy()
pan['UDel_temp_popweight_2'] = pan.UDel_temp_popweight ** 2
core = ['UDel_temp_popweight', 'UDel_temp_popweight_2', 'UDel_precip_popweight', 'UDel_precip_popweight_2']
d = pan.dropna(subset=['growthWDI'] + core).reset_index(drop=True)
print(f'regression sample: {len(d)} country-years, {d.iso_id.nunique()} countries (published: 6584, 166)')

# country FE + country linear and quadratic trends partialled out country by country (Frisch-Waugh-Lovell),
# year FE kept as regressors; a bootstrap draw then only restacks the per-country blocks
yd = pd.get_dummies(d.year, drop_first=True, dtype=float).values
M = np.column_stack([d.growthWDI.values, d[core].values, yd])
t = (d.year - 1960).values.astype(float)
blocks = {}
for c, idx in d.groupby('iso_id').indices.items():
    Z = np.column_stack([np.ones(len(idx)), t[idx], t[idx] ** 2])
    blocks[c] = M[idx] - Z @ np.linalg.lstsq(Z, M[idx], rcond=None)[0]


def ols(R):
    return np.linalg.lstsq(R[:, 1:], R[:, 0], rcond=None)[0][:4]


beta = ols(np.vstack(list(blocks.values())))
# the same regression written out as Stata runs it (year FE, the panel's own _yi_/_y2_ trend columns, country FE)
TR = [c for c in pan.columns if c.startswith(('_yi_', '_y2_'))]
print('NaN in trend columns within the sample:', int(d[TR].isna().sum().sum()))
X = np.column_stack([d[core].values, yd, d[TR].fillna(0).values, pd.get_dummies(d.iso_id, dtype=float).values])
full = np.linalg.lstsq(X, d.growthWDI.values, rcond=None)[0][:4]
print('FWL vs full dummy regression:', beta, full)
assert np.allclose(beta, full, rtol=1e-6, atol=1e-12)

rng = np.random.default_rng(8675309)
keys = list(blocks)
boot = np.array([ols(np.vstack([blocks[keys[i]] for i in rng.integers(len(keys), size=len(keys))]))
                 for _ in range(1000)])
coef = pd.DataFrame(np.vstack([beta, boot]), columns=['temp', 'temp2', 'prec', 'prec2'])
coef.insert(0, 'draw', np.arange(len(coef)))            # draw 0 = point estimate
coef.to_csv(OUT + 'gdp_check_coef.csv', index=False)
B1x, B2x = beta[:2]
B1r, B2r = kt.B1, kt.B2
ok = boot[:, 1] < 0
opt = -boot[ok, 0] / (2 * boot[ok, 1])
print(f'\nexact: temp {B1x:.7f}  temp2 {B2x:.8f}  optimum {-B1x / (2 * B2x):.2f} C | rounded {B1r}/{B2r}: '
      f'optimum {-B1r / (2 * B2r):.2f} C')
print(f'bootstrap sd: temp {boot[:, 0].std():.5f}  temp2 {boot[:, 1].std():.6f}  corr {np.corrcoef(boot[:, 0], boot[:, 1])[0, 1]:.2f};'
      f' optimum 5-95%: {np.percentile(opt, 5):.1f}-{np.percentile(opt, 95):.1f} C ({(~ok).sum()} draws without a maximum)')
for T in (5, 10, 13, 15, 20, 25):
    print(f"  slope g'(T) at {T:2d} C: exact {B1x + 2 * B2x * T:+.5f}  rounded {B1r + 2 * B2r * T:+.5f} per K")

# observed baseline (Burke's: UDel 1980-2010 population-weighted mean) and their CMIP5 country/global warming ratio
b = pan[(pan.year >= 1980) & pan.UDel_temp_popweight.notna() & pan.growthWDI.notna()]
mt = b.groupby('iso').UDel_temp_popweight.agg(['mean', 'size'])
mt.index = [FIX.get(i, i) for i in mt.index]
mt = mt.sort_values('size').groupby(level=0).last()['mean']           # ZAR/COD, ROM/ROU: the longer record
chk = pd.read_csv(dir_list['raw'] + 'burke2015/GDPcap_NOClimateChange_RCP85_SSP5.csv')
chk.index = [FIX.get(i, i) for i in chk.ISO3]
print('baseline vs Burke\'s published meantemp: max |diff|', float((mt - chk.meantemp).abs().max()), 'K over',
      int((mt - chk.meantemp).notna().sum()), 'countries')
tc = pd.read_csv(BUR + 'CCprojections/CountryTempChange_RCP85.csv')
tconv = tc.assign(iso=tc.GMI_CNTRY.map(lambda c: FIX.get(c, c))).groupby('iso').Tconv.mean()


def t_obs(base=GWL_BASE):
    return (mt + (GWL - base) * tconv).dropna()


# ------------------------------------------------------------------------------------------- 2. dg four ways
pairs = pd.read_csv(IMP + 'burke_country_pairs.csv')
summ = pd.read_csv(IMP + 'burke_country_summary.csv')
econ = summ.drop_duplicates('iso3').set_index('iso3')[['name', 'gdp', 'pop_wb']]
pairs = pairs.join(econ, on='iso3')
hat = g(pairs.T_parent + pairs.dT, B1r, B2r) - g(pairs.T_parent, B1r, B2r)
print(f'\nstored dg vs g(window mean) difference: max |diff| {float((hat - pairs.dg).abs().max()):.2e} '
      f'(typical |dg| {float(pairs.dg.abs().median()):.1e}; the rest is the across-year Jensen term)')
Tcom = pairs[pairs.arm == 'glob'].groupby(['model', 'iso3']).T_parent.mean().groupby('iso3').mean()
pairs['T_common'] = pairs.iso3.map(Tcom)
pairs['T_obs'] = pairs.iso3.map(t_obs())
V = {'V0 published': pairs.dg,
     'V1 exact coef': g(pairs.T_parent + pairs.dT, B1x, B2x) - g(pairs.T_parent, B1x, B2x),
     'V2 common level': g(pairs.T_common + pairs.dT, B1x, B2x) - g(pairs.T_common, B1x, B2x),
     'V3 obs level': g(pairs.T_obs + pairs.dT, B1x, B2x, 30) - g(pairs.T_obs, B1x, B2x, 30)}
common = pairs.T_obs.notna() & pairs.gdp.notna()
gall, gcom = econ.gdp.sum(), econ.gdp[econ.index.isin(t_obs().index)].sum()
print(f'countries with an observed baseline and World Bank GDP: {pairs[common].iso3.nunique()} '
      f'({100 * gcom / gall:.1f} % of world GDP); capped at 30 C: {sorted(t_obs()[t_obs() > 30].index)}')


def totals(dg, mask):
    x = pairs.assign(v=dg)[mask].groupby(['arm', 'model', 'iso3']).agg(v=('v', 'mean'), gdp=('gdp', 'first'))
    return (x.v * x.gdp / 1e9).groupby(level=['arm', 'model']).sum()


per_model = pd.DataFrame({'V0 published, all countries': totals(pairs.dg, pairs.gdp.notna()),
                          **{k: totals(v, common) for k, v in V.items()}})
per_model.to_csv(OUT + 'gdp_check_totals.csv')
tab = per_model.groupby(level='arm').agg(['mean', 'std']).loc[ARMS_SHOW].round(1)
print('\nworld dGDP, arm minus parent at matched 2 C (bn constant-2015 US$ per year of growth; mean and sd over models):')
print(tab.to_string())
print(f"  V0 by the window-mean formula (rounded coefficients, no across-year Jensen term): "
      f"{totals(hat, common).loc['glob'].mean():.1f} bn, so V0 -> V1 is the coefficients")
for k in ('V0 published', 'V3 obs level'):
    gw = per_model.loc['glob', k].mean() * 1e9 / gcom
    print(f'  {k}: GDP-weighted dg {100 * gw:+.3f} pp per year; held for 30 years it compounds to {100 * ((1 + gw) ** 30 - 1):+.2f} % of GDP')
for base in (0.5, 0.8):
    tb = t_obs(base)
    v = g(pairs.iso3.map(tb) + pairs.dT, B1x, B2x, 30) - g(pairs.iso3.map(tb), B1x, B2x, 30)
    print(f'  V3 with the 1980-2010 GMST offset at {base} K: glob {totals(v, common).loc["glob"].mean():.1f} bn')

# coefficient uncertainty: dg is linear in (b1, b2), so each draw's world total is b1*A + b2*B
def lin_parts(T0, T1, cap=None):
    T0 = np.minimum(T0, cap) if cap else T0
    T1 = np.minimum(T1, cap) if cap else T1
    A = totals(T1 - T0, common).groupby(level='arm').mean()
    B = totals(T1 ** 2 - T0 ** 2, common).groupby(level='arm').mean()
    return A, B


unc = {}
for name, (T0, T1, cap) in {'V1 exact coef': (pairs.T_parent, pairs.T_parent + pairs.dT, None),
                            'V3 obs level': (pairs.T_obs, pairs.T_obs + pairs.dT, 30)}.items():
    A, B = lin_parts(T0, T1, cap)
    draws = np.outer(boot[:, 0], A) + np.outer(boot[:, 1], B)
    unc[name] = pd.DataFrame({'point': B1x * A + B2x * B, 'p05': np.percentile(draws, 5, axis=0),
                              'p95': np.percentile(draws, 95, axis=0), 'frac_negative': (draws < 0).mean(0)}, index=A.index)
    print(f'\n{name}: model-mean world dGDP with Burke coefficient uncertainty (1000 country-bootstrap draws), bn/yr')
    print(unc[name].loc[ARMS_SHOW].round(2).to_string())

# levels: where each model puts the big economies at 2 C, and what that does to the slope
gl = summ[(summ.arm == 'glob')].dropna(subset=['dgdp'])
top = gl.reindex(gl.dgdp.abs().sort_values(ascending=False).index).iso3.head(12).tolist()
lev = pairs[(pairs.arm == 'glob') & pairs.iso3.isin(top)].groupby(['iso3', 'model']).T_parent.mean().unstack('model')
levels = pd.DataFrame({'name': econ.name.reindex(top).values,
                       'T_min': lev.min(axis=1).reindex(top).values, 'T_mmm': lev.mean(axis=1).reindex(top).values,
                       'T_max': lev.max(axis=1).reindex(top).values, 'T_obs_1980_2010': mt.reindex(top).values,
                       'T_obs_at_2K': t_obs().reindex(top).values}, index=top)
levels['slope_min'] = B1x + 2 * B2x * levels.T_max
levels['slope_max'] = B1x + 2 * B2x * levels.T_min
levels['slope_obs'] = B1x + 2 * B2x * levels.T_obs_at_2K
for k, v in V.items():
    x = pairs.assign(v=v)[(pairs.arm == 'glob') & common].groupby(['model', 'iso3']).v.mean().groupby('iso3').mean()
    levels[f'dgdp_bn {k[:2]}'] = (x * econ.gdp).reindex(top).values / 1e9
levels.to_csv(OUT + 'gdp_check_levels.csv')
print('\nglobal cleanup, 12 largest economies by |dGDP|: absolute population-weighted temperature at matched 2 C '
      '(range over models), observed, slope of g, and dGDP per variant (bn/yr):')
print(levels.round(3).to_string())

# ------------------------------------------------------------------------------------------- 3. dg vs its own null
null = pd.read_csv(IMP + 'burke_country_null.csv').join(econ, on='iso3')
n0 = null.groupby('model').pair.nunique()
print('\nnull pairs per model:', n0.to_dict())
good = n0[n0 >= 3].index
sd0 = null.groupby(['model', 'iso3']).dg.std().rename('sd0')
x = pairs[pairs.model.isin(good)].groupby(['arm', 'model', 'iso3']).dg.agg(['mean', 'size']).join(sd0, on=['model', 'iso3'])
x['z'] = x['mean'] / (x.sd0 / np.sqrt(x['size']))
stouffer = lambda s: s.dropna().sum() / np.sqrt(max(s.notna().sum(), 1))
Zc = x.groupby(level=['arm', 'iso3']).z.apply(stouffer)
tp = pairs[pairs.model.isin(good)].assign(v=pairs.dg * pairs.gdp / 1e9).groupby(['arm', 'model', 'run']).v.sum()
tn = null[null.model.isin(good)].assign(v=null.dg * null.gdp / 1e9).groupby(['model', 'pair']).v.sum()
zt = tp.groupby(level=['arm', 'model']).agg(['mean', 'size']).join(tn.groupby(level='model').std().rename('sd0'), on='model')
zt['z'] = zt['mean'] / (zt.sd0 / np.sqrt(zt['size']))
sig = pd.DataFrame({'countries_|Z|>1.96': Zc.abs().gt(1.96).groupby(level='arm').mean(),
                    'ks_reject_mean': summ.groupby('arm').reject.mean(),
                    'world_total_Z': zt.groupby(level='arm').z.apply(stouffer),
                    'models_negative': zt.groupby(level='arm')['mean'].apply(lambda s: f'{(s < 0).sum()}/{len(s)}')})
sig.to_csv(OUT + 'gdp_check_sig.csv')
named = Zc.loc['glob'].reindex(top).rename(index=econ.name)
print('\ndg against the same-forcing null of dg (Stouffer Z over models with >=3 null pairs; 5 % expected by chance):')
print(sig.loc[ARMS_SHOW].round(3).to_string())
print('global cleanup, Z per large economy:', named.round(1).to_dict())

# ------------------------------------------------------------------------------------------- 4. KS on g(T) vs tas
rows = []
for f in sorted(glob.glob(KS1 + 'tas_*.nc')):
    mod = os.path.basename(f)[4:-3]
    if not os.path.exists(KS1 + f'gdp_{mod}.nc'):
        continue
    a, q = xr.open_dataset(f), xr.open_dataset(KS1 + f'gdp_{mod}.nc')
    runs = [r for r in a.run.values if r in set(q.run.values)]
    a, q = a.sel(run=runs), q.sel(run=runs)
    land = kt.land_mask(a.lat.values, a.lon.values)
    trop = (np.abs(a.lat.values) <= 23.5)[:, None] & land
    same = a.p_asymptotic.values == q.p_asymptotic.values
    agree = (a.p_blockperm.values < 0.05) == (q.p_blockperm.values < 0.05)
    for name, m in [('land', land), ('tropics', trop), ('extratropics', land & ~trop)]:
        rows.append(dict(model=mod, region=name, same_statistic=same[:, m].mean(), same_decision=agree[:, m].mean(),
                         reject_tas=(a.p_blockperm.values[:, m] < 0.05).mean(),
                         reject_gdp=(q.p_blockperm.values[:, m] < 0.05).mean()))
ks = pd.DataFrame(rows)
ks.to_csv(OUT + 'gdp_check_ks.csv', index=False)
print('\nKS on monthly g(T) vs KS on monthly tas, global arm, 1deg (fraction of cells; mean over models):')
print(ks.groupby('region')[['same_statistic', 'same_decision', 'reject_tas', 'reject_gdp']].mean().round(3).to_string())

# ------------------------------------------------------------------------------------------- figure
fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 6.2), gridspec_kw={'width_ratios': [1.15, 1]})
mods = sorted(lev.columns)
col = dict(zip(mods, plt.cm.tab10(np.arange(len(mods)))))
for i, iso in enumerate(top):
    for m in mods:
        if np.isfinite(lev.loc[iso, m]):
            a1.scatter(lev.loc[iso, m], i, color=col[m], s=28, zorder=3, label=m if i == 0 else None)
    a1.scatter(levels.T_obs_at_2K[iso], i, marker='D', color='k', s=40, zorder=4,
               label='observed 1980-2010 + warming to 2 °C' if i == 0 else None)
a1.axvspan(np.percentile(opt, 5), np.percentile(opt, 95), color='0.85', zorder=0, label='optimum, 5-95 % (bootstrap)')
a1.axvline(-B1x / (2 * B2x), color='k', lw=1, label=f'optimum {-B1x / (2 * B2x):.1f} °C (exact coefficients)')
a1.axvline(-B1r / (2 * B2r), color='k', lw=1, ls=':', label=f'optimum {-B1r / (2 * B2r):.1f} °C (rounded, used now)')
a1.set_yticks(range(len(top)), levels.name)
a1.invert_yaxis()
a1.set_xlabel('population-weighted annual temperature at matched 2 °C, ssp370 parent (°C)')
a1.set_title("Where each model puts each economy: the slope of Burke's curve, and so Δg, depends on it", fontsize=10)
a1.legend(fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=4, frameon=False)

var_names = list(per_model.columns)
gp = per_model.loc['glob']
for j, v in enumerate(var_names):
    for m, val in gp[v].items():
        a2.scatter(val, j, color=col.get(m, 'grey'), s=26, zorder=3)
    a2.scatter(gp[v].mean(), j, marker='|', color='k', s=500, zorder=4)
for k, v in unc.items():
    j = var_names.index(k)
    a2.plot([v.loc['glob', 'p05'], v.loc['glob', 'p95']], [j + 0.25] * 2, color='k', lw=2)
a2.axvline(0, color='k', lw=0.8)
a2.set_yticks(range(len(var_names)), ['V0 as published (all countries)', 'V0 as published (countries with\nan observed baseline)',
                                      'V1 exact coefficients', 'V2 common temperature level',
                                      'V3 observation-anchored level\n(Burke-style, 30 °C cap)'])
a2.invert_yaxis()
a2.set_xlabel('world Δ GDP, 126aer − ssp370 at matched 2 °C (bn constant-2015 USD per year of growth)')
a2.set_title('Global cleanup: dots = models, bar = model mean,\nline = 5-95 % from Burke coefficient uncertainty', fontsize=10)
fig.tight_layout()
os.makedirs(FIG, exist_ok=True)
fig.savefig(FIG + 'gdp_check.png', dpi=150, bbox_inches='tight')
print('\nsaved', FIG + 'gdp_check.png')
