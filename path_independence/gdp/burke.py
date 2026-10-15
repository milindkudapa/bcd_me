"""Burke, Hsiang & Miguel (2015) re-estimated from their replication package (raw/burke2015/replication/).

Three of the paper's growth-response regressions (pooled no-lag, pooled 5-lag, rich/poor no-lag), each with the
paper's 1000-draw country bootstrap, and the observational country baselines its projections start from.
`python burke.py` is the self-check: the pooled point estimate must reproduce the paper's coefficients exactly.
"""
import os
import sys
import numpy as np
import pandas as pd

PKG = '/burg-archive/home/mck2199/summer/bcd_me/path_independence'
sys.path.insert(0, '/burg-archive/home/mck2199/summer/bcd_me/code')
from funcs_support import get_params      # noqa: E402

_cwd = os.getcwd()
os.chdir(PKG)                             # get_params reads ./dir_list.csv
dir_list = get_params()
os.chdir(_cwd)

BUR = dir_list['raw'] + 'burke2015/replication/data/input/'
FIX = {'ZAR': 'COD', 'ROM': 'ROU'}        # Burke's country codes -> ISO3
CORE = ['UDel_temp_popweight', 'UDel_temp_popweight_2', 'UDel_precip_popweight', 'UDel_precip_popweight_2']
NAMES = ['temp', 'temp2', 'prec', 'prec2']
SPECS = ('pooled', 'lag5', 'richpoor')
PUBLISHED = (0.0127, -0.0005)             # the rounded pair quoted in the paper (ks_tools.B1, B2)
SEED = 8675309                            # Burke's seed; a different generator, so the draws differ, not their spread
CAP = 30.0                                # Burke holds the response at its 30 degC value above that
GMST_1980_2010 = 0.65                     # K above 1850-1900, Burke's baseline period (HadCRUT5 0.63, Berkeley Earth 0.69)
GMST_2015_2024 = 1.20                     # K above 1850-1900, the first decade of the RAMIP runs


def panel():
    d = pd.read_csv(BUR + 'GrowthClimateDataset.csv', low_memory=False)
    d = d[['iso', 'year', 'growthWDI', 'GDPpctile_WDIppp', 'UDel_temp_popweight', 'UDel_precip_popweight',
           'UDel_precip_popweight_2']].copy()               # precipitation squared as stored, as the do-file uses it
    d['UDel_temp_popweight_2'] = d.UDel_temp_popweight ** 2  # the do-file generates this one (Stata stores it as float)
    return d.sort_values(['iso', 'year']).reset_index(drop=True)


def design(spec, nlag=5):
    """(growth, climate regressors, country, year, column names), rows as Stata's `reg` keeps them"""
    d = panel()
    if spec == 'pooled':
        X = d[CORE].set_axis(NAMES, axis=1)
    elif spec == 'lag5':
        cols = {}
        for k in range(nlag + 1):
            lag = d.groupby('iso')[CORE + ['year']].shift(k)
            ok = (d.year - lag.year) == k                       # a lag exists only where the panel is contiguous
            for c, n in zip(CORE, NAMES):
                cols[f'{n}_L{k}'] = lag[c].where(ok)
        X = pd.DataFrame(cols)
    elif spec == 'richpoor':
        poor = (d.GDPpctile_WDIppp < 50).astype(float).where(d.GDPpctile_WDIppp.notna())
        X = pd.concat([d[CORE].mul(1 - poor, axis=0).set_axis(NAMES, axis=1),
                       d[CORE].mul(poor, axis=0).set_axis([n + 'poor' for n in NAMES], axis=1)], axis=1)
    else:
        raise ValueError(spec)
    keep = d.growthWDI.notna() & X.notna().all(axis=1)
    d, X = d[keep], X[keep]
    return d.growthWDI.values, X.values, d.iso.values, d.year.values, list(X.columns)


def _blocks(y, X, iso, year):
    """country FE and country-specific quadratic trends partialled out country by country (Frisch-Waugh-Lovell);
    year FE stay as regressors. A bootstrap draw then only restacks these blocks."""
    yd = pd.get_dummies(year, drop_first=True, dtype=float).values
    M = np.column_stack([y, X, yd])
    t = (year - 1960).astype(float)
    out = {}
    for c in np.unique(iso):
        idx = np.flatnonzero(iso == c)
        Z = np.column_stack([np.ones(len(idx)), t[idx], t[idx] ** 2])
        out[c] = M[idx] - Z @ np.linalg.lstsq(Z, M[idx], rcond=None)[0]
    return out


def fit(spec, ndraws=1000, seed=SEED, nlag=5):
    """draw 0 is the point estimate; draws 1..ndraws resample countries with replacement (Burke's bootstrap)"""
    y, X, iso, year, names = design(spec, nlag)
    k = X.shape[1]
    blocks = _blocks(y, X, iso, year)
    keys = list(blocks)

    def ols(R):
        return np.linalg.lstsq(R[:, 1:], R[:, 0], rcond=None)[0][:k]

    rng = np.random.default_rng(seed)
    rows = [ols(np.vstack([blocks[c] for c in keys]))]
    for _ in range(ndraws):
        rows.append(ols(np.vstack([blocks[keys[i]] for i in rng.integers(len(keys), size=len(keys))])))
    coef = pd.DataFrame(rows, columns=names)
    if spec == 'lag5':                                           # the cumulative (long-run) response
        for n in NAMES:
            coef[n + '_cum'] = coef[[f'{n}_L{j}' for j in range(nlag + 1)]].sum(axis=1)
    coef.insert(0, 'draw', np.arange(len(coef)))
    coef.attrs = {'n_obs': int(len(y)), 'n_countries': len(keys)}
    return coef


def coefs(coef, spec):
    """(b1, b2) over draws for the pooled or rich response, and (p1, p2) for the poor one where the spec has it"""
    if spec == 'pooled':
        return coef.temp.values, coef.temp2.values, None, None
    if spec == 'lag5':
        return coef.temp_cum.values, coef.temp2_cum.values, None, None
    return coef.temp.values, coef.temp2.values, coef.temppoor.values, coef.temp2poor.values


def moments(T0, T1, axis=-1, cap=CAP):
    """A = mean(T1) - mean(T0) and B = mean(T1^2) - mean(T0^2) over the year axis, after the cap;
    the growth-rate difference for any coefficient draw is then dg = b1 A + b2 B"""
    if cap is not None:
        T0, T1 = np.minimum(T0, cap), np.minimum(T1, cap)
    return T1.mean(axis) - T0.mean(axis), (T1 ** 2).mean(axis) - (T0 ** 2).mean(axis)


def dg(coef, spec, A, B, poor=None):
    """growth-rate difference for every draw: (draws, *A.shape); `poor` is boolean, broadcastable to A"""
    b1, b2, p1, p2 = coefs(coef, spec)
    sh = (-1,) + (1,) * A.ndim
    out = b1.reshape(sh) * A + b2.reshape(sh) * B
    if p1 is not None:
        out = np.where(poor, p1.reshape(sh) * A + p2.reshape(sh) * B, out)
    return out


def baseline():
    """per ISO3: meantemp, the UDel 1980-2010 population-weighted mean Burke's projections start from (his sample
    rule: years with both temperature and growth), and tconv, his CMIP5 RCP8.5 ratio of country to global warming"""
    d = panel()
    b = d[(d.year >= 1980) & d.UDel_temp_popweight.notna() & d.growthWDI.notna()]
    mt = b.groupby('iso').UDel_temp_popweight.agg(['mean', 'size'])
    mt.index = [FIX.get(i, i) for i in mt.index]
    mt = mt.sort_values('size').groupby(level=0).last()           # ZAR/COD and ROM/ROU both occur: the longer record
    tc = pd.read_csv(BUR + 'CCprojections/CountryTempChange_RCP85.csv')
    tc['iso'] = tc.GMI_CNTRY.map(lambda c: FIX.get(c, c))
    return (mt.rename(columns={'mean': 'meantemp', 'size': 'n_years'})
              .join(tc.groupby('iso').Tconv.mean().rename('tconv')))


if __name__ == '__main__':
    c = fit('pooled', ndraws=0)
    assert c.attrs == {'n_obs': 6584, 'n_countries': 166}, c.attrs
    assert np.allclose(c.loc[0, ['temp', 'temp2']], [0.0127177399, -0.000487076481], rtol=1e-6), c.iloc[0]
    assert tuple(np.round(c.loc[0, ['temp', 'temp2']], 4)) == PUBLISHED
    c0 = fit('lag5', ndraws=0, nlag=0)                     # the lag machinery with no lags is the pooled regression
    assert c0.attrs == c.attrs and np.allclose(c0[['temp_cum', 'temp2_cum']].values, c[['temp', 'temp2']].values)
    b = baseline()
    pub = pd.read_csv(dir_list['raw'] + 'burke2015/GDPcap_NOClimateChange_RCP85_SSP5.csv')
    pub.index = [FIX.get(i, i) for i in pub.ISO3]
    assert (b.meantemp - pub.meantemp).abs().max() < 1e-9, 'baseline temperatures differ from the published file'
    for spec in SPECS:
        c = fit(spec, ndraws=100)
        b1, b2, p1, p2 = coefs(c, spec)
        ok = b2[1:] < 0
        print(f'{spec:9s} n={c.attrs["n_obs"]}: temp {b1[0]:+.6f} temp2 {b2[0]:+.7f} optimum {-b1[0] / (2 * b2[0]):.1f} C '
              f'(bootstrap 5-95 %: {np.percentile(-b1[1:][ok] / (2 * b2[1:][ok]), 5):.1f}-'
              f'{np.percentile(-b1[1:][ok] / (2 * b2[1:][ok]), 95):.1f} C)'
              + (f'; poor: temp {p1[0]:+.6f} temp2 {p2[0]:+.7f}' if p1 is not None else ''))
    print('burke.py self-check passed')
