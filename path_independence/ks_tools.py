"""Shared machinery for the block-permutation KS path-independence tests.

Extracted verbatim from ks_path_independence_ramip_gwl.ipynb so the scenario/overshoot
notebooks stop copy-pasting it. `ks_statistic_sorted` is the advisor's kernel from
enso_bias_impacts/code/4-fig_impact_figures.ipynb; the block-permutation wrapper's size is
verified in ks_block_permutation.ipynb (0.057 on same-forcing ensemble pairs, nominal 0.05).
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import xarray as xr
from numba import njit, prange
from scipy import stats as sstats

BLOCK = 12          # months per block: one annual cycle
NDRAWS = 999
WIN = 10            # years per sample -> n = 120 months


@njit
def ks_statistic_sorted(x, y):
    """KS statistic for already-sorted arrays (verbatim from 4-fig_impact_figures.ipynb)."""
    nx = len(x); ny = len(y)
    i = 0; j = 0
    cdfx = 0.0; cdfy = 0.0
    d = 0.0
    while i < nx and j < ny:
        if x[i] <= y[j]:
            i += 1
            cdfx = i / nx
        else:
            j += 1
            cdfy = j / ny
        diff = abs(cdfx - cdfy)
        if diff > d:
            d = diff
    while i < nx:
        i += 1
        cdfx = i / nx
        diff = abs(cdfx - cdfy)
        if diff > d:
            d = diff
    while j < ny:
        j += 1
        cdfy = j / ny
        diff = abs(cdfx - cdfy)
        if diff > d:
            d = diff
    return d


@njit
def ks_block_permutation_test(x, y, block=BLOCK, ndraws=NDRAWS):
    """Block-permutation KS p: keep D, resample the null by shuffling whole annual blocks."""
    nx = len(x); ny = len(y)
    nbx = nx // block; nby = ny // block
    nb = nbx + nby
    D_obs = ks_statistic_sorted(np.sort(x.copy()), np.sort(y.copy()))
    pool = np.empty(nx + ny, dtype=x.dtype)
    pool[:nx] = x
    pool[nx:] = y
    order = np.arange(nb)
    xperm = np.empty(nx, dtype=x.dtype)
    yperm = np.empty(ny, dtype=x.dtype)
    exceed = 0
    for draw in range(ndraws):
        for i in range(nb - 1, 0, -1):
            j = np.random.randint(i + 1)
            tmp = order[i]; order[i] = order[j]; order[j] = tmp
        for b in range(nbx):
            s = order[b] * block
            xperm[b * block:(b + 1) * block] = pool[s:s + block]
        for b in range(nby):
            s = order[nbx + b] * block
            yperm[b * block:(b + 1) * block] = pool[s:s + block]
        xperm.sort(); yperm.sort()
        if ks_statistic_sorted(xperm, yperm) >= D_obs:
            exceed += 1
    return (exceed + 1) / (ndraws + 1)   # Phipson & Smyth 2010


@njit(parallel=True)
def block_perm_map(A, B, block=BLOCK, ndraws=NDRAWS):
    """p-value per pixel; A, B are (pixels, time)."""
    P = A.shape[0]
    out = np.empty(P)
    for px in prange(P):
        out[px] = ks_block_permutation_test(A[px], B[px], block, ndraws)
    return out


def ks_fast(a1, a2):
    """asymptotic p per pixel, for comparison only (kstwo, matching the authors' convention)"""
    P, n1 = a1.shape; n2 = a2.shape[1]
    z = np.concatenate([a1, a2], axis=1)
    order = np.argsort(z, axis=1, kind='stable')
    zs = np.take_along_axis(z, order, axis=1)
    ind = np.where(order < n1, 1.0 / n1, -1.0 / n2)
    cd = np.cumsum(ind, axis=1)
    valid = np.ones_like(cd, dtype=bool)
    valid[:, :-1] = zs[:, 1:] != zs[:, :-1]
    D = np.abs(np.where(valid, cd, 0.0)).max(axis=1)
    en = int(np.round(n1 * n2 / (n1 + n2)))
    uD, inv = np.unique(D, return_inverse=True)
    return np.clip(sstats.kstwo.sf(uD, en), 0, 1)[inv]


def member_files(root, mod, exp, var='tas'):
    """run -> file under {root}/{mod}/. Runs split across time spans: keep the longest."""
    out = {}
    for f in sorted(glob.glob(f'{root}{mod}/{var}_Amon_{mod}_{exp}_r*.zarr')):
        run = re.search(r'_(r\d+i\d+p\d+f\d+)_', f).group(1)
        end = re.search(r'-(\d{4})\d{4}\.zarr$', f)
        end = int(end.group(1)) if end else 0
        if run not in out or end > out[run][1]:
            out[run] = (f, end)
    return {r: v[0] for r, v in out.items()}


def experiments(root, mod, var='tas'):
    """experiment ids present on disk for this model/variable"""
    pat = re.compile(rf'{var}_Amon_{re.escape(mod)}_(.+?)_r\d+i\d+p\d+f\d+_')
    return sorted({m.group(1) for f in glob.glob(f'{root}{mod}/{var}_Amon_{mod}_*.zarr')
                   if (m := pat.search(os.path.basename(f)))})


def gmst(da):
    return da.weighted(np.cos(np.deg2rad(da.lat))).mean([d for d in da.dims if d != 'time'])


def window_vals(fn, var, y_end, win=WIN):
    """monthly values over the `win` years ending at y_end, as (pixels, months), lat, lon, n"""
    ds = xr.open_zarr(fn, consolidated=False)
    yr = ds.time.dt.year
    da = ds[var].sel(time=(yr > y_end - win) & (yr <= y_end))
    v = da.load().transpose('lat', 'lon', 'time')
    return v.values.reshape(-1, v.sizes['time']), v.lat, v.lon, v.sizes['time']


def gwl_window(years, vals, gwl, win=WIN, which='first'):
    """Last year of the first (or last) `win`-yr running-mean window reaching `gwl`.

    `which='last'` is what an overshoot pathway needs: the final window still above the level,
    i.e. the down-crossing side of the hump. NaN if the level is never reached.
    """
    roll = pd.Series(vals, index=years).rolling(win).mean()
    hit = roll[roll >= gwl]
    if not len(hit):
        return np.nan
    return int(hit.index[0] if which == 'first' else hit.index[-1])


def window_at_level(years, vals, gwl, tol=0.1, win=WIN, side='rising'):
    """End year of the `win`-yr window whose running-mean anomaly sits closest to `gwl`.

    `gwl_window` above answers "when does this run first reach the level", which is the right
    question on a monotonic pathway but the wrong one on an overshoot: after the peak the last
    window at-or-above the level can be the end of the record at a much higher temperature. This
    picks the window closest to the target, restricted to the rising or declining branch, and
    returns None when even the best window misses by more than `tol` — so a run that never comes
    back down is dropped rather than silently compared at the wrong warming level.
    """
    roll = pd.Series(vals, index=years).rolling(win).mean().dropna()
    if roll.empty:
        return None
    peak_yr = roll.idxmax()
    branch = roll[roll.index <= peak_yr] if side == 'rising' else roll[roll.index >= peak_yr]
    if branch.empty:
        return None
    best = (branch - gwl).abs().idxmin()
    return int(best) if abs(branch[best] - gwl) <= tol else None


def pair_ks(fn_a, fn_b, var, end_a, end_b, win=WIN, block=BLOCK, ndraws=NDRAWS):
    """block-perm + asymptotic p maps for one pair of windows. None if a window is short."""
    va, lat, lon, n1 = window_vals(fn_a, var, end_a, win)
    vb, *_, n2 = window_vals(fn_b, var, end_b, win)
    if n1 != 12 * win or n2 != 12 * win:
        return None
    return (block_perm_map(va, vb, block, ndraws).reshape(len(lat), len(lon)),
            ks_fast(va, vb).reshape(len(lat), len(lon)), lat, lon)


_LAND = {}


def land_mask(lat, lon):
    """Boolean land mask as a plain numpy array, cached on the exact coordinate bytes.

    Two traps, both of which silently produce a land fraction of exactly 0:
    - keying the cache on grid *shape* hands one model another model's mask;
    - keying it on rounded coordinates does the same for grids that differ in the last few
      float digits (CNRM's and MIROC6's Gaussian lats agree to 6 dp but not exactly).
    Returning numpy rather than xarray also stops `.where()` from aligning on coordinates and
    quietly emptying the array.
    """
    key = (np.asarray(lat).tobytes(), np.asarray(lon).tobytes())
    if key not in _LAND:
        import regionmask
        _LAND[key] = regionmask.defined_regions.natural_earth_v5_0_0.land_110.mask(
            lon, lat).notnull().values
    return _LAND[key]


def land_fraction(ds, dim, alpha=0.05, var='p_blockperm'):
    """Mean over `dim` of the land fraction of pixels with p < alpha."""
    land = land_mask(ds.lat.values, ds.lon.values)
    p = ds[var].transpose(dim, 'lat', 'lon').values
    return float(((p < alpha) & land).sum(axis=(1, 2)).mean() / land.sum())


def demo():
    """iid sanity check: both tests should be roughly uniform, i.e. not reject"""
    r = np.random.default_rng(0)
    x, y = r.standard_normal((4, 120)), r.standard_normal((4, 120))
    p_perm = block_perm_map(x, y, BLOCK, NDRAWS)
    p_asym = ks_fast(x, y)
    assert (p_perm > 0.05).all(), p_perm
    assert (p_asym > 0.05).all(), p_asym
    # a shifted sample must be rejected by both
    p_shift = block_perm_map(x, y + 3.0, BLOCK, NDRAWS)
    assert (p_shift < 0.05).all(), p_shift
    # gwl_window picks the up- and down-crossing sides of an overshoot hump
    yrs = np.arange(2015, 2101)
    hump = 2.5 - 1.5 * ((yrs - 2060) / 40.0) ** 2
    assert gwl_window(yrs, hump, 2.0, which='first') < gwl_window(yrs, hump, 2.0, which='last')
    assert np.isnan(gwl_window(yrs, hump, 3.0, which='first'))
    # window_at_level must land ON the level on both branches, and refuse when the run never
    # returns to it (the failure mode that made the first overshoot test meaningless)
    up = window_at_level(yrs, hump, 2.0, side='rising')
    dn = window_at_level(yrs, hump, 2.0, side='declining')
    roll = pd.Series(hump, index=yrs).rolling(WIN).mean()
    assert up < dn and abs(roll[up] - 2.0) <= 0.1 and abs(roll[dn] - 2.0) <= 0.1
    truncated = hump[yrs <= 2070]                      # cut before it comes back down
    assert window_at_level(yrs[yrs <= 2070], truncated, 1.5, side='declining') is None
    # two grids that differ only in the last float digits must not share a mask, and neither
    # may come back as a land fraction of exactly 0 (the real bug this guards)
    lat = np.linspace(-88.9277353522, 88.9277353522, 128)
    lon = np.linspace(0, 358.59375, 256)
    rng = np.random.default_rng(1)
    mk = lambda la: xr.Dataset(
        {'p_blockperm': (('run', 'lat', 'lon'), rng.uniform(size=(3, 128, 256)))},
        coords={'run': [0, 1, 2], 'lat': la, 'lon': lon})
    a, b = mk(lat), mk(lat + 1e-9)
    fa, fb = land_fraction(a, 'run'), land_fraction(b, 'run')
    assert 0.02 < fa < 0.08 and 0.02 < fb < 0.08, (fa, fb)
    assert len(_LAND) == 2, 'near-identical grids collided in the mask cache'
    print('ks_tools demo OK')


if __name__ == '__main__':
    demo()


def sig_fdr(ps, FDR=0.2):
    ''' Calculate field significance contingent on a given FDR

    Adapted from Wilks 2016: "'The Stippling Shows Statistically
    Significant Grid Points': How Research Results are Routinely
    Overstated and Overinterpreted, and What to Do about it"
    (implementation from ks905383/random funcs.py)

    Parameters
    -----------------------
    ps : xr.DataArray
        a DataArray of significance values over some grid containing
        lat / lon, and any optional number of other dimensions
    FDR : float (default 0.2)
        the desired false discovery rate (FDR)

    Returns
    -----------------------
    sigTests : xr.DataArray
        a DataArray of booleans showing whether a pixel is significant
        based on a given FDR
    '''

    if type(ps) != xr.core.dataarray.DataArray:
        raise TypeError('`ps` must be an xarray DataArray.')

    # stack geographic variables into one
    ps_stack = ps.stack(loc=('lat','lon'))

    # calculate (i/N)*a_fdr, where i is the rank, N is the
    # number of pixels (highest rank gets it; otherwise would
    # be a count of non-nans, which might be faster but uglier),
    # and a_fdr is the significance marker
    sigLevel = (ps_stack.rank('loc')/ps_stack.rank('loc').max('loc'))*FDR

    # significant pixels are those where
    # p < p_fdr = max{p(i):p(i) < (i/N)a_fdr}
    sigTests = ps_stack < ps_stack.where(ps_stack < sigLevel).max('loc')

    # restore nans, which are removed through the <
    sigTests = sigTests.where(~np.isnan(ps_stack))

    # unstack and return
    sigTests = sigTests.unstack()

    return sigTests


def fdr_land_fraction(ds, dim, fdr=0.2, var='p_blockperm'):
    """Mean over `dim` of the land fraction of FDR-significant pixels (Wilks 2016).

    Ocean is NaN'd before `sig_fdr` so N in the FDR ranking is the number of land
    pixels — the field being interrogated is land, same as `land_fraction`.
    """
    land = land_mask(ds.lat.values, ds.lon.values)
    p = ds[var].where(xr.DataArray(land, coords={'lat': ds.lat, 'lon': ds.lon}))
    sig = sig_fdr(p, FDR=fdr)
    return float((sig == 1).sum(('lat', 'lon')).mean(dim) / land.sum())
