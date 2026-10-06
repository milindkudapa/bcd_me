import os, sys, io, contextlib, re, numpy as np, pandas as pd, xarray as xr
sys.path.insert(0, '/burg-archive/home/mck2199/summer/bcd_me/code')
sys.path.insert(0, '/burg-archive/home/mck2199/summer/bcd_me/path_independence')
os.chdir('/burg-archive/home/mck2199/summer/bcd_me/path_independence')   # get_params reads ./dir_list.csv
import ks_tools as kt
from funcs_support import get_params
dir_list = get_params()
RAMIP = dir_list['raw'] + 'ramip/'
AUX = dir_list['aux']
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out') + '/'
PRODUCT = AUX + 'gwl_timeseries.zarr'
GMST = xr.open_zarr(PRODUCT, group='gmst'); GWLDS = xr.open_zarr(PRODUCT, group='gwl')
GWL = 2.0

def gwl_ends(mod, exp):
    ids = GMST.run_id.sel(model=mod, exp=exp).values
    we = GWLDS.window_end.sel(model=mod, exp=exp, gwl=GWL).values
    return {str(r): int(e) for r, e in zip(ids, we) if r and np.isfinite(e)}

def member_idx(mod, exp, run):
    ids = GMST.run_id.sel(model=mod, exp=exp).values
    return int(np.nonzero(ids == run)[0][0])

BOXES = {'eas': [((20, 53), (95, 133))], 'sas': [((5, 35), (65, 95))], 'afr': [((-35, 35), (-20, 60))],
         'nae': [((35, 70), (-20, 45)), ((25, 70), (-150, -45))], 'asia': [((20, 53), (95, 133)), ((5, 35), (65, 95))],
         'safca': [((-35, 12), (-20, 50))], 'sasca': [((5, 35), (65, 95))]}

def box_mask(lat, lon, boxes):
    lon = np.asarray(lon)
    LA, LO = np.meshgrid(lat, np.where(lon > 180, lon - 360, lon), indexing='ij')
    m = np.zeros(LA.shape, bool)
    for (a, b), (c, d) in boxes:
        m |= (LA >= a) & (LA <= b) & (LO >= c) & (LO <= d)
    return m

def notebook_setup():
    """run cell 1 of regional/ks_regional_arms.ipynb (JOBS, ARMS, PARENT, parent_pairs, ...) into a dict"""
    import nbformat
    nb = nbformat.read('/burg-archive/home/mck2199/summer/bcd_me/path_independence/regional/ks_regional_arms.ipynb', 4)
    g = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(nb.cells[1].source, g)
    return g
