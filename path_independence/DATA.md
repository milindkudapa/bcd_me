# Data products

Everything the notebooks in this directory write, what is in it, and what to distrust.

Paths are relative to `summer_data/` (the `dir_list` roots from `code/funcs_support.py`:
`raw/` and `aux/`). All of it is derived — deleting any of it costs compute, not information.

---

## 1. `aux/gwl_timeseries.zarr` — the product

One store, four groups, both archives (RAMIP and ScenarioMIP/overshoot) on a shared set of
coordinates. Written by `pipeline/build_data_product.ipynb`; the `ks` group by
`pipeline/ks_gwl_1deg.ipynb`.

```python
import xarray as xr
P = 'summer_data/aux/gwl_timeseries.zarr'
gmst = xr.open_zarr(P, group='gmst')
gwl  = xr.open_zarr(P, group='gwl')
reg  = xr.open_zarr(P, group='regional')
ks   = xr.open_zarr(P, group='ks')
```

| group | variables | dims |
|---|---|---|
| `gmst` | `gmst_anom`, `run_id`, `collection` | model 10 × exp 15 × member 10 × year 86 |
| `gwl` | `window_end`, `down_end`, `peak` | model 10 × exp 15 × member 10 × gwl 4 |
| `regional` | `tas`, `pr`, `tasmax`, `tasmin` | model 10 × exp 15 × member 10 × year 86 × region 47 |
| `ks` | `p_blockperm`, `run_id` | comparison 2 × var 5 × model 8 × pair 10 × lat 180 × lon 360 |

The `ks` group is the pixel-wise block-permutation p-maps at matched GWL 2.0 on a **common 1°
grid** (xESMF, bilinear for temperature-derived variables, conservative for `pr`, transforms
applied on the native grid *before* regridding). `comparison` is `126aer` (global cleanup) or
`EAS126aer` (East-Asia-only); `var` includes `gdp`, the Burke (2015) growth response
`g(T) = 0.0127T − 0.0005T²` applied to monthly `tas` in °C. `pair` follows the same
index-not-label convention as `member`, with `run_id(model, pair)` carrying the identifier.

**Coordinates**

- `model` — 9 RAMIP models plus `CanESM5` (the ScenarioMIP stand-in for `CanESM5-1`).
- `exp` — 15 experiments on one axis. `collection(exp)` says which archive each came from:
  - *ramip*: `ssp370`, `ssp370-126aer`, `ssp370-LE` (CESM2's reference), `ssp370-ramip`
    (NorESM only), the regional-cleanup arms `ssp370-{AFR,ASIA,EAS,NAE,SAS}126aer`, and
    `ssp370-{SAF,SAS}126ca`.
  - *scenariomip*: `ssp126`, `ssp245`, `ssp585`, `ssp534-over`.
- `member` — **an index, not a label.** See the warning below.
- `year` — 2015–2100.
- `gwl` — 1.5, 2.0, 2.5, 3.0 K above 1850–1900.
- `region` — AR6 land region numbers 0–45, plus `-1` for a global-land aggregate.
  `abbrevs(region)` and `names(region)` come along.

**Variables**

- `gmst_anom` [K] — annual area-weighted global-mean surface air temperature, minus the model's
  own 1850–1900 mean from CMIP6 `historical`.
- `run_id(model, exp, member)` — the CMIP6 member identifier the member index stands for; empty
  string where the combination does not exist.
- `window_end` [year] — last year of the **first** 10-yr running-mean window reaching the GWL.
  The sample it denotes is `window_end-9 … window_end` inclusive. NaN if the run never reaches
  the level within its record.
- `down_end` [year] — last year of the 10-yr window **closest to** the GWL on the *declining*
  branch, required to be within 0.1 K of it. Meaningful only for `ssp534-over`. NaN when the run
  never comes back down — which is most of them, and is the point (see §5).
- `peak` [K] — maximum 10-yr running-mean GMST anomaly of the run.
- `tas`/`tasmax`/`tasmin` [K], `pr` [kg m-2 s-1] — annual means, cos-lat-weighted over each AR6
  land region.

### Read this before slicing

**`member` is an index.** Run identifiers are not comparable across models (`r101i1p1f1` in MRI,
`r1i1p3f1` in GISS), so the dimension is `0…9` and the identifier lives in `run_id`. Members are
ordered **lexically**, so `member=0` is `r10i1p1f1` for MIROC6, *not* `r1i1p1f1`. Selecting the
same member index across models compares unrelated runs.

```python
# find the index for a specific run
i = gmst.run_id.sel(model='MIROC6', exp='ssp370').values.tolist().index('r5i1p1f1')
```

**The array is sparse.** No model runs every experiment; 36.6% of the `gmst` grid is populated
and the rest is NaN. That is structure, not corruption. Count with `.notnull()`, and check
`run_id != ''` before trusting a slice.

**Known substitutions**, also recorded in the store's attrs:

- `CanESM5-1` has no `historical` on pangeo, so its 1850–1900 baseline is borrowed from its
  parent `CanESM5` — a 0.1 K baseline error moves a crossing year by a few years.
- `CanESM5` (not `CanESM5-1`) is what the ScenarioMIP collection contains.
- `EC-Earth3-AerChem` has **no ScenarioMIP data at all** — it is an AerChemMIP-only `source_id`.
- ScenarioMIP `tasmax`/`tasmin` were never downloaded, so those variables are RAMIP-only.

### Examples

```python
# GWL2.0 decade of every ssp370-126aer run
gwl.window_end.sel(exp='ssp370-126aer', gwl=2.0).dropna('member', how='all')

# the aerosol lead: 126aer reaches 2 K ~14 yr before ssp370 in MIROC6
gwl.window_end.sel(model='MIROC6', gwl=2.0, exp=['ssp370', 'ssp370-126aer'])

# South Asia annual tas, both global arms
sas = int(reg.abbrevs.values.tolist().index('SAS'))
reg.tas.sel(model='MIROC6', region=reg.region.values[sas], exp=['ssp370', 'ssp370-126aer'])
```

---

## 2. `aux/data_layer/` — the tidy intermediates the store is built from

| file | contents |
|---|---|
| `gmst_annual_long.csv` | RAMIP GMST, long format: model, exp, run, year, gmst_anom (31,787 rows) |
| `scenariomip_gmst_annual.csv` | same for the ScenarioMIP collection |
| `gwl_crossings.csv`, `gwl_crossing_summary.csv` | RAMIP crossing years per run per GWL |
| `scenariomip_gwl_crossings.csv` | same, plus `up_end`/`down_end`/`peak` |
| `scenariomip_availability.csv`, `scenariomip_inventory.csv` | what pangeo offered vs what landed on disk |
| `regional/{var}_{model}.nc` | RAMIP AR6-regional annual means, 36 files, dims exp×run×year×region |
| `regional_scenariomip/{var}_{model}.nc` | same for ScenarioMIP, 16 files (tas + pr only) |

Prefer the zarr; these are kept because they are what the store is verified against.

---

## 3. `aux/diagnostics/` — KS test output

Pixel-wise p-value maps. The native-grid files below are per-model; their regridded 1°
counterparts live in the zarr's `ks` group (§1) and in `ramip_gwl/ks_1deg/`.

Every file has `p_blockperm` and `p_asymptotic`. The block-permutation p is the one to use — its
size is verified; the asymptotic p is conservative on seasonal-cycle-dominated monthly land data
and is kept only for comparison.

| path | what | dims |
|---|---|---|
| `ramip_gwl/ks/{model}.nc` | `ssp370` vs `ssp370-126aer`, `tas`, matched at GWL2.0 | run×lat×lon |
| `ramip_gwl/ks_{pr,tasmax,tasmin}/{model}.nc` | same for the other variables | run×lat×lon |
| `ramip_gwl/ks_1deg/{var}_{model}.nc` | all five variables (incl. `gdp`) regridded to 1° | run×lat×lon |
| `ramip_gwl/ks_1deg/eas_{var}_{model}.nc` | same vs `ssp370-EAS126aer` | run×lat×lon |
| `ramip_gwl/land_fractions_1deg*.csv` | raw p<0.05 land fractions, 1° vs native, global and EAS arms |  |
| `ramip_gwl/land_fractions_1deg_fdr*.csv` | Wilks (2016) FDR-significant land fractions (α_FDR = 0.2) |  |
| `ramip_gwl/wilks_fdr_families.csv` | raw vs FDR for all four families on one scale |  |
| `ks_scenarios/{scen,over,ramip}_*_{var}.nc` | SSP pairs, overshoot, RAMIP at 4 GWLs (104 files) | pair×lat×lon |
| `ks_blockperm_enspairs/CanESM5-1.nc` | 45 same-forcing member pairs — the null calibration | pair×lat×lon |
| `ramip_gwl/baseline_1850_1900.csv` | per-model pre-industrial GMST and its member spread |  |
| `ramip_gwl/gwl2.0_windows.csv` | the windows the GWL2.0 analysis actually used |  |
| `impacts/burke_gwl2.0.nc` | Burke GDP-growth response: `t_anchor`, `t_other`, `g_*`, `dt`, `dg` | model×run×lat×lon |

Test parameters, identical everywhere: 10-yr windows (n = 120 months), 12-month blocks, 999
permutation draws, α = 0.05.

**The number that makes these interpretable:** the same test on 45 *same-forcing* CanESM5-1
member pairs rejects on **0.057** of land. That is this test's false-positive floor. A land
fraction at 0.057 means "indistinguishable", not "identical", and zero is not the null
expectation. The Wilks (2016) FDR columns need no floor: `ks_tools.sig_fdr` controls the false
discovery rate over the land field, so ~0 *is* the null expectation there — the aerosol and
scenario families sit at ≤0.001, overshoot at 0.077.

`ks_blockperm_ramip/`, `ks_ramip_multi/`, `ks_empirical_null/`, `ramip_regional_means/` and the
loose `ks_pvals_*.nc` are outputs of earlier notebooks, kept for comparison against the
year-matched analysis.

---

## 4. `aux/summary/` — reporting tables

`table1_ramip_members.csv`, `table1b_scenariomip_members.csv` (what data exists),
`table2_gwl_crossings.csv` (crossing years and the aerosol lead),
`table3_landfrac_by_var.csv` (rejection fractions by variable),
`table4_families.csv` (all three comparison families on one scale),
`table5_regional_response.csv` (regional ΔT and Δpr at matched GWL).

Figures land in `summer_data/figures/` as `summary_fig{1..4}_*.png` plus the per-analysis ones.

---

## 5. Caveats that change conclusions

- **Overshoot windows.** `ssp534-over` stores begin at the **2040 branch year**, and most runs
  never return below a given level by 2100. Taking "the last window at or above the level" as the
  down-crossing therefore returns the end of the record at a much higher temperature — CanESM5's
  nominal 1.5 K window sits at 3.5 K, which produces a spurious rejection of 0.73. `down_end`
  guards against this by requiring the window to be within 0.1 K of the target on the declining
  branch, and 17 of 24 model×level combinations legitimately drop as a result. The rising branch
  in the overshoot comparison comes from the same model's `ssp585`/`ssp245`, not from the same
  run.
- **Variant mixing in the overshoot pairing.** Because the two sides come from different
  experiments, member ids do not correspond, and pairing them by index silently compared
  different model configurations: GISS-E2-1-G submits `p1`/`p3`/`p5` (different aerosol
  chemistry) under `f1`/`f2` (different forcing datasets). Six of fifteen pairs were mismatched,
  which roughly doubled GISS's apparent hysteresis (0.179 → **0.096** at 2.0 K once fixed).
  Pairs are now formed within a `p<n>f<n>` variant and the notebook asserts on mismatch. Any new
  cross-experiment comparison needs the same check; same-run comparisons (scenario pairs, RAMIP)
  are immune by construction.
- **RAMIP MIROC6 is not CMIP6 MIROC6.** The deela RAMIP `ssp370` members `r1`–`r3` share ids with
  the CMIP6 ScenarioMIP runs on pangeo but are different simulations (monthly fields differ by
  3.4 K rms, GMST series correlate 0.988), and different *configurations*: the 86-yr climatology
  differs systematically — Arctic (70–90N) **+5 K**, Southern Ocean **−1.2 K**, 1.4 K rms — with
  the same pattern in r1 and r2, while realization noise (r1 vs r2 within one archive) is ~0.
  The deela files also carry a 360-day calendar (MIROC6 is Gregorian). Consequence: any
  RAMIP-vs-ScenarioMIP comparison must exclude MIROC6 (it rejected 52% of land against `ssp585`
  at matched GWL while six other models sat at the floor). Within-archive comparisons — RAMIP
  arms against RAMIP `ssp370`, ScenarioMIP pairs, overshoot — are unaffected. Other RAMIP models
  may be re-runs too (MRI's `r101…` ids are), so cross-archive pairing is always by variant and
  always read against this caveat.
- **Land masks.** Caching a land mask by grid *shape* — or by coordinates rounded to a few
  decimals — hands one model another model's mask, after which `xr.where()` aligns strictly,
  empties the array, and the land fraction reads exactly **0.000** with no error. Use
  `ks_tools.land_fraction`, which keys on `lat.tobytes()` and masks in numpy. `python
  ks_tools.py` asserts two grids differing by 1e-9 stay separate.
- **q01-across-runs maps mislead.** At N = 6–10 runs `quantile(0.01, dim='run')` is essentially
  the minimum, and the minimum of 10 uniform p-values averages ~0.09 under a true null, so those
  maps look dark everywhere. Report land fractions, or the per-pixel fraction of runs rejecting
  (fixed null expectation 0.05 at any N).
- **Burke impacts are differences, not levels.** Country-level coefficients applied per pixel,
  no population weighting, no bias correction, no reanalysis-uncertainty sampling — which
  Schwarzwald et al. find is the *largest* term in absolute damages. The arm-to-arm difference is
  meaningful; the damage level is not.
- **Sample sizes are small in places.** CNRM-ESM2-1 and EC-Earth3-AerChem carry 6 matched runs,
  CESM2's ScenarioMIP entry has 3 members, several overshoot comparisons have n = 1, and the
  multi-GWL RAMIP comparisons are capped at 5 runs (the cap is printed, never silent). Check
  `table1*` before quoting a multi-model number.

---

## 6. Rebuilding

Everything caches; delete a cache to force recomputation. Jobs go through slurm — **never the
login node**:

```bash
sbatch run/run_ks_blockperm.sbatch pipeline/<notebook>.ipynb
```

Order from scratch: `ramip_data_layer` → `scenariomip_download` →
`ks_path_independence_ramip_gwl` → `ks_path_independence_ramip_gwl_vars` →
`ks_location_scenarios` → `ramip_impacts_gdp` → `build_data_product` → `ks_gwl_1deg`
(via `run/run_ks_1deg.sbatch`; needs the ESMF env from `run/build_esmf_env.sbatch` once) →
`ramip_summary`.

Two environment notes: burg needs `export SSL_CERT_FILE=$(python -c 'import certifi;
print(certifi.where())')` for Natural Earth downloads (already in the sbatch), and writing a
pangeo-sourced dataset to zarr needs its inherited encoding cleared —
`for v in ds.variables: ds[v].encoding = {}` — or zarr v3 raises
`Expected a BytesBytesCodec. Got numcodecs.blosc.Blosc`.
