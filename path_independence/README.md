# Path independence

Does the local climate at a given global warming level depend on the pathway taken to get there?
Tested pixel-by-pixel with a block-permutation KS test on monthly data, across three kinds of
pathway difference: aerosol cleanup (RAMIP), emissions scenario (ScenarioMIP), and direction of
travel (overshoot).

Data products and their schemas are documented separately in [DATA.md](DATA.md).

## Layout

```
ks_tools.py     shared test kernels, GWL window pickers, land_fraction, Wilks FDR, 1° pipeline,
                derived variables (dtr) and seasonal windows  (self-checking)
DATA.md         what every output file contains, and its caveats
TECHNICAL.md    the technical report: summary, every experiment with its results, reproduction
AUDIT.md        validation of the aerosol result: power, effect sizes, pooled tests
pipeline/       the current analysis chain (incl. the country-level GDP response)
regional/       the regional cleanup arms: vs ssp585, all seven arms with the in-box test, by season
gdp/            the Burke GDP response revisited: his regressions and bootstrap, observation-anchored levels, the pathway experiment
method/         the test itself: construction and validation; audit/ holds the audit scripts and outputs
earlier/        superseded and prior-phase analyses, kept for comparison
reference/      the advisor's sample code the KS kernel came from
run/            slurm submission scripts, the download of ssp585 tasmax/tasmin, the folder-layout scripts
figures/        the four figures the reports embed (copies; the full set is in summer_data/figures/)
logs/           slurm output and executed notebook copies (gitignored)
dir_list.csv    symlink into ../code/ — get_params() reads it from the cwd
```

## Running

Everything goes through slurm; **nothing runs on the login node**, including short probes.

```bash
sbatch run/run_ks_blockperm.sbatch pipeline/ramip_summary.ipynb
sbatch run/run_ks_blockperm.sbatch ramip_summary.ipynb     # bare names are resolved too
sbatch run/run_ks_1deg.sbatch                              # the 1° regrid + KS notebook (needs ESMF, below)
python ks_tools.py                                          # self-check (submit it, don't run it here)
```

xESMF needs ESMF/esmpy, which is conda-only and needs a newer libstdc++ than the cluster's EL8:
`sbatch run/build_esmf_env.sbatch` builds a small nompi conda env at `summer/.esmf`, links
`esmpy` (and its dist-info) into the working venv, and verifies the import. `run_ks_1deg.sbatch`
sets the two variables that make it load (`ESMFMKFILE`, `LD_LIBRARY_PATH`).

Notebooks execute in place and cache everything they compute, so a rerun after a failure resumes
rather than restarts. Working directory stays at this folder regardless of which subfolder the
notebook lives in, because `get_params()` reads the `dir_list.csv` symlink here.

## pipeline/ — in dependency order

| notebook | produces |
|---|---|
| `ramip_data_layer.ipynb` | annual GMST, GWL crossings, AR6-regional means for every RAMIP experiment |
| `scenariomip_download.ipynb` | `raw/scenariomip/` — tas+pr for ssp126/245/585/534-over from pangeo |
| `ks_path_independence_ramip_gwl.ipynb` | `tas` KS maps, `ssp370` vs `ssp370-126aer` at matched GWL 2.0 |
| `ks_path_independence_ramip_gwl_vars.ipynb` | the same for `pr`, `tasmax`, `tasmin` |
| `ks_location_scenarios.ipynb` | SSP pairs, overshoot rising-vs-declining, RAMIP at 1.5–3.0 K |
| `ramip_impacts_gdp.ipynb` | Burke (2015) GDP-growth response to the aerosol pathway |
| `build_data_product.ipynb` | `aux/gwl_timeseries.zarr` — GMST, crossings and regional timeseries in one store |
| `ks_gwl_1deg.ipynb` | everything regridded to a common 1° grid (xESMF), `gdp` (Burke g(T)) as a fifth variable, the East-Asia-only cleanup arm, Wilks (2016) FDR field significance, and the zarr's `ks` group |
| `ramip_summary.ipynb` | summary tables 1–5 and figures 1–4 |
| `project_figure.ipynb`, `story_figures.ipynb` | the overview and four-figure project summaries |
| `ramip_impacts_gdp_country.ipynb` | the Burke response at country level: GHS-POP population weights, World Bank GDP, Δg, ΔGDP per capita and ΔGDP per country, country-series KS |

## regional/

| notebook | produces |
|---|---|
| `ks_regional_vs_ssp585.ipynb` | `ssp585` vs the East-Asia and global arms (1deg runner); the MIROC6 diagnosis |
| `ks_regional_arms.ipynb` | all seven arms vs own parent and vs ssp585, `tas`/`pr`/`gdp`/`tasmax`/`tasmin`/`dtr`, global and in-box, via `run/run_regional_arms_array.sbatch` |
| `ks_regional_arms_seasonal.ipynb` | the same by season with the pooled in-box test, via `run/run_regional_arms_seasonal_array.sbatch` |

## gdp/

| file | produces |
|---|---|
| `burke.py` | Burke, Hsiang & Miguel's regressions re-estimated from their replication package (`python gdp/burke.py` self-checks against the paper) |
| `country_tas.py` | annual population-weighted country temperature for every RAMIP run (array over models via `run_gdp.sbatch`) |
| `ramip_gdp_burke.ipynb` | matched GWL 2.0 with the exact coefficients, Burke's 1000-draw bootstrap, observation-anchored levels and Δg tested against its own null |
| `ramip_gdp_pathway.ipynb` | the same response year by year along the pathways, compounded to 2050 |

## method/

- `ks_block_permutation.ipynb` — the test's construction and its Monte-Carlo size verification
  under autocorrelation, a cyclostationary ENSO surrogate, real ensemble pairs, and iid.
- `ks_blockperm_diagnostics.ipynb` — why it behaves as it does: observed vs block-permuted ECDFs,
  autocorrelation preservation, permutation vs asymptotic nulls, and the **false-positive floor
  of 0.057** measured on 45 same-forcing member pairs. Every land fraction in this project is read
  against that number, not against zero.

## What it found

At a matched global warming level, the aerosol-cleanup pathway and the choice of SSP leave the
local monthly distribution indistinguishable from same-forcing ensemble noise — land rejection
fractions of 0.050–0.058 across `tas`, `pr`, `tasmax`, `tasmin` against a 0.057 floor. The
year-matched comparison that appeared to show large differences was mostly reading the GMST offset
from aerosol cleanup reaching a given level ~6.5 years earlier.

The same holds on a common 1° grid, for the Burke (2015) GDP-growth response `g(T)` as a fifth
variable (0.054), and for every regional cleanup arm (East Asia, South Asia, Africa, North
America/Europe, combined Asia, the carbonaceous-only arms; 0.039–0.065), against the own `ssp370`
parent and against `ssp585`. Under the FDR test the aerosol and scenario families sit at ≤0.005 of
land with no floor argument needed.

What the per-cell test cannot see, the audit (`AUDIT.md`) and the regional notebooks measure
directly: inside its own emission region a cleanup leaves a small, real fingerprint — East Asia
and South Asia +0.1–0.2 K and +2–4 % precipitation at matched 2 °C, mostly in daytime maxima and
the diurnal range, with a season (monsoon rain over South Asia in JJA, haze-season temperature
along the Indo-Gangetic plain in DJF, year-round over eastern China) — detectable once the
ensemble is pooled (40–50 % of the box rejects) or the box mean is tested, and worth about
−35 ± 22 bn US$ per year of growth in Burke's framework for global cleanup (his coefficient
uncertainty: −78 to +12 bn), because the cleaned world is warmer where people live; along the actual
pathway, where the cleanup also warms the globe faster, the same framework gives −3.4 % of GDP per
capita by 2050 (`gdp/`). One pair of 10-year windows needs a shift of ~1 K per cell to
see any of this.

Overshoot is the exception: the declining branch of `ssp534-over` differs from the rising branch
at the same warming level on 9–18% of land in the models with real ensembles, survives the FDR
test (0.078, driven by MRI-ESM2-0 and GISS at 2.0 K), and rests on few members — the one worth
pushing further. See the caveats in [DATA.md](DATA.md#5-caveats-that-change-conclusions).
