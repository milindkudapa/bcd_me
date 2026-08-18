# Path independence

Does the local climate at a given global warming level depend on the pathway taken to get there?
Tested pixel-by-pixel with a block-permutation KS test on monthly data, across three kinds of
pathway difference: aerosol cleanup (RAMIP), emissions scenario (ScenarioMIP), and direction of
travel (overshoot).

Data products and their schemas are documented separately in [DATA.md](DATA.md).

## Layout

```
ks_tools.py     shared test kernels, GWL window pickers, land_fraction, Wilks FDR  (self-checking)
DATA.md         what every output file contains, and its caveats
pipeline/       the current analysis chain
method/         the test itself: construction and validation
earlier/        superseded and prior-phase analyses, kept for comparison
reference/      the advisor's sample code the KS kernel came from
run/            slurm submission scripts
logs/           slurm output (gitignored)
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
variable (0.054), and for East-Asia-only cleanup (0.043–0.049) — whose GMST barely moves (1.0 yr
lead vs 5.6 yr for global cleanup), and whose local fingerprint over East Asia (~1.6× the global
rejection frequency) still fails Wilks (2016) FDR field significance. Under that FDR test the
aerosol and scenario families collapse to ≤0.001 of land with no floor argument needed.

Overshoot is the exception: the declining branch of `ssp534-over` differs from the rising branch
at the same warming level on 9–18% of land in the models with real ensembles, survives the FDR
test (0.077, driven by MRI-ESM2-0 and GISS at 2.0 K), and rests on few members — the one worth
pushing further. See the caveats in [DATA.md](DATA.md#5-caveats-that-change-conclusions).
