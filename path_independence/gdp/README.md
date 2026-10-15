# gdp/ — the Burke GDP response, revisited

A self-contained redo of the economic layer. Nothing here changes the pipeline notebooks; it reads their inputs
(the RAMIP archive, the GWL windows, the population and World Bank layers) and Burke, Hsiang and Miguel's own
replication package (`raw/burke2015/`), and writes to `diagnostics/gdp/` and `figures/e8_gdp_burke/`.

What it fixes, found by `method/audit/gdp_check.py`:

- Burke's coefficients are re-estimated from his panel, exactly (0.0127177 / −0.00048708, optimum 13.1 °C, not the
  12.7 °C of the rounded pair), with his 1000-draw country bootstrap, so every money number carries his uncertainty;
  the 5-lag cumulative and rich/poor specifications are sensitivity rows.
- Each country is evaluated at an observation-anchored temperature (UDel 1980–2010 mean + his CMIP5 country/global
  warming ratio × the warming to the level), Burke's own construction, instead of each model's raw level; the 30 °C cap
  is applied. The models supply the arm-minus-parent difference and the interannual anomalies.
- Δg is tested against its own same-forcing spread (parent-member pairs) rather than with a KS test on monthly g(T).
- The money is framed three ways: one year of growth at the 2019 economy, and the same difference sustained for 10 and
  30 years.
- A second experiment follows the pathways year by year from 2015 and compounds the growth response, which keeps the
  global-warming channel that matching on warming level removes.

```
burke.py                 Burke's regressions from the replication panel (self-check: python burke.py), response and baselines
country_tas.py           annual population-weighted country temperature for every run (array over models)
ramip_gdp_burke.ipynb    matched GWL 2.0: Δg, significance, money, per country and per arm
ramip_gdp_pathway.ipynb  year-matched along the pathway, compounded to 2050 and to each run's end
run_gdp.sbatch           stage 1: sbatch --array=0-8 gdp/run_gdp.sbatch
                         stage 2: sbatch --dependency=afterok:<id> --export=ALL,STAGE=nb gdp/run_gdp.sbatch
```

Outputs (`diagnostics/gdp/`): `country_tas_{model}.nc`, `countries.csv`, `burke_coef_{pooled,lag5,richpoor}.csv`,
`gdp_totals.csv`, `gdp_totals_by_model_draw.csv`, `gdp_country.csv`, `gdp_significance.csv`, `gdp_variants.csv`,
`pathway_summary.csv`, `pathway_by_model.csv`, `pathway_country_2050.csv`. Figures (`figures/e8_gdp_burke/`):
`gdp_totals`, `gdp_variants`, `gdp_map_{arm}`, `gdp_top20_glob`, `pathway_world`, `pathway_map_2050`.
