# Audit: aerosol path-(in)dependence at a matched warming level

Date: 2026-10-05. Scope: the block-permutation KS pipeline as applied to the RAMIP aerosol
experiments (`pipeline/ks_gwl_1deg.ipynb`, `regional/ks_regional_vs_ssp585.ipynb`,
`regional/ks_regional_arms.ipynb`), its summary statistics (land fractions, Wilks FDR, in-box
test), and the conclusion drawn from them. Every number below was recomputed from the caches or
from the raw zarr windows by the scripts in `method/audit/` (outputs in `method/audit/out/`,
figure `summer_data/figures/ks_regional_arms_audit.png`).

## 1. Verdict

1. **The computations are correct.** Pairing, GWL windows, land and box masks, caches and the
   job array were checked and reproduce. One real bug was found and fixed before the final run:
   `ks_tools.sig_fdr` used a strict `<` in the Benjamini–Hochberg step, which discards the pixel
   that sets the threshold and, with permutation p-values tied at the floor 1/1000, could discard
   every significant pixel. Earlier FDR numbers (`ramip_gwl/land_fractions_1deg_fdr*.csv`) are
   slightly low and should be regenerated.
2. **The headline result reproduces, but it is a detection limit, not a demonstration of
   independence.** Global land fractions sit at the 0.057 floor and FDR is ≈ 0 for every arm. The
   per-pixel test, however, needs an imposed shift of about **1 K** in monthly temperature (or
   about **25 %** in monthly precipitation) to reject at half the land pixels, whereas the regional
   cleanups change in-box means by **≤ 0.25 K and ≤ 4 %** even at full emission divergence. The
   test could not have seen these effects; "indistinguishable" here means "below ~1 K per pixel".
   Given the whole ensemble instead of one pair (7–10 members pooled), the same test rejects
   40–50 % of the cleanup box in the two models tried (§3.4b).
3. **Path dependence at a matched warming level is real, small and local.** A test on the
   box-mean monthly-anomaly series (one series per pair, same block permutation) rejects for the
   East-Asia and combined-Asia cleanups in temperature and for the Africa, Asia and South-Asia
   cleanups in precipitation, at 4–8× the same-forcing null. The effect sizes are 0.1–0.2 K and
   2–4 % of the box mean. Global cleanup leaves a +0.2 K residual over East Asia at 2 °C.
4. **GWL matching carries a timing confound the notebooks do not discuss.** A hot model reaches
   2 °C before the aerosol pathways have separated (CanESM5-1's matched decade ends in 2027,
   three years into the experiment), a cool model after they have fully separated (MIROC6: 2065).
   The in-box response grows with the matched-window year (figure panels c, d); the multi-model
   mean therefore mixes "not yet diverged" with "fully diverged" comparisons.
5. Smaller points that change how numbers should be read: the FDR test cannot flag fewer than
   108 land pixels at 999 draws; `gdp` is statistically the same test as `tas`; the test is
   liberal over the ocean (false-positive rate 0.06–0.10 vs 0.035–0.06 on land); the floor is
   model-dependent (0.035–0.062, the single 0.057 is too high for GISS and NorESM); the
   native-grid in-box null was checked against a 1° one and found equivalent.

## 2. What was checked

| Check | Method | Result |
|---|---|---|
| FDR implementation | synthetic tied p-values; all 119 cached maps old vs corrected | strict `<` dropped floor-tied pixels; fixed, regression test in `ks_tools.demo()` |
| GWL matching precision | realised 10-yr GMST anomaly in each paired window, from the data product | mean \|Δ\| 0.020 K (max 0.069) vs ssp370; 0.025 K (max 0.127, CanESM5) vs ssp585. Level 2.03 K |
| Pairing | audit of `parent_pairs`/`variant_pairs` against the product and `gwl2.0_windows.csv` | no usable run lost; EC-Earth r1–r4 lack a parent window (genuine); NorESM2-LM's `ssp370` never reaches 2 °C by 2059 |
| Same-run pairing bias | CanESM5-1, EAS and SAS at matched GWL: 10 same-member vs 90 cross-member pairs | identical false-positive behaviour (land 0.036 vs 0.039; in-box 0.029 vs 0.030): no initial-condition memory left at 2 °C |
| Emission boxes | against Wilcox et al. (2023), GMD 16, Tables 1–3 | all seven match |
| Caches, jobs | array logs, `sacct`, NaN scan, per-pair FDR | all 927 pair-tests completed; no NaNs; CSV race in the array fixed (only the FINAL run writes tables) |
| Test power | CanESM5-1 same-forcing pairs + imposed shift, raw monthly values and with a common seasonal cycle removed | §3.1 |
| Effect sizes | area-weighted in-box means, arm − parent, GWL-matched and year-matched 2041–2050 | §3.2 |
| Aggregate test | block-permutation KS on box-mean monthly-anomaly series, with a 300-pair same-forcing null per box | §3.4 |
| FDR resolution | BH step-up with p ≥ 0.001 | §3.5 |
| False-positive floor | per-model same-forcing pairs, land with/without Antarctica, ocean | §3.6: 0.035–0.062 by model; ocean 0.06–0.10 |
| `gdp` vs `tas` | cached p-maps | rejection decisions agree on 97.5 % of land pixels; r = 0.87 |
| 1° same-forcing null | CanESM5-1, 45 pairs, regridded exactly like the arms | §3.5: reproduces the native in-box null (EAS FDR 0.019 on both grids); the notebook now uses it |

## 3. Findings

### 3.1 The per-pixel test cannot see effects of the size the regional cleanups produce

Design: pairs of CanESM5-1 `ssp370` members (same forcing, 2041–2050, 1959 land pixels south of
60°S excluded), with a constant shift added to one member. Fraction of land pixels with p < 0.05,
mean of three pairs (`method/audit/power.py`, figure panels a, b).

| imposed shift in monthly `tas` | 0 | 0.25 K | 0.5 K | 0.75 K | 1.0 K | 1.5 K | 2.0 K |
|---|---|---|---|---|---|---|---|
| monthly values, as used | 0.043 | 0.071 | 0.198 | 0.361 | 0.523 | 0.757 | 0.883 |
| … tropics / extratropics | 0.03/0.05 | 0.09/0.07 | 0.37/0.14 | 0.62/0.27 | 0.80/0.43 | 0.95/0.69 | 0.99/0.85 |
| common seasonal cycle removed | 0.049 | 0.078 | 0.257 | 0.510 | 0.720 | 0.924 | 0.976 |

| imposed change in monthly `pr` | 0 | +5 % | +10 % | +20 % | +30 % | +50 % |
|---|---|---|---|---|---|---|
| monthly values, as used | 0.049 | 0.075 | 0.141 | 0.361 | 0.579 | 0.794 |
| common seasonal cycle removed | 0.043 | 0.078 | 0.167 | 0.420 | 0.636 | 0.830 |

The limitation is not the seasonal cycle: removing it buys ~30 % more power at 1 K, no more.
It is n = 120 autocorrelated months against a monthly-anomaly standard deviation of 1.1–2.0 K in
the cleanup regions (48–78 % of the mean for precipitation). A 0.2 K shift is 0.1–0.2 σ.

### 3.2 How large the regional cleanup effects are

In-box land means, arm minus own `ssp370` parent, multi-model mean of per-pair values
(`method/audit/effects.py`, `yearmatched.py`, `eas_effects.py`). "Outside" is all other land.

| arm (box) | models / pairs | GWL-matched ΔT (K) | outside | year-matched 2041–50 ΔT | GWL-matched Δpr (%) | year-matched Δpr |
|---|---|---|---|---|---|---|
| EAS126aer | 8 / 67 | **+0.14** ± 0.22 | 0.00 | +0.25 | +2.3 | +3.2 |
| SAS126aer | 7 / 62 | **+0.10** ± 0.17 | −0.01 | +0.20 | +3.9 | +3.4 |
| ASIA126aer | 2 / 16 | +0.17 ± 0.14 | +0.01 | +0.21 | +3.2 | +2.2 |
| AFR126aer | 8 / 58 | +0.03 ± 0.12 | −0.01 | +0.07 | +1.6 | +1.5 |
| NAE126aer | 8 / 66 | +0.07 ± 0.13 | 0.00 | +0.21 | +0.3 | +0.7 |
| SAF126ca | 2 / 10 | 0.00 ± 0.12 | +0.03 | +0.08 | +2.9 | +0.3 |
| SAS126ca | 3 / 12 | 0.00 ± 0.15 | +0.02 | +0.08 | +0.6 | −0.2 |
| global 126aer, in the EAS box | 8 / 72 | +0.22 | land mean +0.05 | +0.67 (land +0.51) | +3.6 | +6.9 |

Read against §3.1: every regional-arm effect lies in the first column of the power table, where
the per-pixel test rejects at 5–7 % of pixels, i.e. at the null rate. The in-box pixel-FDR signal
found for EAS (0.044) and SAS (0.022) comes from the few pixels where the local response is
several times the box mean; it is consistent with, not in tension with, these numbers.

### 3.3 The timing confound of GWL matching

RAMIP Tier-1 runs cover 2015–2051. The matched 2 °C decade ends in 2027 for CanESM5-1, 2038 for
UKESM1-0-LL, 2040–2042 for GISS-E2-1-G and EC-Earth3-AerChem, 2046 for CESM2, 2051 for
MRI-ESM2-0, 2055 for CNRM-ESM2-1 and 2065–2068 for MIROC6 (the last three use the extended runs).
SSP3-7.0 and SSP1-2.6 aerosol emissions separate gradually after 2015, so at 2 °C a hot model
compares two nearly identical forcings and a cool model two fully separated ones. The per-model
in-box response tracks the window year (figure panels c, d): for SAS126aer, +0.02 K (CanESM5-1,
2027), +0.04 (UKESM, 2038), +0.05 (GISS, 2041), +0.13 (EC-Earth, 2041), +0.06 (CESM2, 2045),
+0.22 (MRI, 2049), +0.22 K (MIROC6, 2065); for EAS126aer −0.03, −0.01, +0.10, +0.19, +0.16, +0.28,
+0.29 K in the same order. The box-mean test rejections (§3.4) come from the same late-window
models. A multi-model mean at a fixed GWL is therefore partly a statement about *when* each model
reaches that level. This is inherent to matching a transient regional forcing on global
temperature; it should be reported, and the comparison repeated at levels where all models have
diverged (§5, B).

### 3.4 An aggregate test does detect path dependence

Block-permutation KS (12-month blocks, 999 draws) on the area-weighted box-mean monthly series,
anomalies relative to the leave-one-out ensemble climatology of the other parent members, one test
per pair; null from all same-forcing parent pairs of every model (300 per box). Fraction of pairs
with p < 0.05 (`method/audit/boxseries.py`, figure panels e, f):

| box | `tas`: arm | `tas`: global cleanup | `tas`: null | `pr`: arm | `pr`: global cleanup | `pr`: null |
|---|---|---|---|---|---|---|
| EAS | **0.21** | **0.25** | 0.05 | 0.09 | **0.29** | 0.02 |
| SAS | 0.07 | 0.08 | 0.03 | 0.15 | 0.17 | 0.09 |
| ASIA | **0.19** | **0.24** | 0.01 | **0.31** | **0.31** | 0.04 |
| AFR | 0.05 | 0.00 | 0.01 | **0.22** | 0.15 | 0.04 |
| NAE | 0.00 | 0.11 | 0.05 | 0.03 | 0.08 | 0.02 |
| SAF (ca) | 0.00 | 0.00 | 0.00 | 0.10 | 0.17 | 0.06 |
| SAS (ca) | 0.00 | 0.07 | 0.03 | 0.17 | 0.17 | 0.09 |

The null rate varies by box (0.00–0.09), because box-mean series carry more interannual memory
than pixel series and 12-month blocks do not remove it (the size study in
`method/ks_block_permutation.ipynb` already shows this for persistent processes); each arm must be
read against its own box's null, as above. With that caveat: East-Asia temperature and
Africa/Asia precipitation differ between pathways at the same global temperature. The effect is
the 0.1–0.2 K / 2–4 % of §3.2.

### 3.4b Pooling the ensemble restores power — and finds the difference

Suggestion A2 of §5, tested on the two models whose matched decade falls after the emission
pathways have separated (`method/audit/pooled.py`): all usable members' windows concatenated on
each side (7–10 members, n = 840–1 200 months per side, 12-month blocks as before), against the
same test on the parents split into two halves (the pooled null) and on a single pair.

| model, arm | design | members / side | land p<0.05 | in-box p<0.05 | outside box | in-box FDR |
|---|---|---|---|---|---|---|
| MRI-ESM2-0, SAS126aer | pooled arm vs parent | 7 | 0.124 | **0.438** | 0.116 | **0.547** |
| | pooled null (parents split) | 3 | 0.044 | 0.063 | 0.044 | 0.000 |
| | single pair | 1 | 0.031 | 0.051 | 0.031 | 0.000 |
| MRI-ESM2-0, EAS126aer | pooled arm vs parent | 10 | 0.097 | **0.409** | 0.082 | **0.471** |
| | pooled null | 5 | 0.021 | 0.003 | 0.022 | 0.000 |
| | single pair | 1 | 0.040 | 0.041 | 0.040 | 0.000 |
| MIROC6, SAS126aer | pooled arm vs parent | 10 | 0.121 | **0.488** | 0.112 | **0.685** |
| | pooled null | 5 | 0.027 | 0.024 | 0.027 | 0.000 |
| | single pair | 1 | 0.056 | 0.071 | 0.056 | 0.000 |
| MIROC6, EAS126aer | pooled arm vs parent | 10 | 0.101 | **0.525** | 0.080 | **0.621** |
| | pooled null | 5 | 0.027 | 0.110 | 0.023 | 0.000 |
| | single pair | 1 | 0.122 | 0.246 | 0.116 | 0.246 |

With the ensemble pooled, 40–50 % of the cleanup box rejects and 47–69 % of it passes FDR,
against ≤ 11 % and 0 for the pooled null; outside the box 8–12 % of land rejects against 2–4 %.
The single-pair design — the one every result table in the project uses — sees 4–7 % in the
same box. At 2 °C the regional cleanups are therefore clearly distinguishable from their parent
once the test is given the whole ensemble: the pixel-level "indistinguishable" verdict was a
sample-size statement. The 8–12 % outside the box is consistent with pattern effects of a few
tenths of a kelvin in some regions (the outside-box *mean* is ≈ 0, §3.2) and should be mapped
before it is interpreted. Only two models were run; the pooled design costs ~10× a single pair.

### 3.5 What the FDR numbers can and cannot resolve

With 999 permutations the smallest p-value is 1/1000. The BH step-up flags rank *i* only if
p(i) ≤ 0.2·i/N, so with p ≥ 0.001 it needs i ≥ 0.005 N: on the 1° grid (N = 21 537 land pixels,
Antarctica included) **at least 108 pixels must sit at the permutation floor before any pixel is
significant**, and the smallest nonzero FDR fraction is 0.0050 — exactly the minimum observed
over the 309 cached pairs (22 nonzero). "FDR ≈ 0" therefore means "fewer than 108 land pixels at
p = 0.001 per map", a blunt instrument; more draws (9 999) would lower the threshold to ~11 pixels.
Inside a box the threshold is 5 pixels (EAS, 987 land pixels at 1°) but **1 pixel on CanESM5-1's
native grid (65 pixels in the SAS box)**, so the native-grid in-box null the notebook first used was not obviously on the same footing
as the 1° in-box FDR values it was compared with. The 1° null computed for this audit
(`method/audit/null1deg.py`, 45 CanESM5-1 pairs regridded exactly like the arms; in-box values
from `inbox_null.py`) settles it:

| box | land px at 1° | `tas` raw / FDR, native 2.8° | `tas` raw / FDR, 1° | `pr` raw / FDR, 1° | `gdp` raw / FDR, 1° | arm vs own parent, in-box FDR (`tas`) |
|---|---|---|---|---|---|---|
| EAS | 987 | 0.057 / 0.019 | 0.059 / 0.019 | 0.062 / 0.008 | 0.055 / 0.010 | **0.044** |
| SAS | 501 | 0.030 / 0.003 | 0.029 / 0.000 | 0.031 / 0.000 | 0.038 / 0.000 | **0.023** |
| AFR | 2 938 | 0.028 / 0.000 | 0.028 / 0.000 | 0.036 / 0.000 | 0.031 / 0.000 | 0.009 |
| NAE | 3 703 | 0.051 / 0.006 | 0.050 / 0.005 | 0.040 / 0.000 | 0.047 / 0.001 | 0.003 |
| ASIA | 1 488 | 0.048 / 0.005 | 0.049 / 0.003 | 0.051 / 0.003 | 0.049 / 0.001 | 0.003 |
| SAF | 1 456 | 0.023 / 0.000 | 0.025 / 0.000 | 0.035 / 0.000 | 0.025 / 0.000 | 0.000 |

The 1° null reproduces the native one to within 0.003: regridding a 2.8° field to 1° turns one
native cell into ~8 correlated cells, so the pixel-count threshold scales with the grid and the
original comparison was fair after all. The notebook now uses the 1° per-variable null anyway.
What the table does show is that the false-positive rate differs between boxes by a factor of
2.5 (SAF 0.025, EAS 0.059 — East Asia has the most decadal variability of the seven, which is
also why its null FDR is the highest), so in-box rates must be read against their own box; and
that the EAS and SAS arms exceed their null FDR (2.3× and 0.023 vs 0.000), AFR marginally, NAE and
ASIA not at all.

Antarctica is in N for every land fraction and FDR (the 0.057 floor was measured without it).
It dilutes N by ~20 % and adds pixels where the test is conservative; harmless for the conclusions
but worth making consistent.

### 3.6 The false-positive floor is model- and region-dependent

The 0.057 floor comes from one model (CanESM5-1) in one decade (2041–2050). At its matched 2 °C
windows (ending ~2027) the same model's floor is 0.040 over all land. Per box (CanESM5-1, matched
windows) the same-forcing rejection rate ranges from 0.022 (sub-Saharan Africa) to 0.068 (North
America/Europe), so in-box comparisons need a box-specific null — the notebook does this. Over the
ocean (60°S–60°N) the same pairs reject 0.082: the dark North Atlantic and Southern Ocean patches
in the maps are mostly the test's own liberality under multi-year SST memory, not a signal.
Per-model floors (`method/audit/floors.py`, 10 same-forcing pairs per model, `tas`, 2041–2050):

| model | land, with Antarctica | land, without | ocean 60°S–60°N | native land px |
|---|---|---|---|---|
| CESM2 | 0.053 | 0.055 | 0.072 | 18 402 |
| CNRM-ESM2-1 | 0.052 | 0.062 | 0.067 | 10 859 |
| CanESM5-1 | 0.050 | 0.062 | 0.084 | 2 689 |
| EC-Earth3-AerChem | 0.047 | 0.049 | 0.062 | 43 413 |
| GISS-E2-1-G | 0.032 | 0.035 | 0.069 | 4 297 |
| MIROC6 | 0.046 | 0.052 | 0.059 | 10 859 |
| MRI-ESM2-0 | 0.052 | 0.061 | 0.080 | 16 955 |
| NorESM2-LM | 0.036 | 0.039 | 0.086 | 4 631 |
| UKESM1-0-LL | 0.047 | 0.050 | 0.101 | 9 164 |

The floor is 0.035–0.062 depending on the model (spread across the 10 pairs 0.01–0.05). The single
0.057 is right for CanESM5-1, CNRM-ESM2-1 and MRI-ESM2-0 and ~0.02 too high for GISS-E2-1-G and
NorESM2-LM: with a per-model floor, GISS's raw land fractions in the arms table (0.027–0.045)
read as "at the floor", not "below it". Over the ocean every model is liberal (0.06–0.10).

### 3.7 `gdp` is not an independent test

A KS test is invariant under monotone transformations. Burke's g(T) is monotone except across
12.7 °C, so `gdp` and `tas` give the same decision at 97.5 % of land pixels (p-map correlation
0.87; the remaining differences are Monte-Carlo noise plus pixels whose monthly temperatures
straddle the optimum). The README's "genuinely different test" holds only for those pixels. The
`gdp` columns are fine to keep but should not be counted as a second variable.

### 3.8 The ssp585 anchor

Against `ssp585` the in-box signals weaken (EAS 0.21 → 0.001, SAS 0.022 → 0.008 FDR). Two
reasons, both consistent with the data: fewer pairs (≤ 25), and SSP5 assumes strong air-quality
controls, so at 2 °C `ssp585`'s regional aerosol burden is likely closer to the cleaned arm than
`ssp370`'s is. The anchor choice sets the aerosol contrast; "vs ssp585" is a different question
(different GHG pathway *and* different aerosol pathway) rather than a stronger version of the
same one. MIROC6 must stay excluded there (DATA.md §5).

### 3.9 Things that are fine

Block permutation size on real pixels (0.051 in the method study, 0.040–0.046 here); the GWL
window definition and its precision; transform-before-regrid; conservative regridding for `pr`;
the pairing within `p<n>f<n>` variants; the exclusion of NorESM2-LM and of EC-Earth r1–r4; the
Phipson–Smyth p-value; the job array and caches after the fixes of 2026-10-05.

## 4. Does this validate path independence for aerosols?

Not in the strong sense. What the evidence supports:

- **Bounded, not absent.** At 2 °C of global warming, a regional aerosol cleanup changes the
  cleanup region's decadal-mean temperature by 0.1–0.2 K and precipitation by 2–4 % relative to
  the uncleaned pathway, and other land by ~0. Global cleanup leaves +0.2 K over East Asia and
  +0.05 K over land as a whole. These are the pattern effects that survive GWL matching.
- **Detectable with the right test.** The box-mean series test sees them (EAS, ASIA temperature;
  AFR, ASIA, SAS precipitation), and the per-pixel test sees them as soon as the ensemble is pooled
  (40–50 % of the box, in-box FDR 0.5–0.7, in MRI-ESM2-0 and MIROC6). One pair of 10-year windows
  cannot, by a factor of ~5 in shift.
- **Small compared with what matters elsewhere in the project.** The overshoot comparison
  rejects 9–18 % of land and survives FDR at 0.077; the aerosol effects are an order of magnitude
  smaller in the test's own units.
- **Conditional on timing.** The numbers are for the decade each model reaches 2 °C, which for
  the hot half of the ensemble precedes most of the emission divergence.

The defensible statement is: *at a fixed global warming level, regional aerosol pathway
differences of the RAMIP size leave local monthly climate indistinguishable at the pixel level and
shift regional means by at most ~0.2 K and ~4 %, with the largest and most consistent effects in
the cleanup region itself.* The current README sentence that the East-Asia fingerprint "fails
Wilks FDR" should be replaced by this.

## 5. Suggested next experiments (RAMIP data on disk unless noted)

**A. Make the measurement match the signal (do first)**

1. Report effect sizes with uncertainty — in-box and AR6-region ΔT, Δpr per arm with the
   ensemble spread — alongside p-values. Frame independence as an equivalence test (TOST): "the
   regional response lies within ±0.2 K" is a positive, falsifiable claim; "we failed to reject"
   is not.
2. Pool members: test the ensemble-pooled window (k members × 120 months per side, blocks =
   member-years). Power scales with √n; ten members cut the detectable shift from ~1 K to
   ~0.3 K per pixel. Demonstrated in §3.4b for two models: 40–50 % of the box rejects against
   4–7 % for a single pair. Run it for every model and arm (and for the overshoot family, where
   it would turn the 9–18 % into a map).
3. Use the box-mean / region-mean series test of §3.4 for all 46 AR6 regions, with the
   box-specific null; this is the cheapest way to map where pathways differ.
4. Raise `NDRAWS` to 9 999 for anything that feeds the FDR; 999 draws makes FDR a "≥ 108 pixels"
   detector. Cost ×10 (the largest array task ran 73 min at 999).
5. Replace the single CanESM5-1 floor with per-model, per-box nulls: `floors.py` gives the
   per-model values (0.035–0.062) from 10 pairs each; the 1° in-box null exists and the notebook
   uses it.
6. Remove the common seasonal cycle before testing (small gain, no cost) and, for the ocean,
   use 24-month blocks or exclude it explicitly.

**B. Remove the timing confound**

7. Report the emission (or GMST-lead) divergence at each matched window, and repeat the test at
   2.5 and 3.0 °C where the extended runs allow (CESM2 to 2079, GISS `126aer`/`AFR` to 2070,
   MIROC6 and UKESM1-0-LL to 2100, NorESM2-LM to 2059–2064): there every model has fully
   diverged. CanESM5-1 reaches 3 °C around 2045, inside Tier 1.
8. Alternatively fix the decade (2041–2050) and match the parent's window to the arm's GMST
   rather than to a round level — same design, model-specific level, maximal divergence.

**C. New comparisons inside RAMIP**

9. Arm vs arm at matched GWL (EAS126aer vs SAS126aer, NAE vs AFR, …), same run id: two
   *different* paths to the same level, the cleanest version of the path question. All pairs are
   on disk.
10. Additivity: ASIA126aer vs EAS126aer + SAS126aer (CNRM, GISS, MIROC6) — the non-linearity
    RAMIP was designed to expose (Wilcox et al. 2023); and global 126aer vs the sum of the four
    regional arms (9 models).
11. Species: SAS126aer vs SAS126ca and AFR126aer vs SAF126ca (3 models; boxes differ for AFR) —
    whether the sulfate component carries the response.
12. Remote responses: apply the in-box machinery to literature teleconnection regions (NAE
    cleanup → Sahel rainfall and North Atlantic; EAS → North Pacific and western North America;
    SAS → monsoon core), with their own nulls.

**D. Variables and extremes**

13. `tasmax`/`tasmin` are on disk for every arm and were only tested for global cleanup; the
    direct shortwave effect of aerosols shows first in `tasmax` and the diurnal range.
14. Download `rsds` (and `clt`, `psl`): surface shortwave is the forcing fingerprint and will show
    path dependence where temperature cannot; `psl` tests the circulation response. The fixed-SST
    `piClim-370-*` companions would separate fast (atmospheric) from slow (SST-mediated)
    responses; the fast part is what differs at fixed GWL.
15. Tail statistics (counts of months above a parent-defined 95th percentile, block-bootstrapped):
    KS is weakest in the tails, which is where impacts live.

**E. Framing**

16. Treat the 7 arms × 2 anchors × 3 variables as a family: pre-specify the primary comparison
    (own parent, `tas`, in-box) and report the rest as secondary, or fit one hierarchical model of
    the in-box effect across models and arms.

## 6. Reproducibility

`method/audit/` (run from any directory; each script imports `common.py`, which chdirs to
`path_independence/` so `funcs_support.get_params` finds `dir_list.csv`; ESMF variables as in
`run/run_regional_arms_array.sbatch`):

| script | what it computes | output |
|---|---|---|
| `quick.py` | data spans, GMST mismatch, `gdp` vs `tas`, FDR resolution | stdout |
| `power.py` | §3.1 power curves | `out/power.csv` |
| `effects.py`, `yearmatched.py`, `eas_effects.py` | §3.2 in-box effects, GWL- and year-matched | `out/effects*.csv` |
| `crossmember.py` | same-member vs cross-member pairing, per-box nulls (CanESM5-1) | `out/crossmember.csv` |
| `boxseries.py`, `boxseries_eas.py` | §3.4 box-mean series test and null | `out/boxseries*.csv` |
| `null1deg.py` | §3.5 1° same-forcing null | `diagnostics/regional_arms/null_1deg_{tas,pr,gdp}.nc` |
| `floors.py` | §3.6 per-model floors | `out/floors.csv` |
| `pooled.py` | pooled-member test (suggestion A2) | `out/pooled.csv` |
| `audit_fig.py` | figure | `figures/ks_regional_arms_audit.png` |

All of it ran on an allocated compute node (slurm job 10225324); the three heavy scripts take
~1 h together at 32 threads.
