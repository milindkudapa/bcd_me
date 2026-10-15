"""Move figures and superseded diagnostics into per-experiment folders (dry run unless --apply).

Figures: summer_data/figures/<group>/ ; diagnostics of the earlier/ notebooks: aux/diagnostics/earlier/.
The notebooks' out_dir / cache paths are patched separately (see the commit that added this script).
"""
import os, sys, glob, shutil, re
sys.path.insert(0, '/burg-archive/home/mck2199/summer/bcd_me/code')
os.chdir('/burg-archive/home/mck2199/summer/bcd_me/path_independence')
from funcs_support import get_params
AUX = get_params()['aux']; FIG = AUX + '../figures/'; DIAG = AUX + 'diagnostics/'
APPLY = '--apply' in sys.argv

GROUPS = {
    'e0_foundation': ['gwl_timeseries_panelA.png', 'gwl_crossing_panelB.png', 'ssp_gwl_comparison.png', 'ssp_gwl_comparison.pdf',
                      'bcdme_fig2_local.png', 'bcdme_fig2_faithful.png'],
    'e1_global_cleanup': ['ks_blockperm_maps_ramip_gwl2.0.png', 'ks_blockperm_maps_ramip_gwl2.0_pr.png', 'ks_blockperm_maps_ramip_gwl2.0_tasmax.png',
                          'ks_blockperm_maps_ramip_gwl2.0_tasmin.png', 'ramip_fig2_regional_boxplots_gwl2.0.png', 'ramip_gwl2.0_windows.png',
                          'ks_blockperm_landfrac_vars_gwl2.0.png', 'ramip_gmst_all_exps.png', 'ks_blockperm_fracmaps_1deg_gwl2.0.png',
                          'ks_blockperm_landfrac_1deg_gwl2.0.png'],
    'e2_east_asia': ['ks_blockperm_fracmaps_1deg_gwl2.0_eas.png', 'ks_blockperm_landfrac_1deg_gwl2.0_eas.png'],
    'e3_regional_arms': ['ks_regional_arms_landfrac.png', 'ks_regional_arms_inbox.png'] + [f'ks_regional_arms_fracmaps_{v}.png' for v in ('tas', 'pr', 'gdp', 'tasmax', 'tasmin', 'dtr')],
    'e4_vs_ssp585': ['ks_regional_vs_ssp585_fracmaps.png', 'ks_regional_vs_ssp585_landfrac.png', 'diag_miroc_archives.png',
                     'diag_miroc_climatology.png', 'diag_miroc_ssp370_pangeo_vs_deela.png'],
    'e5_scenarios_overshoot': ['ks_scenarios_fracmaps_pr.png', 'ks_scenarios_fracmaps_tas.png', 'ks_scenarios_landfrac.png'],
    'e7_gdp': ['ramip_impacts_burke_gwl2.0.png', 'ramip_impacts_burke_regions_gwl2.0.png'],
    'method': ['ks_block_permutation_size.png', 'ks_blockperm_pixel_diagnostics.png', 'ks_blockperm_ensemble_pairs.png'],
    'summary': ['summary_fig1_gmst.png', 'summary_fig2_landfrac_vars.png', 'summary_fig3_regional_response.png', 'summary_fig4_panel.png',
                'story_fig1_design.png', 'story_fig2_test.png', 'story_fig3_result.png', 'story_fig4_impacts.png', 'project_overview.png',
                'ks_wilks_fdr_families.png'],
    'audit': ['ks_regional_arms_audit.png'],
    'earlier': ['ks_blockperm_maps_ramip_multi.png', 'ks_pvalue_maps_ramip_multi.png', 'ks_sample_pixels_ramip_multi.png',
                'ramip_fig2_regional_boxplots.png', 'ks_empirical_null_maps_ssp370_GWL2.png', 'ks_pvalue_maps_ssp370_GWL2.png',
                'ks_pvalue_maps_ramip_giss.png'],
}
DIAG_EARLIER = ['ks_blockperm_ramip', 'ks_ramip_multi', 'ks_empirical_null', 'ramip_regional_means',
                'ks_pvals_ramip_giss_370_vs_126aer.nc', 'ks_pvals_ssp370_vs_others_GWL2_allruns.nc', 'ks_pvals_ssp370_vs_others_GWL2.nc',
                'gwl_test_data_gwl2_faithful.nc', 'gwl_test_rawchange_gwl2_faithful.nc']

moves = []
known = {f for fs in GROUPS.values() for f in fs}
for grp, files in GROUPS.items():
    for f in files:
        if os.path.exists(FIG + f):
            moves.append((FIG + f, FIG + grp + '/' + f))
loose = sorted(f for f in os.listdir(FIG) if os.path.isfile(FIG + f) and f not in known)
for item in DIAG_EARLIER:
    if os.path.exists(DIAG + item):
        moves.append((DIAG + item, DIAG + 'earlier/' + item))
print(f'{len(moves)} moves; loose figures not in any group: {loose}')
for src, dst in moves:
    print(f'  {os.path.relpath(src, AUX + "..")}  ->  {os.path.relpath(dst, AUX + "..")}')
if APPLY:
    for src, dst in moves:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.move(src, dst)
    print('applied')
