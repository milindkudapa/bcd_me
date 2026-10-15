"""Point every notebook at the per-experiment output folders (dry run unless --apply); pairs with reorganize_outputs.py."""
import sys, re, nbformat, glob
APPLY = '--apply' in sys.argv
FIG_GROUP = {
    'pipeline/ramip_data_layer.ipynb': 'e1_global_cleanup', 'pipeline/ks_path_independence_ramip_gwl.ipynb': 'e1_global_cleanup',
    'pipeline/ks_path_independence_ramip_gwl_vars.ipynb': 'e1_global_cleanup', 'pipeline/ks_gwl_1deg.ipynb': 'e1_global_cleanup',
    'pipeline/ks_location_scenarios.ipynb': 'e5_scenarios_overshoot', 'pipeline/ramip_impacts_gdp.ipynb': 'e7_gdp',
    'pipeline/ramip_summary.ipynb': 'summary', 'pipeline/project_figure.ipynb': 'summary', 'pipeline/story_figures.ipynb': 'summary',
    'method/ks_block_permutation.ipynb': 'method', 'method/ks_blockperm_diagnostics.ipynb': 'method',
    'regional/ks_regional_vs_ssp585.ipynb': 'e4_vs_ssp585', 'regional/ks_regional_arms.ipynb': 'e3_regional_arms',
    'earlier/figure_ssps_gwls_local.ipynb': 'earlier', 'earlier/ks_empirical_null.ipynb': 'earlier',
    'earlier/ks_path_independence_maps.ipynb': 'earlier', 'earlier/ks_path_independence_ramip_multi.ipynb': 'earlier',
}
DIAG_EARLIER = ['ks_blockperm_ramip/', 'ks_ramip_multi/', 'ks_empirical_null/', 'ramip_regional_means/',
                'ks_pvals_ramip_giss_370_vs_126aer.nc', 'ks_pvals_ssp370_vs_others_GWL2_allruns.nc', 'ks_pvals_ssp370_vs_others_GWL2.nc']
OUT_RE = re.compile(r"out_dir\s*=\s*dir_list\['aux'\]\s*\+\s*'\.\./figures/'")
for fn in sorted(glob.glob('pipeline/*.ipynb') + glob.glob('method/*.ipynb') + glob.glob('regional/*.ipynb') + glob.glob('earlier/*.ipynb')):
    nb = nbformat.read(fn, 4); changes = []
    for c in nb.cells:
        if c.cell_type != 'code':
            continue
        s = c.source
        if fn in FIG_GROUP and OUT_RE.search(s):
            s = OUT_RE.sub(f"out_dir = dir_list['aux'] + '../figures/{FIG_GROUP[fn]}/'\nimport os; os.makedirs(out_dir, exist_ok=True)", s)
            changes.append(f'out_dir -> figures/{FIG_GROUP[fn]}/')
        if fn == 'pipeline/ks_gwl_1deg.ipynb' and 'gwl{GWL}{suf}.png' in s:
            s = s.replace("out_dir + f'ks_blockperm_fracmaps_1deg_gwl{GWL}{suf}.png'", "(out_dir_eas if suf else out_dir) + f'ks_blockperm_fracmaps_1deg_gwl{GWL}{suf}.png'")
            s = s.replace("out_dir + f'ks_blockperm_landfrac_1deg_gwl{GWL}{suf}.png'", "(out_dir_eas if suf else out_dir) + f'ks_blockperm_landfrac_1deg_gwl{GWL}{suf}.png'")
            if 'out_dir_eas =' not in s:
                s = "out_dir_eas = dir_list['aux'] + '../figures/e2_east_asia/'\nimport os; os.makedirs(out_dir_eas, exist_ok=True)\n" + s
            changes.append('EAS figures -> figures/e2_east_asia/')
        for item in DIAG_EARLIER:
            if f"'diagnostics/{item}" in s or f'diagnostics/{item}' in s:
                s = s.replace(f'diagnostics/{item}', f'diagnostics/earlier/{item}')
                changes.append(f'diagnostics/{item} -> diagnostics/earlier/')
        c.source = s
    if changes:
        print(f'{fn}: ' + '; '.join(sorted(set(changes))))
        if APPLY:
            nbformat.write(nb, fn)
# the audit figure script
s = open('method/audit/audit_fig.py').read()
if "figures/ks_regional_arms_audit.png" in s:
    print('method/audit/audit_fig.py: figure -> figures/audit/')
    if APPLY:
        open('method/audit/audit_fig.py', 'w').write(s.replace("figures/ks_regional_arms_audit.png", "figures/audit/ks_regional_arms_audit.png"))
print('applied' if APPLY else 'dry run')
