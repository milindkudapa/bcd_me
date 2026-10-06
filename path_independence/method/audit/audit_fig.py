import os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out') + '/'
FIG = '/burg-archive/home/mck2199/summer/summer_data/figures/ks_regional_arms_audit.png'
pw = pd.read_csv(T + 'power.csv'); ef = pd.read_csv(T + 'effects.csv'); ym = pd.read_csv(T + 'effects_ym.csv'); ea = pd.read_csv(T + 'effects_eas.csv')
bs = pd.read_csv(T + 'boxseries.csv'); bn = pd.read_csv(T + 'boxseries_null.csv')
if os.path.exists(T + 'boxseries_eas.csv'): bs = pd.concat([bs, pd.read_csv(T + 'boxseries_eas.csv')])
e_g = ea[ea.kind == 'gwl'].assign(arm='eas', box='eas'); e_y = ea[ea.kind == 'ym'].assign(arm='eas', box='eas')
ef = pd.concat([ef, e_g[ef.columns.intersection(e_g.columns)]]); ym = pd.concat([ym, e_y[ym.columns.intersection(e_y.columns)]])
LAB = {'eas': 'EAS126aer', 'sas': 'SAS126aer', 'asia': 'ASIA126aer', 'afr': 'AFR126aer', 'nae': 'NAE126aer', 'safca': 'SAF126ca', 'sasca': 'SAS126ca'}
fig, axs = plt.subplots(2, 3, figsize=(18, 9.5))
# (a, b) power
for k, (var, unit, scale, bands) in enumerate([('tas', 'imposed shift in monthly tas (K)', 1, [(0, 0.25, 'tab:green', 'regional arms: in-box ΔT at matched GWL (0–0.25 K)'), (0.2, 0.7, 'tab:orange', 'global cleanup: in-box ΔT (0.2–0.7 K)')]),
                                               ('pr', 'imposed multiplicative change in monthly pr (%)', 100, [(0, 4, 'tab:green', 'regional arms: in-box Δpr (0–4 %)'), (1, 7, 'tab:orange', 'global cleanup: in-box Δpr (1–7 %)')])]):
    ax = axs[k, 0]
    for mode, ls, tag in [('raw', '-', 'monthly values (as used)'), ('anom', '--', 'common seasonal cycle removed')]:
        g = pw[(pw['var'] == var) & (pw['mode'] == mode)].groupby('shift')[['all', 'tropics', 'extratropics']].mean()
        for col, c in [('all', 'k'), ('tropics', 'tab:red'), ('extratropics', 'tab:blue')]:
            ax.plot(g.index * scale, g[col], ls, color=c, label=f'{col} land, {tag}')
    for lo, hi, c, lab in bands: ax.axvspan(lo, hi, color=c, alpha=0.15, label=lab)
    ax.axhline(0.05, color='0.5', lw=0.8); ax.axhline(0.5, color='0.5', lw=0.8, ls=':')
    ax.set_xlabel(unit); ax.set_ylabel('fraction of land pixels with p<0.05'); ax.set_ylim(0, 1)
    ax.set_title(f'({"ab"[k]}) power of the per-pixel block-permutation KS, {var}\nCanESM5-1 same-forcing member pairs + imposed change, n=120', fontsize=10); ax.legend(fontsize=6.5, loc='upper left' if k == 0 else 'lower right')
# (c, d) effect size vs matched-window year
cols = {'eas': 'tab:red', 'sas': 'tab:orange', 'asia': 'tab:purple', 'afr': 'tab:green', 'nae': 'tab:blue'}
for k, (var, ylab) in enumerate([('tas', 'in-box ΔT, arm − parent (K)'), ('pr', 'in-box Δpr, arm − parent (% of parent)')]):
    ax = axs[k, 1]
    for arm, c in cols.items():
        d = ef[(ef['var'] == var) & (ef.arm == arm)].groupby('model').agg(yr=('arm_end', 'mean'), d=('d_in', 'mean'), n=('run', 'size'))
        ax.scatter(d.yr, d.d, color=c, s=12 + 4 * d.n, label=f'{LAB[arm]}: per model at its matched 2 °C decade')
        ax.axhline(ym[(ym['var'] == var) & (ym.arm == arm)].d_in.mean(), color=c, ls='--', lw=1)
    ax.axhline(0, color='k', lw=0.6); ax.set_xlabel('end year of the matched 2 °C decade'); ax.set_ylabel(ylab)
    ax.set_title(f'({"cd"[k]}) in-box response vs when the model reaches 2 °C, {var}\n(dot size = pairs; dashed = year-matched 2041–2050 multi-model mean)', fontsize=10); ax.legend(fontsize=6.5)
# (e, f) box-mean series test
order = ['eas', 'sas', 'asia', 'afr', 'nae', 'safca', 'sasca']
for k, var in enumerate(['tas', 'pr']):
    ax = axs[k, 2]
    rej = bs[bs['var'] == var].groupby(['arm', 'box']).p.apply(lambda s: (s < .05).mean())
    nul = bn[bn['var'] == var].groupby('box').p.apply(lambda s: (s < .05).mean())
    x = np.arange(len(order))
    ax.bar(x - 0.2, [rej.get((a, a), np.nan) for a in order], 0.4, color='tab:red', label='regional arm vs own parent')
    ax.bar(x + 0.2, [rej.get(('glob', a), np.nan) for a in order], 0.4, color='tab:orange', label='global cleanup (126aer) vs parent, same box')
    ax.hlines([nul[a] for a in order], x - 0.45, x + 0.45, color='k', ls='--', lw=1.2, label='same-forcing parent pairs, same box (null)')
    ax.set_xticks(x, [LAB[a] for a in order], fontsize=8); ax.set_ylabel('fraction of pairs with p<0.05')
    ax.set_title(f'({"ef"[k]}) block-permutation KS on the box-mean monthly anomaly series, {var}\nGWL-matched; one series per pair, n=120', fontsize=10); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(FIG, dpi=150, bbox_inches='tight'); print('saved', FIG)
