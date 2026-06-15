"""
Supplementary Figure: Scatter plots for CKD (binary) outcome
Analogous to Figure 2 (eGFR scatter plots)
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats as scipy_stats

proc = "data/processed/"
figs = "figures/"

# ── Data ──────────────────────────────────────────────────────────────────────
ivw     = pd.read_csv(proc + "mr_results_ivw.csv")
all_r   = pd.read_csv(proc + "mr_results_all.csv")
harm_ckd = pd.read_csv(proc + "harmonised_ckd.csv")

# IVW estimates for CKD
ckd_ivw = ivw[ivw['outcome'] == 'CKD (binary)'].set_index('exposure')

# FDR significance
ivw_fdr = ckd_ivw['pval_fdr'].to_dict()

# ── Exposure order + colors (DAL-role based, same as Figure 1/2) ──────────────
EXP_ORDER_SCATTER = [
    ('Tea intake',                 'Tea',          '#2E8B57'),
    ('Salad/raw vegetable intake', 'Salad/raw veg','#2E8B57'),
    ('Dried fruit intake',         'Dried fruit',  '#2E8B57'),
    ('Cooked vegetable intake',    'Cooked veg',   '#2E8B57'),
    ('Fresh fruit intake',         'Fresh fruit',  '#2E8B57'),
    ('Oily fish intake',           'Oily fish',    '#B8860B'),
    ('Processed meat intake',      'Proc. meat',   '#C0392B'),
    ('Pork intake',                'Pork',         '#C0392B'),
    ('Beef intake',                'Beef',         '#C0392B'),
    ('Lamb/mutton intake',         'Lamb/mutton',  '#C0392B'),
]

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 10.5), facecolor='white')
gs  = gridspec.GridSpec(3, 4, figure=fig, hspace=.44, wspace=.36,
                        top=.94, bottom=.05, left=.07, right=.97)

fig.text(.5, .975,
         'Solid line = IVW slope; dashed line = MR-Egger slope.  '
         'Tea (green): the only FDR-significant CKD association.',
         ha='center', va='top', fontsize=8.5, color='#555')

for idx, (exp_full, exp_short, col) in enumerate(EXP_ORDER_SCATTER):
    ax = fig.add_subplot(gs[idx // 4, idx % 4])

    grp  = harm_ckd[harm_ckd['exposure'] == exp_full].copy()
    be   = grp['beta_exposure'].values
    bo   = grp['beta_outcome_aligned'].values   # log-OR for CKD
    se_o = grp['se_outcome'].values
    n_s  = len(grp)

    # IVW slope (from stored results)
    ivw_b = ckd_ivw.loc[exp_full, 'beta']   # log-OR scale

    # MR-Egger (weighted)
    w      = 1 / se_o**2
    sw     = np.sum(w)
    sbe    = np.sum(w * be)
    sbo    = np.sum(w * bo)
    sbe2   = np.sum(w * be**2)
    sbebo  = np.sum(w * be * bo)
    D      = sw * sbe2 - sbe**2
    eg_b   = (sw * sbebo - sbe * sbo) / D
    eg_a   = (sbo - eg_b * sbe) / sw

    # Plot
    ax.errorbar(be, bo, yerr=1.96 * se_o,
                fmt='none', ecolor='#CCCCCC', elinewidth=.6, zorder=1)
    ax.scatter(be, bo, color=col, s=22, zorder=3,
               edgecolors='white', linewidths=.4)

    xr = np.array([be.min() * 1.15, be.max() * 1.15])
    ax.plot(xr, ivw_b * xr,    color=col,   lw=1.6, zorder=4)
    ax.plot(xr, eg_a + eg_b * xr, color='#333', lw=1.2, ls='--', zorder=4)

    ax.axhline(0, color='#AAAAAA', lw=.6)
    ax.axvline(0, color='#AAAAAA', lw=.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6.5)

    # FDR dagger
    dag = '†' if ivw_fdr.get(exp_full, 1.) < 0.05 else ''

    ax.set_title(exp_short, fontsize=8.5, fontweight='bold', color=col, pad=2)
    ax.set_xlabel(f"SNP→{exp_short[:10]}", fontsize=6.5, labelpad=2)
    ax.set_ylabel('SNP→log-OR (CKD)', fontsize=6.5, labelpad=2)

    # OR and p in annotation
    or_val = np.exp(ivw_b)
    p_fdr  = ivw_fdr.get(exp_full, np.nan)
    ax.text(.04, .97,
            f"N={n_s}  OR={or_val:.3f}{dag}",
            transform=ax.transAxes, fontsize=7, va='top', color=col)

# ── Empty 11th + 12th subplot (3×4 grid has 12 slots, we use 10) ─────────────
for slot in [(2, 2), (2, 3)]:
    fig.add_subplot(gs[slot]).axis('off')

# ── Save ──────────────────────────────────────────────────────────────────────
out = figs + 'FigureS_scatter_CKD.png'
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}")
plt.close()
