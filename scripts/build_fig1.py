"""
Figure 1 – Extended Forest Plot
Robust rewrite: all text placed via figure-level transforms (no mixed
axes-fraction / data-coord annotations that break on resize).
Target: DIN A4 portrait, all text readable, left whitespace minimal,
no top legend/title, bottom legend only.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FixedFormatter
import matplotlib.transforms as mtransforms

# ── Data ──────────────────────────────────────────────────────────────────────
proc = "data/processed/"
ivw  = pd.read_csv(proc + "mr_results_ivw.csv")
harm = pd.read_csv(proc + "harmonised_egfr.csv")

EXP_ORDER = [
    'Tea intake',
    'Salad/raw vegetable intake',
    'Dried fruit intake',
    'Cooked vegetable intake',
    'Fresh fruit intake',
    'Oily fish intake',
    'Processed meat intake',
    'Pork intake',
    'Beef intake',
    'Lamb/mutton intake',
]
EXP_SHORT = {
    'Tea intake':                  'Tea',
    'Salad/raw vegetable intake':  'Salad/raw veg',
    'Dried fruit intake':          'Dried fruit',
    'Cooked vegetable intake':     'Cooked veg',
    'Fresh fruit intake':          'Fresh fruit',
    'Oily fish intake':            'Oily fish',
    'Processed meat intake':       'Processed meat',
    'Pork intake':                 'Pork',
    'Beef intake':                 'Beef',
    'Lamb/mutton intake':          'Lamb/mutton',
}
COLORS = {
    'Tea intake':                  '#2E8B57',  # base-forming green
    'Salad/raw vegetable intake':  '#2E8B57',
    'Dried fruit intake':          '#2E8B57',
    'Cooked vegetable intake':     '#2E8B57',
    'Fresh fruit intake':          '#2E8B57',
    'Oily fish intake':            '#B8860B',  # mild acid gold
    'Processed meat intake':       '#C0392B',  # acid-forming red
    'Pork intake':                 '#C0392B',
    'Beef intake':                 '#C0392B',
    'Lamb/mutton intake':          '#C0392B',
}
DAL_ROLE = {
    'Tea intake':                  'Weakly basic',
    'Salad/raw vegetable intake':  'Base-forming',
    'Dried fruit intake':          'Base-forming',
    'Cooked vegetable intake':     'Base-forming',
    'Fresh fruit intake':          'Base-forming',
    'Oily fish intake':            'Mild acid-forming',
    'Processed meat intake':       'Acid-forming',
    'Pork intake':                 'Acid-forming',
    'Beef intake':                 'Acid-forming',
    'Lamb/mutton intake':          'Acid-forming',
}
DAL_COLOR = {
    'Weakly basic':     '#2E8B57',
    'Base-forming':     '#2E8B57',
    'Mild acid-forming':'#B8860B',
    'Acid-forming':     '#C0392B',
}

egfr_df = ivw[ivw['outcome']=='eGFR (log-transformed)'].set_index('exposure')
ckd_df  = ivw[ivw['outcome']=='CKD (binary)'].set_index('exposure')
inst    = harm.groupby('exposure').agg(N=('rsid','count'), MeanF=('F_stat','mean'))
inst['MeanF'] = inst['MeanF'].round(1)

N = len(EXP_ORDER)

# ── P-value formatter ─────────────────────────────────────────────────────────
def fmt_p(p):
    if p < 0.0001:
        e = int(np.floor(np.log10(p)))
        m = p / 10**e
        sup_map = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')
        exp_str = '⁻' + str(abs(e)).translate(sup_map)
        return f"{m:.1f}×10{exp_str}"
    elif p < 0.001:
        return f"{p:.4f}"
    else:
        return f"{p:.3f}"

# ── Figure setup ──────────────────────────────────────────────────────────────
# A4 portrait at 150 dpi → 1240 × 1754 px
# We use 16.54 × 11.69 inch (landscape A4) — wide enough for all columns,
# tall enough for two panels
FIG_W = 16.54
FIG_H = 11.69
FONT  = 8.5      # base font size (pt)

fig = plt.figure(figsize=(FIG_W, FIG_H), facecolor='white')

# ── Panel geometry (in figure fractions) ─────────────────────────────────────
# Left text block:  0.00 – 0.38  (exposure name, DAL, N, MeanF)
# Forest plot:      0.38 – 0.64
# Right text block: 0.64 – 1.00  (est, CI, P, P_FDR)
#
# Vertical:  Panel A top=0.96 bot=0.54  Panel B top=0.50 bot=0.08
# Legend at bottom: 0.01–0.06

FOREST_L = 0.38   # left edge of forest axes (figure fraction)
FOREST_W = 0.27   # width of forest axes
PANEL_H  = 0.40   # height of each panel
AX_E_BOT = 0.545  # bottom of Panel A forest axes
AX_C_BOT = 0.07   # bottom of Panel B forest axes

ax_e = fig.add_axes([FOREST_L, AX_E_BOT, FOREST_W, PANEL_H])
ax_c = fig.add_axes([FOREST_L, AX_C_BOT, FOREST_W, PANEL_H])

# ── Row y-positions (data coords 0..N-1, top row = N-1) ──────────────────────
Y = list(range(N-1, -1, -1))   # Y[0]=Tea at top

# ── Helper: figure-x from axes-fraction ──────────────────────────────────────
def ax_to_fig_x(ax, frac):
    """Convert axes-fraction x to figure-fraction x."""
    bbox = ax.get_position()
    return bbox.x0 + frac * bbox.width

def data_to_fig_y(ax, y_data):
    """Convert data-y to figure-fraction y (approximate, uses bbox)."""
    ylim = ax.get_ylim()
    bbox = ax.get_position()
    frac = (y_data - ylim[0]) / (ylim[1] - ylim[0])
    return bbox.y0 + frac * bbox.height

# ── Draw forest panel ─────────────────────────────────────────────────────────
def draw_forest(ax, df, outcome, label, title):
    ax.set_facecolor('white')
    ax.set_ylim(-0.7, N - 0.3)

    if outcome == 'egfr':
        XMIN, XMAX = -0.045, 0.052
        NULL = 0.0
        xlabel = 'β per SD increase in dietary exposure  [95% CI]'
    else:
        XMIN, XMAX = 0.30, 2.20
        NULL = 1.0
        xlabel = 'Odds Ratio for CKD  [95% CI]'

    ax.set_xlim(XMIN, XMAX)

    # Alternating row shading
    for i, y in enumerate(Y):
        c = '#F4F6F8' if i % 2 == 0 else 'white'
        ax.axhspan(y - 0.45, y + 0.45, color=c, zorder=0)

    # Null line
    ax.axvline(NULL, color='#555', lw=1.0, ls='--', zorder=1, alpha=0.8)

    # Light reference lines
    if outcome == 'egfr':
        for ref in [-0.025, 0.025]:
            ax.axvline(ref, color='#DDD', lw=0.6, zorder=0)
    else:
        for ref in [0.5, 2.0]:
            ax.axvline(ref, color='#DDD', lw=0.6, zorder=0)

    # CI lines + markers
    for exp, y in zip(EXP_ORDER, Y):
        row = df.loc[exp]
        col = COLORS[exp]
        fdr = row['pval_fdr'] < 0.05

        if outcome == 'egfr':
            est, lo, hi = row['beta'], row['ci_lo'], row['ci_hi']
        else:
            est = np.exp(row['beta'])
            lo  = max(np.exp(row['ci_lo']), XMIN * 1.005)
            hi  = min(np.exp(row['ci_hi']), XMAX * 0.998)

        ax.plot([lo, hi], [y, y], color=col, lw=2.2,
                solid_capstyle='round', zorder=3)
        cap = 0.18
        ax.plot([lo, lo], [y-cap, y+cap], color=col, lw=2.0, zorder=3)
        ax.plot([hi, hi], [y-cap, y+cap], color=col, lw=2.0, zorder=3)
        mk = 'D' if fdr else 'o'
        ms = 7.5 if fdr else 6.5
        ax.plot(est, y, marker=mk, ms=ms, color=col,
                markeredgecolor='white', markeredgewidth=0.9, zorder=5)

    # Axes styling
    for spine in ['top','right','left']:
        ax.spines[spine].set_visible(False)
    ax.tick_params(left=False)
    ax.set_yticks([])
    ax.set_xlabel(xlabel, fontsize=FONT, labelpad=5)
    ax.tick_params(axis='x', labelsize=FONT - 0.5)

    if outcome == 'egfr':
        ax.set_xticks([-0.04, -0.02, 0.0, 0.02, 0.04])
        ax.set_xticklabels(['-0.040','-0.020','0.000','+0.020','+0.040'])
    else:
        ax.set_xscale('log')
        ax.set_xlim(0.30, 2.20)
        ax.xaxis.set_major_locator(FixedLocator([0.3, 0.5, 1.0, 1.5, 2.0]))
        ax.xaxis.set_major_formatter(FixedFormatter(['0.30','0.50','1.00','1.50','2.00']))

    # Panel label + title (above axes, in axes-fraction coords)
    ax.text(-0.02, 1.04, label, transform=ax.transAxes,
            fontsize=12, fontweight='bold', va='bottom', ha='right')
    ax.text(0.0, 1.04, title, transform=ax.transAxes,
            fontsize=FONT + 0.5, va='bottom', color='#222')

    return Y

# ── Draw text columns using figure-level text ─────────────────────────────────
def draw_text_cols(fig, ax, df, outcome):
    """
    Place all text using fig.text() with figure-fraction coordinates.
    This is robust to figure size changes.
    """
    ax_bbox = ax.get_position()

    # Figure-fraction x for left columns (right-justified to their slots)
    # Slots: [DAL role | Mean F | N | Exposure name] ... [forest] ... [est | CI | P | P_FDR]
    # We define column centers as figure fractions
    cx = {
        'dal':      ax_bbox.x0 - 0.155,
        'meanf':    ax_bbox.x0 - 0.105,
        'n':        ax_bbox.x0 - 0.072,
        'exp':      ax_bbox.x0 - 0.006,
        'est':      ax_bbox.x1 + 0.028,
        'ci':       ax_bbox.x1 + 0.095,
        'p':        ax_bbox.x1 + 0.158,
        'pfdr':     ax_bbox.x1 + 0.198,
    }

    ylim = ax.get_ylim()

    def fig_y(y_data):
        frac = (y_data - ylim[0]) / (ylim[1] - ylim[0])
        return ax_bbox.y0 + frac * ax_bbox.height

    for exp, y in zip(EXP_ORDER, Y):
        row = df.loc[exp]
        fdr  = row['pval_fdr'] < 0.05
        fw   = 'bold' if fdr else 'normal'
        ecol = '#8B0000' if fdr else '#333333'
        n    = inst.loc[exp, 'N']
        mf   = inst.loc[exp, 'MeanF']
        dal  = DAL_ROLE[exp]
        fy   = fig_y(y)

        # ── LEFT ────────────────────────────────────────────────────────────
        fig.text(cx['dal'],  fy, dal,          ha='center', va='center',
                 fontsize=FONT - 1.5, color=DAL_COLOR[dal], transform=fig.transFigure)
        fig.text(cx['meanf'],fy, f'{mf:.1f}',  ha='center', va='center',
                 fontsize=FONT - 0.5, color='#333', transform=fig.transFigure)
        fig.text(cx['n'],    fy, str(n),        ha='center', va='center',
                 fontsize=FONT - 0.5, color='#333', transform=fig.transFigure)
        fig.text(cx['exp'],  fy, EXP_SHORT[exp],ha='right',  va='center',
                 fontsize=FONT,     color=COLORS[exp], fontweight=fw,
                 transform=fig.transFigure)

        # ── RIGHT ────────────────────────────────────────────────────────────
        if outcome == 'egfr':
            est_s = f"{row['beta']:+.4f}"
            ci_s  = f"[{row['ci_lo']:+.4f}, {row['ci_hi']:+.4f}]"
        else:
            est_s = f"{np.exp(row['beta']):.3f}"
            ci_s  = f"[{np.exp(row['ci_lo']):.3f}, {np.exp(row['ci_hi']):.3f}]"
        p_s   = fmt_p(row['pval'])
        fdr_s = fmt_p(row['pval_fdr'])
        star  = ' *' if fdr else ''

        fig.text(cx['est'],  fy, est_s,          ha='left',   va='center',
                 fontsize=FONT, fontweight=fw, color=ecol, transform=fig.transFigure)
        fig.text(cx['ci'],   fy, ci_s,            ha='center', va='center',
                 fontsize=FONT - 0.5, color='#333', transform=fig.transFigure)
        fig.text(cx['p'],    fy, p_s + star,      ha='center', va='center',
                 fontsize=FONT, fontweight=fw, color=ecol, transform=fig.transFigure)
        fig.text(cx['pfdr'], fy, fdr_s + star,    ha='center', va='center',
                 fontsize=FONT, fontweight=fw, color=ecol, transform=fig.transFigure)

    # ── HEADERS ─────────────────────────────────────────────────────────────
    hy = fig_y(N - 0.1)  # just above top row
    headers_l = [
        (cx['dal'],   'DAL role',  'center'),
        (cx['meanf'], 'Mean F',   'center'),
        (cx['n'],     'N',         'center'),
        (cx['exp'],   'Exposure',  'right'),
    ]
    headers_r = [
        (cx['est'],   'β' if outcome=='egfr' else 'OR', 'left'),
        (cx['ci'],    '95% CI',   'center'),
        (cx['p'],     'P',         'center'),
        (cx['pfdr'],  'P_FDR',    'center'),
    ]
    for xf, lbl, ha in headers_l + headers_r:
        fig.text(xf, hy, lbl, ha=ha, va='bottom',
                 fontsize=FONT, fontweight='bold', color='#111',
                 transform=fig.transFigure)

    # Separator line just below headers
    line_y = fig_y(N - 0.38)
    line = plt.Line2D([ax_bbox.x0 - 0.20, ax_bbox.x1 + 0.21],
                       [line_y, line_y],
                       transform=fig.transFigure, color='#AAAAAA', lw=0.8,
                       clip_on=False)
    fig.add_artist(line)


# ── Build ─────────────────────────────────────────────────────────────────────
draw_forest(ax_e, egfr_df, 'egfr', 'A',
            'Dietary exposures → eGFR  [CKDGen, N = 480,698, European ancestry]')
draw_forest(ax_c, ckd_df,  'ckd',  'B',
            'Dietary exposures → CKD (binary)  [CKDGen; 41,395 cases / 439,303 controls]')

# Text columns need axes to be drawn first
fig.canvas.draw()
draw_text_cols(fig, ax_e, egfr_df, 'egfr')
draw_text_cols(fig, ax_c, ckd_df,  'ckd')

# ── Legend (bottom) ───────────────────────────────────────────────────────────
legend_handles = [
    Line2D([0],[0], marker='D', color='#555', ms=7, lw=0,
           markeredgecolor='white', markeredgewidth=0.8, label='FDR < 0.05  (◆)'),
    Line2D([0],[0], marker='o', color='#888', ms=6, lw=0,
           markeredgecolor='white', markeredgewidth=0.8, label='Not significant  (●)'),
    mpatches.Patch(color='#2E8B57', alpha=0.85, label='Base-forming / weakly basic'),
    mpatches.Patch(color='#B8860B', alpha=0.85, label='Mild acid-forming'),
    mpatches.Patch(color='#C0392B', alpha=0.85, label='Acid-forming'),
]
fig.legend(handles=legend_handles, loc='lower center', ncol=5,
           fontsize=FONT - 0.5, frameon=True,
           bbox_to_anchor=(0.5, 0.005),
           fancybox=False, edgecolor='#CCC',
           handlelength=1.4, columnspacing=1.0)

# ── Save ──────────────────────────────────────────────────────────────────────
out = 'figures/Figure1_forest_extended.png'
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Saved: {out}")
plt.close()
