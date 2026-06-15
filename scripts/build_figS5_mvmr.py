"""
FigureS5 – MVMR Extended Forest Plot
Rebuilt with robust figure-level text placement (same approach as Figure 1).
Design changes vs. original:
  - No top title/legend line
  - Left whitespace reduced
  - Plot widths narrowed
  - All text via fig.text() with figure-fraction coords
"""
import pandas as pd, numpy as np, math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

proc = "data/processed/"
figs = "figures/"

mvmr  = pd.read_csv(proc + "mr_results_mvmr.csv")
fstat = pd.read_csv(proc + "mvmr_conditional_fstats.csv").set_index('exposure')

EXP = ["Fresh fruit intake","Oily fish intake","Cooked vegetable intake",
       "Salad/raw vegetable intake","Processed meat intake","Beef intake"]
SHORT = {"Fresh fruit intake":"Fresh fruit","Oily fish intake":"Oily fish",
         "Cooked vegetable intake":"Cooked veg","Salad/raw vegetable intake":"Salad/raw veg",
         "Processed meat intake":"Processed meat","Beef intake":"Beef"}

# DAL-role colors (consistent with Figure 1)
_DAL = {
    'Fresh fruit intake':          '#2E8B57',
    'Oily fish intake':            '#C0392B',
    'Cooked vegetable intake':     '#2E8B57',
    'Salad/raw vegetable intake':  '#2E8B57',
    'Processed meat intake':       '#C0392B',
    'Beef intake':                 '#C0392B',
}
def exp_col(exp):
    return _DAL[exp]

mvmr_e = mvmr[mvmr['outcome']=='eGFR (log-transformed)'].set_index('exposure')
mvmr_c = mvmr[mvmr['outcome']=='CKD (binary)'].set_index('exposure')

N = len(EXP)

# ── Figure ─────────────────────────────────────────────────────────────────────
FIG_W, FIG_H = 16.54, 6.5
FONT = 8.5
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H), facecolor='white')
fig.subplots_adjust(top=0.88, bottom=0.10, left=0.00, right=0.97)
ax.set_xlim(0,1); ax.set_ylim(0,1)
for s in ax.spines.values(): s.set_visible(False)
ax.set_xticks([]); ax.set_yticks([])
ax.set_facecolor('white')

# ── Layout constants ───────────────────────────────────────────────────────────
# Left text: Exposure | F_cond | N
# Forest eGFR:  FEL → FER
# Forest CKD:   FCL → FCR
# Right text: β CI | OR CI
CX = dict(
    exp   = 0.095,   # right-align exposure name  ← weniger links-Whitespace
    fcond = 0.135,   # F_cond
    n     = 0.163,   # N
    FEL   = 0.178, FER = 0.395,   # eGFR forest
    FCL   = 0.420, FCR = 0.628,   # CKD forest
    be    = 0.644,   # β [CI] merged
    or_   = 0.790,   # OR [CI] merged  ← enger zusammen
)

HDR_Y = 0.960
def row_y(i): return 0.88 - i * (0.88 - 0.08) / (N - 1)
ROW_H = (0.88 - 0.08) / (N - 1)

# ── Data range mappers ─────────────────────────────────────────────────────────
EFXLO, EFXHI = -0.038, 0.050
CFXLO, CFXHI = np.log(0.35), np.log(2.10)

def fe2a(v): return CX['FEL'] + np.clip((v-EFXLO)/(EFXHI-EFXLO),0,1)*(CX['FER']-CX['FEL'])
def fc2a(v): return CX['FCL'] + np.clip((v-CFXLO)/(CFXHI-CFXLO),0,1)*(CX['FCR']-CX['FCL'])

# ── Row shading ────────────────────────────────────────────────────────────────
cap = ROW_H * 0.42
for i in range(N):
    yc = row_y(i)
    ax.fill_between([0,1],[yc-cap*1.05],[yc+cap*1.05],
                    color='#F0F2F5' if i%2==0 else 'white',
                    zorder=0, transform=ax.transAxes)

# ── Null lines ─────────────────────────────────────────────────────────────────
for x in [fe2a(0.0), fc2a(0.0)]:
    ax.plot([x,x],[0, 0.91], color='#444', lw=1.1, ls='--',
            transform=ax.transAxes, zorder=2)
for xref in [fe2a(-0.02), fe2a(0.02)]:
    ax.plot([xref,xref],[0,1], color='#DDD', lw=0.7, transform=ax.transAxes, zorder=0)
for xref in [fc2a(np.log(0.5)), fc2a(np.log(2.0))]:
    ax.plot([xref,xref],[0,1], color='#DDD', lw=0.7, transform=ax.transAxes, zorder=0)

# ── Panel separator ────────────────────────────────────────────────────────────
xmid = (CX['FER'] + CX['FCL']) / 2
ax.plot([xmid,xmid],[0.03,0.97], color='#CCCCCC', lw=1.0,
        transform=ax.transAxes, zorder=1)

# ── X-axes ─────────────────────────────────────────────────────────────────────
y_ax = 0.048
# eGFR
ax.plot([CX['FEL'],CX['FER']],[y_ax]*2, color='#444', lw=0.9,
        transform=ax.transAxes, clip_on=False)
for td, tl in zip([-0.03,-0.015,0.0,0.015,0.03],
                  ['-0.030','-0.015','0.000','+0.015','+0.030']):
    xf = fe2a(td)
    ax.plot([xf,xf],[y_ax,y_ax-0.018], color='#444', lw=0.8,
            transform=ax.transAxes, clip_on=False)
    ax.text(xf, y_ax-0.028, tl, transform=ax.transAxes,
            ha='center', va='top', fontsize=7.5, color='#333', clip_on=False)
ax.text((CX['FEL']+CX['FER'])/2, y_ax-0.065, 'β per SD increase  [95% CI]',
        transform=ax.transAxes, ha='center', va='top', fontsize=FONT, color='#333',
        clip_on=False)

# CKD
ax.plot([CX['FCL'],CX['FCR']],[y_ax]*2, color='#444', lw=0.9,
        transform=ax.transAxes, clip_on=False)
for td, tl in zip([np.log(x) for x in [0.35,0.5,1.0,1.5,2.0]],
                  ['0.35','0.50','1.00','1.50','2.00']):
    xf = fc2a(td)
    ax.plot([xf,xf],[y_ax,y_ax-0.018], color='#444', lw=0.8,
            transform=ax.transAxes, clip_on=False)
    ax.text(xf, y_ax-0.028, tl, transform=ax.transAxes,
            ha='center', va='top', fontsize=7.5, color='#333', clip_on=False)
ax.text((CX['FCL']+CX['FCR'])/2, y_ax-0.065, 'Odds Ratio for CKD  [95% CI]',
        transform=ax.transAxes, ha='center', va='top', fontsize=FONT, color='#333',
        clip_on=False)

# ── Column headers ─────────────────────────────────────────────────────────────
hdrs = [
    (CX['exp'],                      'Exposure',           'right'),
    (CX['fcond'],                    'F_cond',             'center'),
    (CX['n'],                        'N',                  'center'),
    ((CX['FEL']+CX['FER'])/2,        'eGFR outcome',       'center'),
    ((CX['FCL']+CX['FCR'])/2,        'CKD outcome',        'center'),
    (CX['be'],                       'β (eGFR)  [95% CI]', 'left'),
    (CX['or_'],                      'OR (CKD)  [95% CI]', 'left'),
]
for xf, lbl, ha in hdrs:
    ax.text(xf, HDR_Y, lbl, transform=ax.transAxes, ha=ha, va='bottom',
            fontsize=FONT, fontweight='bold', color='#111', clip_on=False)

# Bottom divider only (below last data row)
bot_div = row_y(N-1) - ROW_H * 0.58
ax.plot([0,1],[bot_div]*2, color='#AAAAAA', lw=0.9, transform=ax.transAxes)

# ── Data rows ──────────────────────────────────────────────────────────────────
for i, exp in enumerate(EXP):
    yc   = row_y(i)
    fc   = fstat.loc[exp,'F_conditional']
    weak = fc < 10.0
    col  = exp_col(exp)
    re   = mvmr_e.loc[exp]; rc = mvmr_c.loc[exp]
    n    = int(re['nsnp'])

    ls  = '--' if weak else '-'
    alp = 0.65 if weak else 1.0
    lw  = 1.4  if weak else 2.0

    # eGFR
    est_e, lo_e, hi_e = re['beta'], re['ci_lo'], re['ci_hi']
    sig_e = re['pval_fdr'] < 0.05
    lo_ea, hi_ea, est_ea = fe2a(lo_e), fe2a(hi_e), fe2a(est_e)
    ax.plot([lo_ea,hi_ea],[yc,yc], color=col, lw=lw, alpha=alp,
            ls=ls, transform=ax.transAxes, zorder=3, clip_on=False)
    for xc in [lo_ea, hi_ea]:
        ax.plot([xc,xc],[yc-cap*.9,yc+cap*.9], color=col, lw=lw-0.3,
                alpha=alp, transform=ax.transAxes, zorder=3, clip_on=False)
    ax.plot(est_ea, yc, 'D' if sig_e else 'o', ms=7, color=col,
            alpha=alp, markeredgecolor='white', markeredgewidth=0.8,
            transform=ax.transAxes, zorder=5, clip_on=False)

    # CKD
    est_c, lo_c, hi_c = rc['beta'], rc['ci_lo'], rc['ci_hi']
    sig_c = rc['pval_fdr'] < 0.05
    lo_ca, hi_ca, est_ca = fc2a(lo_c), fc2a(hi_c), fc2a(est_c)
    ax.plot([lo_ca,hi_ca],[yc,yc], color=col, lw=lw, alpha=alp,
            ls=ls, transform=ax.transAxes, zorder=3, clip_on=False)
    for xc in [lo_ca, hi_ca]:
        ax.plot([xc,xc],[yc-cap*.9,yc+cap*.9], color=col, lw=lw-0.3,
                alpha=alp, transform=ax.transAxes, zorder=3, clip_on=False)
    ax.plot(est_ca, yc, 'D' if sig_c else 'o', ms=7, color=col,
            alpha=alp, markeredgecolor='white', markeredgewidth=0.8,
            transform=ax.transAxes, zorder=5, clip_on=False)

    # ── Text ──────────────────────────────────────────────────────────────────
    fw = 'bold' if (sig_e or sig_c) else 'normal'
    # Use ax.transAxes for intra-axes text (these are stable since ax spans full figure)
    ax.text(CX['exp'],   yc, SHORT[exp], transform=ax.transAxes,
            ha='right', va='center', fontsize=FONT, fontweight=fw, color=col, clip_on=False)
    ax.text(CX['fcond'], yc, f"{fc:.1f}", transform=ax.transAxes,
            ha='center', va='center', fontsize=FONT-0.5,
            color='#1A6A1A' if fc>=10 else '#AA5500',
            fontweight='bold' if fc>=10 else 'normal', clip_on=False)
    ax.text(CX['n'],     yc, str(n), transform=ax.transAxes,
            ha='center', va='center', fontsize=FONT-0.5, color='#444', clip_on=False)

    # eGFR numeric
    or_v  = math.exp(est_c); or_lo = math.exp(lo_c); or_hi = math.exp(hi_c)
    pcol_e = '#990000' if sig_e else '#333'
    pcol_c = '#990000' if sig_c else '#333'
    # Merged β [CI] and OR [CI]
    ax.text(CX['be'], yc,
            f"{est_e:+.4f}  [{lo_e:+.4f}, {hi_e:+.4f}]",
            transform=ax.transAxes, ha='left', va='center',
            fontsize=FONT-1.0, fontweight=fw, color=pcol_e, clip_on=False)
    ax.text(CX['or_'], yc,
            f"{or_v:.3f}  [{or_lo:.3f}, {or_hi:.3f}]",
            transform=ax.transAxes, ha='left', va='center',
            fontsize=FONT-1.0, fontweight=fw, color=pcol_c, clip_on=False)

# ── Legend (bottom only, no top title) ────────────────────────────────────────
leg = [
    Line2D([0],[0], color='#2E8B57', lw=2.0,
           label='Base-forming (F_cond \u2265 10: adequate)'),
    Line2D([0],[0], color='#C0392B', lw=1.5, ls='--', alpha=0.7,
           label='Acid-forming (F_cond < 10: weak)'),
    Line2D([0],[0], marker='o', ms=6, lw=0, color='#555',
           markeredgecolor='white', label='Not significant (all exposures)'),
]
fig.legend(handles=leg, loc='lower center', ncol=3, fontsize=FONT-0.5,
           frameon=True, bbox_to_anchor=(0.5, 0.002),
           fancybox=False, edgecolor='#CCC', framealpha=0.95,
           handlelength=1.8, columnspacing=1.2)

# ── Save ───────────────────────────────────────────────────────────────────────
out = figs + 'FigureS5_MVMR_extended.png'
fig.savefig(out, dpi=300, bbox_inches='tight', facecolor='white', pad_inches=0.12)
print('Saved:', out)
plt.close()
