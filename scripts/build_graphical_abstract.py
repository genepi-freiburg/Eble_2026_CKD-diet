#!/usr/bin/env python3
"""
Build a richly-structured graphical abstract for the
Mendelian Randomization study (dietary components -> kidney function).

Layout (top -> bottom):
  1. Title block
  2. Two-column header: Exposure GWAS  |  Outcome GWAS
  3. Diamond: Instrument selection / filtering / harmonisation
  4. Filter-cascade strip with exact SNP/pair counts
  5. Discovery block with three sub-pipelines:
       (1) Univariable MR  (2) Multivariable MR  (3) GI-PRAL composite
  6. Sensitivity-analyses block (pleiotropy/robustness, eGFRcys, GI-PRAL LOO)
  7. Validation bar: FinnGen R12 replication

Outputs:
  figures/graphical_abstract.png  (300 dpi)
  figures/graphical_abstract.pdf  (vector)

Run from repo root:
  python3 scripts/build_graphical_abstract.py
"""

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon, Rectangle

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "figures"
OUT_DIR.mkdir(exist_ok=True)

# --- Palette ------------------------------------------------------------------
C_BG          = "#EAF2F6"   # canvas background (light blue)
C_EXP_DARK    = "#1F6F73"   # exposure header (teal)
C_EXP_LIGHT   = "#DCEBEB"
C_OUT_DARK    = "#4F6B2A"   # outcome header (olive)
C_OUT_LIGHT   = "#E5ECDA"
C_DIAMOND     = "#7A4B91"   # purple diamond (filtering)
C_DIAMOND_LT  = "#E9DDF2"
C_STRIP       = "#D6D2DC"   # filter cascade strip
C_DISC_DARK   = "#26365A"   # discovery navy header
C_DISC_LIGHT  = "#E4E6EC"
C_SENS_DARK   = "#8E4A2F"   # sensitivity terracotta
C_SENS_LIGHT  = "#F4E1D6"
C_VAL_DARK    = "#D4791F"   # validation orange
C_VAL_LIGHT   = "#FBE6CE"
C_TEXT        = "#1A1A1A"
C_TEXT_LIGHT  = "#FFFFFF"
C_ARROW       = "#3A3A3A"

# --- Canvas -------------------------------------------------------------------
FIG_W, FIG_H = 8.5, 18.0
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_aspect("auto")
ax.axis("off")
# Background
ax.add_patch(Rectangle((0, 0), 100, 100, facecolor=C_BG, zorder=-10))

# --- Helpers ------------------------------------------------------------------
def rbox(x, y, w, h, *, fc, ec=None, lw=1.0, rounding=1.4, zorder=1):
    """Rounded rectangle."""
    if ec is None:
        ec = fc
    p = FancyBboxPatch(
        (x, y), w, h,
        boxstyle=f"round,pad=0.02,rounding_size={rounding}",
        linewidth=lw, edgecolor=ec, facecolor=fc, zorder=zorder,
    )
    ax.add_patch(p)
    return p

def arrow(x1, y1, x2, y2, *, color=C_ARROW, lw=1.8, mscale=18, zorder=3):
    a = FancyArrowPatch(
        (x1, y1), (x2, y2),
        arrowstyle="-|>", mutation_scale=mscale,
        linewidth=lw, color=color, zorder=zorder,
    )
    ax.add_patch(a)

def text(x, y, s, *, size=10, weight="normal", color=C_TEXT,
         ha="center", va="center", style="normal", zorder=5):
    ax.text(x, y, s, fontsize=size, fontweight=weight, color=color,
            ha=ha, va=va, fontstyle=style, zorder=zorder)

# =============================================================================
# 1) Title
# =============================================================================
text(50, 97.8,
     "Causal Effects of Dietary Components on Kidney Function",
     size=16.5, weight="bold")
text(50, 95.8,
     "Two-Sample Mendelian Randomization · Study Flowchart",
     size=11.0, style="italic", color="#333333")

# =============================================================================
# 2) Two-column header: Exposure | Outcome
# =============================================================================
HDR_TOP  = 94.0
HDR_BOT  = 76.0   # taller boxes (height 18) for 6 body lines + title
GAP_COL  = 3.0
COL_W    = (100 - 2*4 - GAP_COL) / 2   # margins 4 on each side
LX = 4
RX = LX + COL_W + GAP_COL

# --- Exposure column ---------------------------------------------------------
rbox(LX, HDR_TOP - 2.6, COL_W, 2.6, fc=C_EXP_DARK, ec=C_EXP_DARK, rounding=1.0)
text(LX + COL_W/2, HDR_TOP - 1.3, "EXPOSURE GWAS",
     size=10.5, weight="bold", color=C_TEXT_LIGHT)
rbox(LX, HDR_BOT, COL_W, HDR_TOP - 2.6 - HDR_BOT, fc=C_EXP_LIGHT,
     ec=C_EXP_DARK, lw=1.2, rounding=1.2)
exp_body_x = LX + 2.5
exp_body_y = HDR_TOP - 4.5
text(exp_body_x, exp_body_y, "10 dietary exposures",
     size=10.5, weight="bold", ha="left", va="top")
exp_lines = [
    "(UK Biobank, European ancestry)",
    "MRC-IEU GWAS pipeline (ukb-b-*)",
    "N = 421,764 – 461,981 per exposure",
    "",
    "Fruits · vegetables · tea · oily fish",
    "processed meat · lamb · pork · beef",
]
for i, line in enumerate(exp_lines):
    text(exp_body_x, exp_body_y - 1.9 - 1.55*i, line,
         size=8.7, ha="left", va="top")

# --- Outcome column ----------------------------------------------------------
rbox(RX, HDR_TOP - 2.6, COL_W, 2.6, fc=C_OUT_DARK, ec=C_OUT_DARK, rounding=1.0)
text(RX + COL_W/2, HDR_TOP - 1.3, "OUTCOME GWAS",
     size=10.5, weight="bold", color=C_TEXT_LIGHT)
rbox(RX, HDR_BOT, COL_W, HDR_TOP - 2.6 - HDR_BOT, fc=C_OUT_LIGHT,
     ec=C_OUT_DARK, lw=1.2, rounding=1.2)
out_body_x = RX + 2.5
out_body_y = HDR_TOP - 4.5
text(out_body_x, out_body_y, "Primary: CKDGen 2019",
     size=10.5, weight="bold", ha="left", va="top")
out_lines = [
    "Wuttke et al. — European ancestry",
    "eGFRcrea  (N = 567,460)",
    "CKD       (N = 480,698)",
    "",
    "Replication: FinnGen R12",
    "N = 493,235  (12,787 cases)",
]
for i, line in enumerate(out_lines):
    text(out_body_x, out_body_y - 1.9 - 1.55*i, line,
         size=8.7, ha="left", va="top")

# =============================================================================
# 3) Diamond: Instrument selection
# =============================================================================
DIA_CY = 70.0
DIA_HW = 20.0   # half-width (wider so text fits inside)
DIA_HH = 6.5    # half-height (taller so 3 text lines fit comfortably)
diamond_pts = [
    (50,            DIA_CY + DIA_HH),
    (50 + DIA_HW,   DIA_CY),
    (50,            DIA_CY - DIA_HH),
    (50 - DIA_HW,   DIA_CY),
]
ax.add_patch(Polygon(diamond_pts, closed=True,
                     facecolor=C_DIAMOND, edgecolor=C_DIAMOND, linewidth=1.0, zorder=2))
text(50, DIA_CY + 2.0, "Instrument selection,",
     size=9.0, weight="bold", color=C_TEXT_LIGHT, zorder=6)
text(50, DIA_CY,       "filtering &",
     size=9.0, weight="bold", color=C_TEXT_LIGHT, zorder=6)
text(50, DIA_CY - 2.0, "harmonisation",
     size=9.0, weight="bold", color=C_TEXT_LIGHT, zorder=6)

# Arrows from header columns into diamond
edge_t = 0.30
NW = (50 - DIA_HW + edge_t*(DIA_HW),            DIA_CY + edge_t*DIA_HH)
NE = (50 + DIA_HW - edge_t*(DIA_HW),            DIA_CY + edge_t*DIA_HH)
arrow(LX + COL_W/2,  HDR_BOT - 0.3,  NW[0], NW[1],
      color=C_EXP_DARK, lw=1.8)
arrow(RX + COL_W/2,  HDR_BOT - 0.3,  NE[0], NE[1],
      color=C_OUT_DARK, lw=1.8)

# =============================================================================
# 4) Filter-cascade strip
# =============================================================================
STRIP_TOP = 61.5
STRIP_H   = 7.0
STRIP_X   = 4
STRIP_W   = 92
rbox(STRIP_X, STRIP_TOP - STRIP_H, STRIP_W, STRIP_H,
     fc=C_STRIP, ec="#9C99A6", lw=0.8, rounding=1.2)
text(50, STRIP_TOP - 2.3,
     "326 SNPs  →  PheWAS screening (TG · BMI · T2D · SBP · CAD)  →  −46 SNPs",
     size=9.5, weight="bold")
text(50, STRIP_TOP - 4.8,
     "280 instruments  →  harmonisation (−13)  →  267  →  Steiger filtering (−3)  →  "
     "264 SNP–exposure pairs",
     size=9.5, weight="bold")

# Arrow from diamond bottom to strip
arrow(50, DIA_CY - DIA_HH - 0.3, 50, STRIP_TOP + 0.3,
      color=C_DIAMOND, lw=2.0, mscale=20)

# =============================================================================
# 5) Discovery block
# =============================================================================
DISC_TOP  = STRIP_TOP - STRIP_H - 2.0   # ~ 52.5
DISC_BOT  = 29.0
DISC_X    = 4
DISC_W    = 92
DISC_H    = DISC_TOP - DISC_BOT

# Outer container (navy edge)
rbox(DISC_X, DISC_BOT, DISC_W, DISC_H,
     fc="#F4F5F8", ec=C_DISC_DARK, lw=1.6, rounding=1.4)
# Dark header strip
rbox(DISC_X, DISC_TOP - 3.0, DISC_W, 3.0,
     fc=C_DISC_DARK, ec=C_DISC_DARK, rounding=1.0)
text(50, DISC_TOP - 1.5, "DISCOVERY",
     size=12.0, weight="bold", color=C_TEXT_LIGHT)

# Arrow from strip down into Discovery
arrow(50, STRIP_TOP - STRIP_H - 0.3, 50, DISC_TOP + 0.3,
      color=C_DISC_DARK, lw=2.0, mscale=20)

# Three sub-panels
SUB_PAD   = 2.0
SUB_TOP   = DISC_TOP - 3.0 - SUB_PAD
SUB_BOT   = DISC_BOT + SUB_PAD
SUB_W     = (DISC_W - 4*SUB_PAD) / 3
SUB_H     = SUB_TOP - SUB_BOT

sub_x = [DISC_X + SUB_PAD + i*(SUB_W + SUB_PAD) for i in range(3)]

def sub_panel(i, num, title_lines, body_lines):
    x = sub_x[i]
    rbox(x, SUB_BOT, SUB_W, SUB_H,
         fc=C_DISC_LIGHT, ec="#8C95AE", lw=0.9, rounding=1.2)
    cy = SUB_TOP - 2.0
    ax.add_patch(plt.Circle((x + 3.0, cy), 1.3,
                            facecolor=C_DISC_DARK, edgecolor=C_DISC_DARK, zorder=4))
    text(x + 3.0, cy, str(num), size=9.0, weight="bold",
         color=C_TEXT_LIGHT, zorder=5)
    for j, tl in enumerate(title_lines):
        text(x + SUB_W/2 + 2.0, cy - 0.1 - 1.35*j, tl,
             size=9.2, weight="bold", ha="center")
    body_y0 = cy - 1.35*len(title_lines) - 1.4
    for j, bl in enumerate(body_lines):
        text(x + 2.0, body_y0 - 1.35*j, bl,
             size=7.8, ha="left", va="top")

sub_panel(
    0, 1,
    ["Univariable", "two-sample MR"],
    [
        "• IVW (primary)",
        "• Benjamini–Hochberg",
        "  FDR correction",
        "• 10 exposures ×",
        "  2 outcomes (eGFR, CKD)",
        "• Steiger-filtered",
        "  264 SNP–exposure pairs",
    ],
)
sub_panel(
    1, 2,
    ["Multivariable MR", "(co-exposures)"],
    [
        "• MVMR-IVW",
        "• Sanderson–Windmeijer",
        "  conditional F-statistic",
        "• Instrument-strength",
        "  check per exposure",
        "• Direct vs total",
        "  effect decomposition",
    ],
)
sub_panel(
    2, 3,
    ["GI-PRAL", "composite"],
    [
        "• Dietary acid load (DAL)",
        "  omnibus test vs eGFR",
        "• PRAL weights from",
        "  Remer & Manz 1995",
        "• Weighted combination",
        "  of 10 acid-load foods",
        "• β = −0.324,  P = 0.032",
    ],
)

# =============================================================================
# 6) Sensitivity Analyses block
# =============================================================================
SENS_TOP  = DISC_BOT - 2.0     # ~ 27.0
SENS_BOT  = 6.5
SENS_X    = 4
SENS_W    = 92
SENS_H    = SENS_TOP - SENS_BOT

# Outer container
rbox(SENS_X, SENS_BOT, SENS_W, SENS_H,
     fc="#FBF1EA", ec=C_SENS_DARK, lw=1.6, rounding=1.4)
# Dark header strip
rbox(SENS_X, SENS_TOP - 3.0, SENS_W, 3.0,
     fc=C_SENS_DARK, ec=C_SENS_DARK, rounding=1.0)
text(50, SENS_TOP - 1.5, "SENSITIVITY ANALYSES",
     size=12.0, weight="bold", color=C_TEXT_LIGHT)

# Arrow from Discovery into Sensitivity
arrow(50, DISC_BOT - 0.3, 50, SENS_TOP + 0.3,
      color=C_SENS_DARK, lw=2.0, mscale=20)

# Two sub-panels (Pleiotropy/robustness | Outcome- and food-level checks)
SENS_PAD = 2.0
SENS_SUB_TOP = SENS_TOP - 3.0 - SENS_PAD
SENS_SUB_BOT = SENS_BOT + SENS_PAD
SENS_SUB_W   = (SENS_W - 3*SENS_PAD) / 2
SENS_SUB_H   = SENS_SUB_TOP - SENS_SUB_BOT
sens_x = [SENS_X + SENS_PAD + i*(SENS_SUB_W + SENS_PAD) for i in range(2)]

def sens_panel(i, title, body_lines):
    x = sens_x[i]
    rbox(x, SENS_SUB_BOT, SENS_SUB_W, SENS_SUB_H,
         fc=C_SENS_LIGHT, ec="#B97B5C", lw=0.9, rounding=1.2)
    text(x + SENS_SUB_W/2, SENS_SUB_TOP - 1.6, title,
         size=9.6, weight="bold")
    for j, bl in enumerate(body_lines):
        text(x + 2.0, SENS_SUB_TOP - 3.4 - 1.55*j, bl,
             size=8.0, ha="left", va="top")

sens_panel(
    0,
    "Pleiotropy & robustness",
    [
        "• Weighted median  ·  Weighted mode",
        "• MR-Egger  (+ intercept test)",
        "• Cochran's Q  (heterogeneity)",
        "• Steiger filtering",
        "  (3 SNP–exposure pairs removed)",
        "• Leave-one-out  (Tea, Salad,",
        "  Dried fruit)  ×  eGFR & CKD",
    ],
)
sens_panel(
    1,
    "Alternative outcomes & composite LOO",
    [
        "• eGFRcys  (Pattaro 2016)",
        "  ieu-a-1106,  N = 33,152",
        "  116 / 264 SNPs;  all P > 0.20",
        "",
        "• GI-PRAL leave-one-out",
        "  Full:  β = −0.324,  P = 0.032",
        "  w/o Pork:  β = −0.179,  P = 0.132",
    ],
)

# =============================================================================
# 7) Validation bar
# =============================================================================
VAL_TOP = 5.5
VAL_BOT = 0.5
VAL_X   = 4
VAL_W   = 92
rbox(VAL_X, VAL_BOT, VAL_W, VAL_TOP - VAL_BOT,
     fc=C_VAL_LIGHT, ec=C_VAL_DARK, lw=1.6, rounding=1.2)
# Dark header strip
rbox(VAL_X, VAL_TOP - 3.0, VAL_W, 3.0,
     fc=C_VAL_DARK, ec=C_VAL_DARK, rounding=1.0)
text(50, VAL_TOP - 1.5,
     "VALIDATION — FinnGen R12 Replication",
     size=11.5, weight="bold", color=C_TEXT_LIGHT)

text(50, VAL_TOP - 3.5,
     "Same IVW + sensitivity pipeline applied independently  ·  "
     "Replicated = FDR < 0.05 AND directionally consistent",
     size=8.4, ha="center")

# Arrow from Sensitivity into Validation
arrow(50, SENS_BOT - 0.3, 50, VAL_TOP + 0.3,
      color=C_VAL_DARK, lw=2.0, mscale=20)

# (Footer note removed — folded into validation strip above)

# --- Save --------------------------------------------------------------------
fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
png_path = OUT_DIR / "graphical_abstract.png"
pdf_path = OUT_DIR / "graphical_abstract.pdf"
fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor=C_BG)
fig.savefig(pdf_path,           bbox_inches="tight", facecolor=C_BG)
plt.close(fig)

print(f"Wrote: {png_path}")
print(f"Wrote: {pdf_path}")
