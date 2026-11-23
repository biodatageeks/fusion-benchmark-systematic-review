"""PRISMA 2020 flow diagram (Figure 2).

Database branch: counts from Supplementary Table S2 (PubMed and Scopus search).
Other-methods branch: benchmarks identified by backward and forward citation searching
of the included benchmarks (Creason 2021, Liu 2016, Tan 2015, Rodriguez-Martin 2017,
Jia 2013, Benelli 2012, Nicorici 2014). Four of them (Jia 2013, Tan 2015, Liu 2016,
Rodriguez-Martin 2017) are also retrieved by the database queries of Supplementary Table S2
but had not been selected from the database records; they were included after full-text
assessment.

Outputs: results/prisma/prisma_flow.png and .pdf
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "prisma"

FILL = "#8e60a0"
EDGE = "#775086"
TEXT = "white"

# Database search (Supplementary Table S2)
PUBMED, SCOPUS = 647, 802
DUPLICATES = 472
SCREENED, SCREEN_EXCLUDED = 977, 826
SOUGHT, NOT_RETRIEVED = 151, 0
ASSESSED, FULLTEXT_EXCLUDED = 151, 140
EXCLUSION_REASONS = [
    ("No benchmark or validation component", 47),
    ("No RNA-seq/RNA-based NGS data", 5),
    ("Mainly DNA-seq/WGS/SV calling", 9),
    ("Review without original comparative analysis", 8),
    ("Case study without tool assessment", 31),
    ("No clear performance metrics", 19),
    ("Insufficient description of data/methods", 12),
    ("Duplicate or extended version of same analysis", 5),
    ("Other reason", 4),
]
DATABASE_INCLUDED = 11
# Backward and forward citation searching
CITATION_IDENTIFIED = CITATION_ASSESSED = CITATION_INCLUDED = 7
# of these, also retrieved by the database queries but not selected from the database records
CITATION_ALSO_IN_DATABASE = 4


def pick_font() -> str:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Arial Rounded MT Bold", "Arial", "Helvetica", "DejaVu Sans"):
        if name in available:
            return name
    return "DejaVu Sans"


WEIGHT = "bold"


def box(ax, x0, y0, x1, y1, text, size=13, weight=None, align="center", lines=None):
    weight = weight or WEIGHT
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle="round,pad=0,rounding_size=1.6",
                                facecolor=FILL, edgecolor=EDGE, linewidth=1.6))
    if align == "center":
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, text, ha="center", va="center", color=TEXT,
                fontsize=size, fontweight=weight, linespacing=1.35)
    else:  # heading plus left-aligned list
        ax.text((x0 + x1) / 2, y1 - 1.6, text, ha="center", va="top", color=TEXT, fontsize=size,
                fontweight=weight)
        ax.text(x0 + 1.2, y1 - 5.2, "\n".join(lines), ha="left", va="top", color=TEXT,
                fontsize=size - 1.5, fontweight=weight, linespacing=1.5)


def side_label(ax, y0, y1, text):
    ax.add_patch(FancyBboxPatch((1, y0), 5, y1 - y0, boxstyle="round,pad=0,rounding_size=1.2",
                                facecolor=FILL, edgecolor=EDGE, linewidth=1.6))
    ax.text(3.5, (y0 + y1) / 2, text, rotation=90, ha="center", va="center", color=TEXT,
            fontsize=15, fontweight=WEIGHT)


def arrow(ax, x0, y0, x1, y1):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14,
                                 color=EDGE, linewidth=1.6, shrinkA=0, shrinkB=0))


def main() -> None:
    global WEIGHT
    font = pick_font()
    plt.rcParams["font.family"] = font
    WEIGHT = "normal" if "Bold" in font else "bold"
    fig, ax = plt.subplots(figsize=(17, 10))
    ax.set_xlim(0, 137)
    ax.set_ylim(0, 100)
    ax.axis("off")

    main_x0, main_x1 = 9, 41
    excl_x0, excl_x1 = 47, 92
    oth_x0, oth_x1 = 100, 136
    mid = (main_x0 + main_x1) / 2
    oth_mid = (oth_x0 + oth_x1) / 2

    # Headers
    box(ax, main_x0, 91, excl_x1, 98.5, "Identification of studies via databases", size=16)
    box(ax, oth_x0, 91, oth_x1, 98.5, "Identification of studies\nvia other methods", size=14)

    # Side labels
    side_label(ax, 72, 90, "Identification")
    side_label(ax, 30, 70, "Screening")
    side_label(ax, 1.5, 28, "Inclusion")

    # Database branch
    box(ax, main_x0, 76, main_x1, 88, f"Records identified from:\nPubMed (n={PUBMED})\nScopus (n={SCOPUS})")
    box(ax, excl_x0, 77, excl_x1, 87, f"Duplicate records removed\nbefore screening (n={DUPLICATES})")
    box(ax, main_x0, 62, main_x1, 68.5, f"Records screened (n={SCREENED})")
    box(ax, excl_x0, 62, excl_x1, 68.5, f"Records excluded (n={SCREEN_EXCLUDED})")
    box(ax, main_x0, 49.5, main_x1, 57, f"Reports sought for retrieval\n(n={SOUGHT})")
    box(ax, excl_x0, 49.5, excl_x1, 57, f"Reports not retrieved (n={NOT_RETRIEVED})")
    box(ax, main_x0, 35.5, main_x1, 43, f"Reports assessed for\neligibility (n={ASSESSED})")
    reasons = [f"{name} (n={n})" for name, n in EXCLUSION_REASONS]
    box(ax, excl_x0, 12, excl_x1, 46.5, f"Reports excluded (n={FULLTEXT_EXCLUDED})", align="list", lines=reasons)
    assert sum(n for _, n in EXCLUSION_REASONS) == FULLTEXT_EXCLUDED
    assert PUBMED + SCOPUS - DUPLICATES == SCREENED and SCREENED - SCREEN_EXCLUDED == SOUGHT
    assert ASSESSED - FULLTEXT_EXCLUDED == DATABASE_INCLUDED

    arrow(ax, main_x1, 82, excl_x0, 82)
    arrow(ax, mid, 76, mid, 68.5)
    arrow(ax, main_x1, 65.25, excl_x0, 65.25)
    arrow(ax, mid, 62, mid, 57)
    arrow(ax, main_x1, 53.25, excl_x0, 53.25)
    arrow(ax, mid, 49.5, mid, 43)
    arrow(ax, main_x1, 39.25, excl_x0, 39.25)

    # Other methods: backward and forward citation searching
    box(ax, oth_x0, 72.5, oth_x1, 88.5,
        f"Reports identified by backward\nand forward citation searching\nof included studies (n={CITATION_IDENTIFIED})\n\n")
    ax.text(oth_mid, 75.6, f"of which also among the database records\nbut not selected there (n={CITATION_ALSO_IN_DATABASE})",
            ha="center", va="center", color=TEXT, fontsize=10.5, fontweight=WEIGHT, linespacing=1.3)
    box(ax, oth_x0, 35.5, oth_x1, 43, f"Reports assessed for\neligibility (n={CITATION_ASSESSED})")
    arrow(ax, oth_mid, 72.5, oth_mid, 43)

    # Inclusion
    total = DATABASE_INCLUDED + CITATION_INCLUDED
    box(ax, main_x0, 1.5, main_x1, 15.5,
        f"Studies included in review (n={total})\n"
        f"via databases (n={DATABASE_INCLUDED})\nvia citation searching (n={CITATION_INCLUDED})")
    arrow(ax, mid, 35.5, mid, 15.5)
    ax.plot([oth_mid, oth_mid], [35.5, 8.5], color=EDGE, linewidth=1.6)
    arrow(ax, oth_mid, 8.5, main_x1, 8.5)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for ext, kw in (("png", {"dpi": 300}), ("pdf", {})):
        path = OUTPUT_DIR / f"prisma_flow.{ext}"
        fig.savefig(path, bbox_inches="tight", facecolor="white", **kw)
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
