#!/usr/bin/env python
"""Validation figures for the README (docs/figures/). Needs the Boileau et al. 2015 data and matplotlib.

    python scripts/make_figures.py --data /path/to/SuppMaterial_1DBenchmark_Boileau_etal

Colours: arterial1d blue, published FVM orange (dashed), reference (3-D / in vitro / theory /
published FVM for ADAN56) grey. Validated as a categorical pair (light mode).
"""
import argparse
from pathlib import Path

import numpy as np

from arterial1d import MMHG, benchmarks

BLUE, ORANGE, GREY, GRID, INK = "#2a78d6", "#eb6834", "#8a8984", "#e6e5e0", "#52514e"
OUT = Path(__file__).resolve().parent.parent / "docs" / "figures"


def style(ax):
    ax.grid(color=GRID, lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK)
    ax.tick_params(colors=INK, labelsize=8)


def ref_line(ax, x, y, label=None):
    ax.plot(x, y, color=GREY, lw=3.2, alpha=0.55, label=label, solid_capstyle="round")


def ours(ax, x, y, label=None):
    ax.plot(x, y, color=BLUE, lw=1.6, label=label)


def fvm(ax, x, y, label=None):
    ax.plot(x, y, color=ORANGE, lw=1.3, ls=(0, (3, 3)), label=label, zorder=3)   # on top: it often coincides


def header(fig, ax, title):
    """Title and a one-row legend at fixed distances (inches) from the top, then tight layout below."""
    h = fig.get_figheight()
    fig.suptitle(title, fontsize=10, color=INK, x=0.01, y=1 - 0.12 / h, ha="left", va="top")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.01, 1 - 0.42 / h), ncol=len(labels),
               frameon=False, fontsize=8, labelcolor=INK)
    fig.tight_layout(rect=(0, 0, 1, 1 - 0.72 / h))


def case1(plt, root):
    o = benchmarks.single_pulse(root, traces=True)
    fig, axes = plt.subplots(2, 1, figsize=(8, 4.6), sharex=True)
    for ax, kind in zip(axes, ("inviscid", "viscous")):
        tr = o[kind]["trace"]
        for i, t in enumerate(benchmarks.PULSE_TIMES):
            ours(ax, tr["x"], tr["P"][t], "arterial1d" if i == 0 else None)
            fvm(ax, tr["x_FVM"], tr["P_FVM"][t], "published FVM" if i == 0 else None)
        ax.set_ylabel("pressure (Pa)", fontsize=9, color=INK)
        ax.set_title(f"{kind}: pressure along the tube at t = 0.1, 0.3, ..., 1.5 s", fontsize=9, color=INK, loc="left")
        style(ax)
    xs = np.linspace(0.3, 9.0, 100)
    p0 = o["inviscid"]["peak_Pa"][0]
    ref_line(axes[1], xs, p0 * np.exp(-o["eq15"]["rate_per_m"] * xs), "Eq. 15 peak decay")
    axes[1].set_xlabel("distance from inlet (m)", fontsize=9, color=INK)
    header(fig, axes[1], f"Case 1, single pulse (dx = {o['dx'] * 1e3:.1f} mm, Courant {o['courant']:.2f})")
    fig.savefig(OUT / "case1_single_pulse.png", dpi=150)


def cases_2_4(plt, root):
    rows = [("carotid", "mid", "case 2 carotid, mid"), ("thoracic_aorta", "mid", "case 3 thoracic aorta, mid"),
            ("bifurcation", "aorta_mid", "case 4 aorta, mid"), ("bifurcation", "aorta_end", "case 4 aorta, end"),
            ("bifurcation", "iliac_mid", "case 4 iliac, mid")]
    runs = {c: benchmarks.run_case(c, root, traces=True) for c in ("carotid", "thoracic_aorta", "bifurcation")}
    fig, axes = plt.subplots(len(rows), 2, figsize=(8, 9.5), sharex=False)
    for i, (case, site, label) in enumerate(rows):
        tr = runs[case]["sites"][site]["trace"]
        e = runs[case]["sites"][site]["vs_3D"]
        for j, (q, scale, unit) in enumerate((("P", 1 / MMHG, "mmHg"), ("Q", 1e6, "mL/s"))):
            ax = axes[i, j]
            ref_line(ax, tr["t"], tr["ref"][q] * scale, "3-D (Xiao et al. 2014)" if i == 0 else None)
            ours(ax, tr["t"], tr["mine"][q] * scale, "arterial1d" if i == 0 else None)
            fvm(ax, tr["t"], tr["FVM"][q] * scale, "published FVM" if i == 0 else None)
            ax.set_title(f"{label}: {'pressure' if q == 'P' else 'flow'} (RMS error vs 3-D {e[q]['RMS']:.2f}%)",
                         fontsize=8, color=INK, loc="left")
            ax.set_ylabel(unit, fontsize=8, color=INK)
            style(ax)
    for ax in axes[-1]:
        ax.set_xlabel("time in cycle (s)", fontsize=8, color=INK)
    header(fig, axes[0, 0], "Cases 2-4, against 3-D fluid-structure solutions")
    fig.savefig(OUT / "cases2-4_vs_3d.png", dpi=150)


def case5(plt, root):
    o = benchmarks.network37(root, traces=True)
    sites = list(o["sites"])
    fig, axes = plt.subplots(len(sites), 2, figsize=(8, 13))
    for i, s in enumerate(sites):
        tr, e = o["sites"][s]["trace"], o["sites"][s]["vs_exp"]
        for j, (q, scale, unit) in enumerate((("P", 1 / MMHG, "mmHg"), ("Q", 1e6, "mL/s"))):
            ax = axes[i, j]
            ref_line(ax, tr["t"], tr["ref"][q] * scale, "in vitro (Matthys et al. 2007)" if i == 0 else None)
            ours(ax, tr["t"], tr["mine"][q] * scale, "arterial1d" if i == 0 else None)
            fvm(ax, tr["t"], tr["FVM"][q] * scale, "published FVM" if i == 0 else None)
            ax.set_title(f"{s}: {'pressure' if q == 'P' else 'flow'} (RMS error vs in vitro {e[q]['RMS']:.1f}%)",
                         fontsize=8, color=INK, loc="left")
            ax.set_ylabel(unit, fontsize=8, color=INK)
            style(ax)
    for ax in axes[-1]:
        ax.set_xlabel("time in cycle (s)", fontsize=8, color=INK)
    header(fig, axes[0, 0], "Case 5, 37-artery network, against in-vitro measurements")
    fig.savefig(OUT / "case5_in_vitro.png", dpi=150)


def adan56(plt, root):
    o = benchmarks.adan56(root, traces=True)
    fig, axes = plt.subplots(4, 3, figsize=(8, 8.5), sharex=True)
    for ax, (site, v) in zip(axes.flat, o.items()):
        tr = v["trace"]
        ref_line(ax, tr["t"], tr["P_ref"] / MMHG, "published FVM (reference)")
        ours(ax, tr["t"], tr["P"] / MMHG, "arterial1d")
        ax.set_title(f"{site} ({100 * v['P_rms_over_PP']:.2f}% of PP)", fontsize=8, color=INK, loc="left")
        style(ax)
    for ax in axes[:, 0]:
        ax.set_ylabel("mmHg", fontsize=8, color=INK)
    for ax in axes[-1]:
        ax.set_xlabel("time in cycle (s)", fontsize=8, color=INK)
    header(fig, axes[0, 0], "Case 6, ADAN56: mid-vessel pressure against the published FVM solution")
    fig.savefig(OUT / "case6_adan56.png", dpi=150)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", help=f"supplementary data folder (default: ${benchmarks.ENV})")
    a = ap.parse_args()
    root = benchmarks.reference_root(a.data)
    if root is None:
        raise SystemExit(f"Boileau et al. 2015 data not found; pass --data or set {benchmarks.ENV}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    for f in (case1, cases_2_4, case5, adan56):
        f(plt, root)
        print("done:", f.__name__)


if __name__ == "__main__":
    main()
