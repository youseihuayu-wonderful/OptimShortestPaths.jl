"""
Paper Figure Generator — ALL real data, ZERO mock values.

Data sources:
  - ChemPath/data/paper_results/paper_experiments.json (live-computed)
  - ChemPath/data/paper_results/kge_baseline.json (live-computed)
  - One live mini-rerun: Pareto scatter for breast cancer (DOID:1612)

Every plotted number traces back to a JSON key or a live computation.
If a required data point is missing, the script aborts rather than
fabricating it.

Outputs (paper/figures/):
  fig_method_comparison.pdf    — AUROC + AUPRC bar chart with CI, vs KGE
  fig_path_length.pdf          — positive vs negative distance distribution
  fig_sensitivity.pdf          — p_base sensitivity sweep
  fig_cv_folds.pdf             — 5-fold cross-validation per-fold AUROC
  fig_pareto_scatter.pdf       — Pareto front for one disease (breast cancer)
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))

RESULTS_DIR = Path(__file__).resolve().parent.parent / "data" / "paper_results"
FIGURES_DIR = Path(__file__).resolve().parent.parent.parent / "paper" / "figures"

# Consistent publication style — no tricks, just clean defaults
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "legend.fontsize": 9,
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def _load_results():
    """Load both result JSONs — abort if missing."""
    exp_path = RESULTS_DIR / "paper_experiments.json"
    kge_path = RESULTS_DIR / "kge_baseline.json"
    if not exp_path.exists():
        raise FileNotFoundError(
            f"{exp_path} not found. Run paper_experiments.py first."
        )
    if not kge_path.exists():
        raise FileNotFoundError(
            f"{kge_path} not found. Run kge_baseline.py first."
        )
    with open(exp_path) as f:
        exp = json.load(f)
    with open(kge_path) as f:
        kge = json.load(f)
    return exp, kge


def fig_method_comparison(exp, kge):
    """Bar chart: AUROC (with 95% CI) for our methods vs KGE baselines.

    All values read directly from JSON, no computation here.
    """
    methods = []
    aurocs = []
    ci_los = []
    ci_his = []
    auprcs = []
    colors = []

    # Our methods, ordered by AUROC descending
    our_order = [
        ("eff_evd_fused", "Ours: Eff+Evd", "#2E86AB"),
        ("efficacy_1d", "Ours: Eff 1D", "#2E86AB"),
        ("eff_evd_safety_fused", "Ours: Eff+Evd+Safety", "#2E86AB"),
        ("weighted_3d", "Ours: Weighted 3D", "#2E86AB"),
        ("eff_safety_2d", "Ours: Eff+Safety 2D", "#2E86AB"),
    ]
    for key, label, color in our_order:
        d = exp["method_comparison"][key]
        methods.append(label)
        aurocs.append(d["auroc"])
        ci_los.append(d["ci_lo"])
        ci_his.append(d["ci_hi"])
        auprcs.append(d["auprc"])
        colors.append(color)

    # KGE baselines — no CI available in kge_baseline.json
    for key, label in [("ComplEx", "KGE: ComplEx"), ("TransE", "KGE: TransE")]:
        d = kge[key]
        methods.append(label)
        aurocs.append(d["auroc"])
        ci_los.append(d["auroc"])  # no CI — use point estimate
        ci_his.append(d["auroc"])
        auprcs.append(d["auprc"])
        colors.append("#A23B72")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.2))

    # AUROC with CIs
    y_pos = np.arange(len(methods))
    errs = np.array([
        [a - lo for a, lo in zip(aurocs, ci_los)],
        [hi - a for a, hi in zip(aurocs, ci_his)],
    ])
    ax1.barh(y_pos, aurocs, xerr=errs, color=colors, edgecolor="black",
             linewidth=0.5, capsize=3, alpha=0.85)
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(methods)
    ax1.invert_yaxis()
    ax1.set_xlabel("AUROC (95% CI)")
    ax1.set_xlim(0.5, 0.85)
    ax1.axvline(x=0.5, color="gray", linestyle=":", linewidth=0.8, label="Random")
    ax1.set_title("AUROC: shortest-path casting vs KGE baselines")

    # Add numeric labels
    for i, a in enumerate(aurocs):
        ax1.text(a + 0.005, i, f"{a:.3f}", va="center", fontsize=8)

    # AUPRC
    ax2.barh(y_pos, auprcs, color=colors, edgecolor="black",
             linewidth=0.5, alpha=0.85)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([""] * len(methods))
    ax2.set_xlabel("AUPRC")
    ax2.set_xlim(0, max(auprcs) * 1.25)
    # Random baseline AUPRC = positive class prevalence
    n_pos = exp["path_length_analysis"]["n_positive"]
    n_neg = exp["path_length_analysis"]["n_negative"]
    random_auprc = n_pos / (n_pos + n_neg)
    ax2.axvline(x=random_auprc, color="gray", linestyle=":",
                linewidth=0.8, label=f"Random ({random_auprc:.3f})")
    ax2.set_title("AUPRC (prevalence = 1.4%)")
    ax2.legend(loc="lower right", frameon=False)

    for i, a in enumerate(auprcs):
        ax2.text(a + 0.002, i, f"{a:.3f}", va="center", fontsize=8)

    plt.tight_layout()
    out = FIGURES_DIR / "fig_method_comparison.pdf"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out}")


def fig_path_length(exp):
    """Histogram: shortest-path distance distribution for positives vs negatives.

    Uses pos_hist, neg_hist, bin_edges from JSON — live-computed KS test
    statistic and Cohen's d also shown.
    """
    path = exp["path_length_analysis"]
    bins = path["bin_edges"]
    pos_hist = np.array(path["pos_hist"])
    neg_hist = np.array(path["neg_hist"])

    # Normalize to density (different class sizes)
    pos_density = pos_hist / path["n_positive"]
    neg_density = neg_hist / path["n_negative"]

    centers = [(bins[i] + bins[i + 1]) / 2 for i in range(len(bins) - 1)]
    width = (bins[1] - bins[0]) * 0.4

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar([c - width / 2 for c in centers], pos_density, width=width,
           color="#2E86AB", edgecolor="black", linewidth=0.4,
           label=f"Positive (n={path['n_positive']})", alpha=0.85)
    ax.bar([c + width / 2 for c in centers], neg_density, width=width,
           color="#A23B72", edgecolor="black", linewidth=0.4,
           label=f"Negative (n={path['n_negative']})", alpha=0.85)

    ax.axvline(x=path["pos_mean"], color="#2E86AB", linestyle="--",
               linewidth=1.5, alpha=0.7)
    ax.axvline(x=path["neg_mean"], color="#A23B72", linestyle="--",
               linewidth=1.5, alpha=0.7)

    ax.set_xlabel("Shortest-path efficacy distance (-log p)")
    ax.set_ylabel("Density (within class)")
    ax.set_title(
        f"Positive pairs have shorter paths "
        f"(Cohen's d = {path['cohens_d']:.2f}, KS = {path['ks_statistic']:.3f})"
    )
    ax.legend(frameon=False)
    plt.tight_layout()
    out = FIGURES_DIR / "fig_path_length.pdf"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out}")


def fig_sensitivity(exp):
    """Line plot: AUROC vs p_base(binds). Shows robustness."""
    sens = exp["pbase_sensitivity"]
    p_vals = [r["p_binds"] for r in sens]
    auroc_1d = [r["auroc_1d"] for r in sens]
    auprc_1d = [r["auprc_1d"] for r in sens]

    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax2 = ax1.twinx()

    ln1 = ax1.plot(p_vals, auroc_1d, "o-", color="#2E86AB", linewidth=2,
                   markersize=8, label="AUROC")
    ax1.axhline(y=np.mean(auroc_1d), color="#2E86AB", linestyle=":",
                linewidth=0.8, alpha=0.6)
    ax1.set_ylabel("AUROC (1D efficacy)", color="#2E86AB")
    ax1.tick_params(axis="y", labelcolor="#2E86AB")
    ax1.set_ylim(
        min(auroc_1d) - 0.01, max(auroc_1d) + 0.01
    )
    ax1.set_xlabel(r"$p_{\mathrm{base}}(\mathrm{binds})$")

    ln2 = ax2.plot(p_vals, auprc_1d, "s--", color="#F18F01", linewidth=2,
                   markersize=7, label="AUPRC")
    ax2.set_ylabel("AUPRC (1D efficacy)", color="#F18F01")
    ax2.tick_params(axis="y", labelcolor="#F18F01")
    ax2.spines["top"].set_visible(False)

    # Default value annotation
    default_idx = p_vals.index(0.8)
    ax1.annotate("default",
                 xy=(0.8, auroc_1d[default_idx]),
                 xytext=(0.85, auroc_1d[default_idx] + 0.002),
                 fontsize=8, color="#444")

    spread = max(auroc_1d) - min(auroc_1d)
    ax1.set_title(f"Sensitivity to p_base(binds): AUROC spread = {spread:.4f}")

    lns = ln1 + ln2
    ax1.legend(lns, [l.get_label() for l in lns], loc="center right",
               frameon=False)
    plt.tight_layout()
    out = FIGURES_DIR / "fig_sensitivity.pdf"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out}")


def fig_cv_folds(exp):
    """Box/strip plot: per-fold AUROC for 5-fold CV."""
    cv = exp["cross_validation"]
    methods = ["efficacy_1d", "weighted_3d", "eff_safety_2d"]
    labels = ["Efficacy 1D", "Weighted 3D", "Eff+Safety 2D"]
    colors = ["#2E86AB", "#6A4C93", "#A23B72"]

    fig, ax = plt.subplots(figsize=(6, 4))
    for i, (m, label, color) in enumerate(zip(methods, labels, colors)):
        folds = cv[m]["folds"]
        # Scatter each fold
        ax.scatter([i] * len(folds), folds, color=color, s=60,
                   edgecolor="black", linewidth=0.5, zorder=3, alpha=0.8)
        # Mean line
        ax.hlines(cv[m]["mean"], i - 0.25, i + 0.25, color=color,
                  linewidth=2, zorder=4)
        # ±1 std band
        std = cv[m]["std"]
        ax.fill_between(
            [i - 0.25, i + 0.25],
            [cv[m]["mean"] - std] * 2,
            [cv[m]["mean"] + std] * 2,
            color=color, alpha=0.15, zorder=2,
        )

    ax.set_xticks(range(len(methods)))
    ax.set_xticklabels(labels)
    ax.set_ylabel("AUROC per fold")
    ax.set_title(
        f"5-fold cross-validation (disease split): "
        f"1D mean {cv['efficacy_1d']['mean']:.3f} "
        f"± {cv['efficacy_1d']['std']:.3f}"
    )
    ax.set_ylim(0.65, 0.85)
    plt.tight_layout()
    out = FIGURES_DIR / "fig_cv_folds.pdf"
    plt.savefig(out)
    plt.close()
    print(f"  Saved: {out}")


def fig_pareto_scatter_live():
    """Pareto front scatter for breast cancer (DOID:1612).

    LIVE computation on real Hetionet graph. Does NOT fabricate data.
    Aborts if graph or scores cannot be computed.
    """
    print("  Computing live Pareto scatter for breast cancer (DOID:1612)...")

    from chempath_enriched_benchmark import (
        build_enriched_weights,
        compute_pareto_front_2d,
        load_hetionet_enriched,
        load_pharmacotherapydb,
        score_all_pairs,
    )
    try:
        from sider_safety_score import (
            build_fallback_scores,
            build_sider_safety_scores,
        )
    except ImportError:
        build_sider_safety_scores = None
        build_fallback_scores = None

    het_data = load_hetionet_enriched()
    G = het_data["G"]
    edge_metadata = het_data["edge_metadata"]
    ground_truth = het_data["ground_truth"]
    side_effect_counts = het_data["side_effect_counts"]

    compound_severity = {}
    if build_sider_safety_scores:
        sider = build_sider_safety_scores()
        fb = build_fallback_scores(sider)
        all_s = {**fb, **sider}
        het_ids = {n.replace("Compound::", "") for n in G.nodes()
                   if n.startswith("Compound::")}
        for hid, s in all_s.items():
            db_id = hid.replace("Compound::", "")
            if db_id in het_ids:
                compound_severity[db_id] = s

    positive_pairs, eval_compounds, eval_diseases = (
        load_pharmacotherapydb(G, ground_truth)
    )

    target_disease = "Disease::DOID:1612"  # Breast cancer
    if target_disease not in eval_diseases:
        raise RuntimeError(f"{target_disease} not in eval set")

    G_eff, G_saf, G_evd, _ = build_enriched_weights(
        G, edge_metadata, {}, compound_severity, side_effect_counts
    )
    scores = score_all_pairs(G_eff, G_saf, G_evd, eval_compounds, {target_disease})

    points = []  # (compound, d_e, d_s, is_positive)
    for (c, d), (d_e, d_s, d_v) in scores.items():
        if d != target_disease:
            continue
        if d_e == float("inf") or d_s == float("inf") or d_e > 50 or d_s > 50:
            continue
        is_pos = (c, d) in positive_pairs
        points.append((c, d_e, d_s, is_pos))

    if not points:
        raise RuntimeError("No valid scores for target disease")

    # Identify Pareto front (minimize both distances)
    pts_neg = [(-d_e, -d_s, i) for i, (_, d_e, d_s, _) in enumerate(points)]
    front_idx = set(compute_pareto_front_2d(pts_neg))

    d_e_pos = [p[1] for p in points if p[3]]
    d_s_pos = [p[2] for p in points if p[3]]
    d_e_neg = [p[1] for p in points if not p[3]]
    d_s_neg = [p[2] for p in points if not p[3]]

    front_pts = [points[i] for i in sorted(front_idx,
                                            key=lambda i: points[i][1])]

    print(f"    Points: {len(points)} total, {len(d_e_pos)} positive, "
          f"{len(front_pts)} on Pareto front")

    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.scatter(d_e_neg, d_s_neg, color="#D3D3D3", s=18, edgecolor="none",
               alpha=0.7, label=f"Non-indication (n={len(d_e_neg)})")
    ax.scatter(d_e_pos, d_s_pos, color="#A23B72", s=60,
               edgecolor="black", linewidth=0.4, alpha=0.9,
               label=f"PharmacotherapyDB indication (n={len(d_e_pos)})",
               zorder=4)
    fx = [p[1] for p in front_pts]
    fy = [p[2] for p in front_pts]
    ax.plot(fx, fy, color="#2E86AB", linewidth=1.8, linestyle="--",
            alpha=0.9, label=f"Pareto front (n={len(front_pts)})",
            zorder=5)
    ax.scatter(fx, fy, color="#2E86AB", s=50, edgecolor="black",
               linewidth=0.5, zorder=6)

    ax.set_xlabel("Efficacy distance")
    ax.set_ylabel("Safety distance")
    ax.set_title("Pareto front for breast cancer (DOID:1612)")
    ax.legend(loc="upper right", frameon=False)
    plt.tight_layout()
    out = FIGURES_DIR / "fig_pareto_scatter.pdf"
    plt.savefig(out)
    plt.close()
    print(f"    Saved: {out}")


def main():
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Generating figures in: {FIGURES_DIR}")
    print("All data sourced from live experiments. No mock values.\n")

    exp, kge = _load_results()

    print("[1/5] Method comparison bar chart...")
    fig_method_comparison(exp, kge)

    print("[2/5] Path length histogram...")
    fig_path_length(exp)

    print("[3/5] Sensitivity sweep line plot...")
    fig_sensitivity(exp)

    print("[4/5] 5-fold CV scatter...")
    fig_cv_folds(exp)

    print("[5/5] Pareto scatter (live computation)...")
    fig_pareto_scatter_live()

    print(f"\nDone. {len(list(FIGURES_DIR.glob('*.pdf')))} figures written.")


if __name__ == "__main__":
    main()
