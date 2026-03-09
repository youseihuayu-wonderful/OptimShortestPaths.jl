#!/usr/bin/env python3
"""
Generate publication-quality figures for OptimShortestPaths.jl paper.

All benchmark data comes from actual measured runs.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 12,
    "axes.linewidth": 1.2,
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#333333",
    "xtick.color": "#333333",
    "ytick.color": "#333333",
    "text.color": "#333333",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.grid": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 1.0,
    "ytick.major.width": 1.0,
})


# ===== Figure 1: DMY vs Dijkstra Runtime Scaling ===========================
def figure1():
    nodes = [100, 500, 1000, 5000, 10000, 50000, 100000]
    networkx_ms = [0.2, 1.1, 2.1, 15.4, 35.2, 334.1, 785.2]
    dijkstra_julia_ms = [0.02, 0.12, 0.25, 1.59, 3.22, 32.38, 73.73]
    dmy_julia_ms = [0.09, 0.55, 1.03, 6.19, 14.73, 127.69, 297.90]

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.loglog(nodes, networkx_ms, "s-", color="#d62728", markersize=7,
              linewidth=2, label="NetworkX Dijkstra (Python)")
    ax.loglog(nodes, dijkstra_julia_ms, "o-", color="#1f77b4", markersize=7,
              linewidth=2, label="Dijkstra binary heap (Julia)")
    ax.loglog(nodes, dmy_julia_ms, "^-", color="#2ca02c", markersize=7,
              linewidth=2, label="DMY (Julia)")

    ax.set_xlabel("Graph size (nodes)", fontsize=14)
    ax.set_ylabel("Runtime (ms)", fontsize=14)
    ax.set_title("Single-Source Shortest Path: Runtime Scaling",
                 fontsize=16, fontweight="bold", pad=12)
    ax.legend(fontsize=11, frameon=True, edgecolor="#cccccc",
              fancybox=False, loc="upper left")

    # Annotations
    ax.annotate("Julia compiled: 10\u201312\u00d7 faster than Python",
                xy=(50000, 334.1), xytext=(2000, 600),
                fontsize=10.5, fontstyle="italic", color="#d62728",
                arrowprops=dict(arrowstyle="->", color="#d62728", lw=1.2))

    ax.annotate("DMY: 4\u00d7 overhead vs Dijkstra\n(constant factor)",
                xy=(50000, 127.69), xytext=(300, 80),
                fontsize=10.5, fontstyle="italic", color="#2ca02c",
                arrowprops=dict(arrowstyle="->", color="#2ca02c", lw=1.2))

    ax.tick_params(which="both", length=5)
    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, "dmy_vs_dijkstra_scaling.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ===== Figure 2: Validation Across Datasets =================================
def figure2():
    datasets = ["Curated\n(18 drugs)", "ChEMBL\n(1553 drugs)", "Hetionet\n(47K nodes)"]
    auroc = [0.806, 0.909, 0.794]
    mrr = [0.691, 0.733, 0.210]

    x = np.arange(len(datasets))
    width = 0.32

    fig, ax = plt.subplots(figsize=(10, 6))

    bars1 = ax.bar(x - width / 2, auroc, width, color="#1f77b4", edgecolor="white",
                   linewidth=0.8, label="AUROC", zorder=3)
    bars2 = ax.bar(x + width / 2, mrr, width, color="#ff7f0e", edgecolor="white",
                   linewidth=0.8, label="MRR", zorder=3)

    # Value labels
    for bar in bars1:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=11,
                fontweight="bold", color="#1f77b4")
    for bar in bars2:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                f"{bar.get_height():.3f}", ha="center", va="bottom", fontsize=11,
                fontweight="bold", color="#ff7f0e")

    # Random baseline
    ax.axhline(y=0.5, color="#999999", linestyle="--", linewidth=1.2, zorder=2)
    ax.text(len(datasets) - 0.55, 0.51, "random baseline", fontsize=10,
            color="#999999", fontstyle="italic")

    ax.set_ylim(0, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(datasets, fontsize=12)
    ax.set_ylabel("Score", fontsize=14)
    ax.set_title("Prediction Quality Across Biological Datasets",
                 fontsize=16, fontweight="bold", pad=12)
    ax.legend(fontsize=12, frameon=True, edgecolor="#cccccc", fancybox=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, "dataset_validation.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ===== Figure 3: Drug Repurposing Case Study ================================
def figure3():
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.set_xlim(-0.5, 12.5)
    ax.set_ylim(-1.5, 3.0)
    ax.axis("off")

    # Node definitions: (label, center_x, color)
    nodes = [
        ("Sorafenib", 1.5, "#4393c3"),       # Drug – blue
        ("VEGFR2", 4.5, "#5aae61"),           # Target – green
        ("VEGF\nSignaling", 7.5, "#f4a582"),  # Pathway – orange
        ("Hepatocellular\nCarcinoma", 10.5, "#d6604d"),  # Disease – red
    ]
    box_w, box_h = 2.2, 1.1
    edge_weights = ["w = 0.22", "w = 0.69", "w = 0.22"]

    # Draw boxes
    for label, cx, color in nodes:
        rect = mpatches.FancyBboxPatch(
            (cx - box_w / 2, 1.0 - box_h / 2), box_w, box_h,
            boxstyle="round,pad=0.15", facecolor=color, edgecolor="#333333",
            linewidth=1.5, alpha=0.90, zorder=3)
        ax.add_patch(rect)
        ax.text(cx, 1.0, label, ha="center", va="center",
                fontsize=12, fontweight="bold", color="white", zorder=4)

    # Draw arrows and weights
    for i in range(len(nodes) - 1):
        x_start = nodes[i][1] + box_w / 2
        x_end = nodes[i + 1][1] - box_w / 2
        ax.annotate("", xy=(x_end, 1.0), xytext=(x_start, 1.0),
                     arrowprops=dict(arrowstyle="-|>", color="#333333",
                                     lw=2, mutation_scale=18), zorder=2)
        mid_x = (x_start + x_end) / 2
        ax.text(mid_x, 0.20, edge_weights[i], ha="center", va="center",
                fontsize=11, color="#555555", fontstyle="italic")

    # Legend for node types
    type_labels = [("Drug", "#4393c3"), ("Target", "#5aae61"),
                   ("Pathway", "#f4a582"), ("Disease", "#d6604d")]
    legend_patches = [mpatches.Patch(facecolor=c, edgecolor="#333333",
                                     label=t, linewidth=0.8) for t, c in type_labels]
    ax.legend(handles=legend_patches, loc="upper right", fontsize=10,
              frameon=True, edgecolor="#cccccc", fancybox=False, ncol=4)

    # Bottom summary
    ax.text(6.0, -1.0,
            "Total path cost: 1.13 (DMY) vs 1.13 (Dijkstra) \u2014 verified identical",
            ha="center", va="center", fontsize=12, color="#333333",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0",
                      edgecolor="#cccccc", linewidth=1))

    ax.set_title("Case Study: Sorafenib \u2192 Hepatocellular Carcinoma Repurposing Path",
                 fontsize=15, fontweight="bold", pad=10)

    fig.tight_layout()
    path = os.path.join(FIGURES_DIR, "case_study_path.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {path}")


# ===== Main =================================================================
if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    print("\nAll figures generated successfully.")
