"""
ChemPath Drug Screening Dashboard

An interactive tool for exploring drug-target interactions using
graph-based shortest-path optimization on real ChEMBL bioactivity data.

Run with: uv run streamlit run chempath/ui/dashboard.py
"""

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import math
from pathlib import Path

from chempath.data.chembl_client import load_saved_data, ANTICANCER_TARGETS
from chempath.graph.builder import build_drug_target_graph, get_graph_summary
from chempath.graph.optimizer import (
    rank_compounds_for_target,
    compute_pareto_front,
    find_knee_point,
    compute_selectivity,
    compare_compounds_across_targets,
)
from chempath.chemistry.properties import compute_properties, compute_enriched_properties
from chempath.graph.validation import (
    validate_target,
    compute_multi_objective_scores,
    compute_multi_objective_pareto,
)

DATA_PATH = Path(__file__).parent.parent.parent / "data" / "chembl_real.json"

# ---------------------------------------------------------------------------
# Theme & constants
# ---------------------------------------------------------------------------
COLORS = {
    "primary": "#0066FF",
    "secondary": "#6C63FF",
    "success": "#00C48C",
    "warning": "#FFB020",
    "danger": "#FF4D4F",
    "bg_card": "#F8F9FC",
    "text_muted": "#8C8C8C",
    "high": "#00C48C",
    "medium": "#FFB020",
    "low": "#FF4D4F",
}

CONFIDENCE_COLORS = {"HIGH": COLORS["high"], "MEDIUM": COLORS["medium"], "LOW": COLORS["low"]}

PLOTLY_TEMPLATE = "plotly_white"

# Shared Plotly layout defaults for consistency
_BASE_LAYOUT = dict(
    template=PLOTLY_TEMPLATE,
    font=dict(family="Inter, system-ui, sans-serif", size=13),
    margin=dict(l=48, r=24, t=56, b=48),
    hoverlabel=dict(bgcolor="white", font_size=12, bordercolor="#E0E0E0"),
    paper_bgcolor="white",
    plot_bgcolor="white",
)


# ---------------------------------------------------------------------------
# Data loading (cached)
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="ChemPath | Drug Screening Dashboard",
    page_icon="https://img.icons8.com/fluency/48/test-tube.png",
    layout="wide",
    initial_sidebar_state="expanded",
)


@st.cache_data(show_spinner=False)
def load_data():
    return load_saved_data(DATA_PATH)


@st.cache_resource(show_spinner=False)
def build_graph(toxicity_penalty: float):
    data = load_data()
    return build_drug_target_graph(data, toxicity_penalty=toxicity_penalty, verbose=False)


def target_lookup():
    return {t["name"]: t["chembl_id"] for t in ANTICANCER_TARGETS}


def target_descriptions():
    return {t["name"]: t.get("description", t["name"]) for t in ANTICANCER_TARGETS}


# ---------------------------------------------------------------------------
# Custom CSS for professional look
# ---------------------------------------------------------------------------
st.markdown("""
<style>
    /* Remove default Streamlit padding */
    .block-container { padding-top: 1.5rem; }

    /* Sidebar styling */
    section[data-testid="stSidebar"] > div:first-child {
        padding-top: 1rem;
    }

    /* Card-like metric containers */
    div[data-testid="stMetric"] {
        background: #F8F9FC;
        border: 1px solid #E8EAF0;
        border-radius: 10px;
        padding: 12px 16px;
    }
    div[data-testid="stMetric"] label {
        color: #6B7280 !important;
        font-size: 0.78rem !important;
        text-transform: uppercase;
        letter-spacing: 0.04em;
    }
    div[data-testid="stMetric"] div[data-testid="stMetricValue"] {
        font-size: 1.5rem !important;
        font-weight: 700 !important;
        color: #1A1A2E !important;
    }

    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 0px;
        border-bottom: 2px solid #E8EAF0;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 8px 24px;
        font-weight: 500;
    }

    /* Expander styling */
    .streamlit-expanderHeader {
        font-weight: 600;
        font-size: 0.95rem;
    }

    /* Info boxes */
    .info-box {
        background: linear-gradient(135deg, #F0F4FF 0%, #E8F0FE 100%);
        border-left: 4px solid #0066FF;
        border-radius: 0 8px 8px 0;
        padding: 16px 20px;
        margin-bottom: 16px;
        font-size: 0.9rem;
        line-height: 1.5;
        color: #1A1A2E;
    }

    /* Legend items */
    .legend-row {
        display: flex;
        gap: 24px;
        margin: 8px 0 16px 0;
        flex-wrap: wrap;
    }
    .legend-item {
        display: flex;
        align-items: center;
        gap: 6px;
        font-size: 0.82rem;
        color: #555;
    }
    .legend-dot {
        width: 10px;
        height: 10px;
        border-radius: 50%;
        display: inline-block;
    }
</style>
""", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("## ChemPath")
    st.caption("Drug Screening Dashboard  |  v0.1")
    st.markdown("---")

    # Target selection with descriptions
    target_opts = target_lookup()
    target_descs = target_descriptions()
    selected_target_name = st.selectbox(
        "Target Protein",
        list(target_opts.keys()),
        help="Select the protein target you want to screen compounds against.",
    )
    selected_target_id = target_opts[selected_target_name]
    st.caption(f"_{target_descs.get(selected_target_name, '')}_")

    st.markdown("---")
    st.markdown("##### Optimization Parameters")

    strategy = st.radio(
        "Ranking Strategy",
        ["balanced", "efficacy", "safety"],
        horizontal=True,
        help=(
            "**Balanced** (default): Optimizes both potency and safety via weighted shortest-path.  \n"
            "**Efficacy**: Ranks purely by IC50 (binding strength).  \n"
            "**Safety**: Ranks by toxicity first, then IC50."
        ),
    )
    strategy_labels = {
        "balanced": "Balanced (IC50 + Toxicity)",
        "efficacy": "Efficacy-first (IC50 only)",
        "safety": "Safety-first (Toxicity priority)",
    }

    top_n = st.slider("Show Top N Compounds", 5, 50, 15, help="Number of top-ranked compounds to display in tables and charts.")

    tox_penalty = st.slider(
        "Toxicity Penalty",
        0.0, 3.0, 1.0, 0.1,
        help=(
            "Controls how heavily toxicity risk penalizes graph edge weights.  \n"
            "**0.0** = ignore toxicity entirely. **1.0** = default. **3.0** = strongly penalize toxic compounds."
        ),
    )

    st.markdown("---")
    # Data summary
    data = load_data()
    data_compounds_count = len(data["compounds"])
    data_targets_count = len(data["targets"])
    data_activities_count = len(data["bioactivities"])
    phase_compounds = sum(1 for c in data["compounds"] if c.get("phase", 0) > 0)

    st.markdown("##### Data Summary")
    st.markdown(
        f"**{data_compounds_count:,}** compounds  &middot;  "
        f"**{data_targets_count}** targets  &middot;  "
        f"**{data_activities_count:,}** interactions"
    )
    st.markdown(f"**{phase_compounds}** with clinical phase data")
    st.caption("Source: ChEMBL REST API (experimental IC50)")

# ---------------------------------------------------------------------------
# Build graph
# ---------------------------------------------------------------------------
with st.spinner("Building optimization graph..."):
    G = build_graph(tox_penalty)
    summary = get_graph_summary(G)

# Compute rankings once — reused across tabs
all_recs = rank_compounds_for_target(G, selected_target_id, strategy=strategy)
total_compounds_for_target = len(all_recs)
display_recs = all_recs[:top_n]

# ---------------------------------------------------------------------------
# Hero header
# ---------------------------------------------------------------------------
st.markdown(f"# Drug Screening: **{selected_target_name}**")
st.markdown(
    f'<div class="info-box">'
    f"<strong>What is this?</strong> &nbsp; ChemPath uses <strong>shortest-path graph optimization</strong> "
    f"to rank {total_compounds_for_target:,} compounds tested against {selected_target_name} "
    f"({target_descs.get(selected_target_name, '')}).  "
    f"Each compound-target interaction is an edge weighted by binding affinity (IC50) and toxicity risk. "
    f"The <strong>{strategy_labels[strategy]}</strong> strategy is active."
    f"</div>",
    unsafe_allow_html=True,
)

# KPI row
kpi1, kpi2, kpi3, kpi4, kpi5 = st.columns(5)
kpi1.metric("Total Compounds", f"{total_compounds_for_target:,}")

high_count = sum(1 for r in all_recs if r.confidence == "HIGH")
kpi2.metric("HIGH Confidence", high_count)

sub_nm = sum(1 for r in all_recs if r.ic50_nm < 1)
kpi3.metric("Sub-nM Potency", sub_nm)

lipinski_pass = sum(1 for r in all_recs if r.toxicity == 0)
kpi4.metric("Lipinski-Clean", lipinski_pass, help="Compounds with zero Lipinski violations (MW<500, LogP<5, HBD<5, HBA<10)")

if all_recs:
    best = all_recs[0]
    kpi5.metric("Top Compound", best.compound_name)


# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab_ranking, tab_pareto, tab_network, tab_selectivity, tab_properties, tab_validation, tab_methods = st.tabs([
    "Rankings",
    "Pareto Optimization",
    "Interaction Network",
    "Selectivity Profile",
    "Drug-likeness",
    "Validation",
    "Methodology & Benchmarks",
])


# ===========================================================================
# TAB 1 — Rankings
# ===========================================================================
with tab_ranking:
    if not all_recs:
        st.warning(f"No compounds found for {selected_target_name}.")
    else:
        # --- Explanation ---
        with st.expander("How to read this table", expanded=False):
            st.markdown("""
**Rank** — Position in the optimized ordering (1 = best candidate).

**IC50 (nM)** — Half-maximal inhibitory concentration. Lower = more potent binding.
- < 1 nM: Exceptionally potent
- 1-100 nM: Potent
- 100-1000 nM: Moderate
- > 1000 nM: Weak

**Toxicity Risk** — Estimated from molecular properties (Lipinski violations, MW, LogP).
0.00 = no structural risk flags.

**Weight** — Graph edge weight used for optimization. Lower = better candidate.
Combines IC50 potency and toxicity penalty.

**Confidence** — Composite score based on:
- Data source (experimental vs predicted)
- Clinical phase (Phase 4 = approved drug)
- Potency tier
- Toxicity level
""")

        # --- Color-coded table ---
        st.markdown(f"### Top {len(display_recs)} Compounds")
        st.markdown(
            f'<div class="legend-row">'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["high"]}"></span> HIGH confidence</span>'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["medium"]}"></span> MEDIUM confidence</span>'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["low"]}"></span> LOW confidence</span>'
            f"</div>",
            unsafe_allow_html=True,
        )

        # Compute enriched properties for displayed compounds
        _enriched_cache = {}
        for r in display_recs:
            smiles = G.nodes[r.compound_id].get("smiles", "")
            if smiles:
                _enriched_cache[r.compound_id] = compute_enriched_properties(smiles)

        df = pd.DataFrame([
            {
                "Rank": r.rank,
                "Compound": r.compound_name,
                "IC50 (nM)": round(r.ic50_nm, 3),
                "Tox Risk": round(r.toxicity, 3),
                "QED": _enriched_cache[r.compound_id].qed_score if r.compound_id in _enriched_cache else 0,
                "LogP": _enriched_cache[r.compound_id].estimated_logp if r.compound_id in _enriched_cache else 0,
                "MW": _enriched_cache[r.compound_id].estimated_mw if r.compound_id in _enriched_cache else 0,
                "SA": _enriched_cache[r.compound_id].sa_score if r.compound_id in _enriched_cache else 0,
                "TPSA": _enriched_cache[r.compound_id].estimated_tpsa if r.compound_id in _enriched_cache else 0,
                "Weight": round(r.weight, 4),
                "Confidence": r.confidence,
            }
            for r in display_recs
        ])

        def highlight_confidence(val):
            colors = {"HIGH": f"background-color: {COLORS['high']}22; color: {COLORS['high']}; font-weight:600",
                       "MEDIUM": f"background-color: {COLORS['medium']}22; color: {COLORS['medium']}; font-weight:600",
                       "LOW": f"background-color: {COLORS['low']}22; color: {COLORS['low']}; font-weight:600"}
            return colors.get(val, "")

        styled = df.style.map(highlight_confidence, subset=["Confidence"]).format({
            "IC50 (nM)": lambda x: f"{x:.3f}" if x < 1 else (f"{x:.1f}" if x < 100 else f"{x:.0f}"),
            "Tox Risk": "{:.3f}",
            "QED": "{:.3f}",
            "LogP": "{:.2f}",
            "MW": "{:.0f}",
            "SA": "{:.1f}",
            "TPSA": "{:.0f}",
            "Weight": "{:.4f}",
        })
        st.dataframe(styled, hide_index=True, width="stretch")

        # --- Side-by-side charts ---
        chart_left, chart_right = st.columns(2)

        with chart_left:
            st.markdown("#### IC50 Distribution")
            log_ic50 = [math.log10(r.ic50_nm) if r.ic50_nm > 0 else 0 for r in all_recs]
            fig_hist = go.Figure()
            fig_hist.add_trace(go.Histogram(
                x=log_ic50,
                nbinsx=40,
                marker_color=COLORS["primary"],
                opacity=0.85,
                hovertemplate="log10(IC50)=%{x:.1f}<br>Count=%{y}<extra></extra>",
            ))
            # Add potency zone annotations
            fig_hist.add_vrect(x0=-2, x1=0, fillcolor=COLORS["high"], opacity=0.08,
                               annotation_text="Sub-nM", annotation_position="top left",
                               annotation_font_size=10, annotation_font_color=COLORS["high"])
            fig_hist.add_vrect(x0=0, x1=2, fillcolor=COLORS["medium"], opacity=0.05,
                               annotation_text="Potent", annotation_position="top left",
                               annotation_font_size=10, annotation_font_color=COLORS["medium"])
            fig_hist.update_layout(
                **_BASE_LAYOUT,
                xaxis_title="log10(IC50 nM)",
                yaxis_title="Number of Compounds",
                height=350,
                showlegend=False,
                bargap=0.05,
            )
            st.plotly_chart(fig_hist, width="stretch")

        with chart_right:
            st.markdown("#### Confidence Breakdown")
            conf_counts = {"HIGH": 0, "MEDIUM": 0, "LOW": 0}
            for r in all_recs:
                conf_counts[r.confidence] += 1

            fig_donut = go.Figure(data=[go.Pie(
                labels=list(conf_counts.keys()),
                values=list(conf_counts.values()),
                hole=0.55,
                marker_colors=[COLORS["high"], COLORS["medium"], COLORS["low"]],
                textinfo="label+percent",
                textfont_size=13,
                hovertemplate="%{label}: %{value} compounds (%{percent})<extra></extra>",
            )])
            fig_donut.update_layout(
                **_BASE_LAYOUT,
                height=350,
                showlegend=False,
                annotations=[dict(text=f"{total_compounds_for_target}", x=0.5, y=0.5,
                                   font_size=22, font_color=COLORS["primary"],
                                   showarrow=False)],
            )
            st.plotly_chart(fig_donut, width="stretch")

        # --- Scatter: IC50 vs Toxicity overview ---
        st.markdown("#### Potency vs Safety Overview")
        st.caption("Each dot is a compound. Hover to see details. The best candidates are in the bottom-left (low IC50, low toxicity).")
        scatter_df = pd.DataFrame([
            {"Compound": r.compound_name, "IC50 (nM)": r.ic50_nm,
             "Toxicity Risk": r.toxicity, "Confidence": r.confidence,
             "Weight": r.weight, "Rank": r.rank}
            for r in all_recs
        ])
        fig_scatter = px.scatter(
            scatter_df,
            x="IC50 (nM)", y="Toxicity Risk",
            color="Confidence",
            color_discrete_map=CONFIDENCE_COLORS,
            hover_name="Compound",
            hover_data={"Rank": True, "Weight": ":.4f", "IC50 (nM)": ":.2f", "Toxicity Risk": ":.3f"},
            log_x=True,
            category_orders={"Confidence": ["HIGH", "MEDIUM", "LOW"]},
        )
        fig_scatter.update_traces(marker=dict(size=7, line=dict(width=0.5, color="white")))
        fig_scatter.update_layout(
            **_BASE_LAYOUT,
            height=420,
            xaxis_title="IC50 (nM)  [log scale, lower = more potent]",
            yaxis_title="Toxicity Risk Score  [lower = safer]",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5, title=None),
        )
        # Ideal zone annotation
        fig_scatter.add_annotation(
            x=math.log10(0.5), y=0.0,
            text="Ideal Zone",
            showarrow=False,
            font=dict(size=11, color=COLORS["high"]),
            bgcolor="rgba(0,196,140,0.08)",
            bordercolor=COLORS["high"],
            borderwidth=1,
            borderpad=6,
        )
        st.plotly_chart(fig_scatter, width="stretch")


# ===========================================================================
# TAB 2 — Pareto Optimization
# ===========================================================================
with tab_pareto:
    mo_scores = compute_multi_objective_scores(G, selected_target_id)
    if not mo_scores:
        st.warning("No data for Pareto analysis.")
    else:
        mo_pareto = compute_multi_objective_pareto(mo_scores)

        # Also compute 2-objective Pareto for comparison
        recs_pareto = rank_compounds_for_target(G, selected_target_id, strategy="balanced")
        pareto_2obj = compute_pareto_front(recs_pareto)

        st.markdown("### Multi-Objective Pareto Optimization")
        st.markdown(
            f'<div class="info-box">'
            f"<strong>4-Objective Pareto front</strong> optimizes IC50 (potency), toxicity risk, "
            f"synthetic accessibility, and QED (drug-likeness) simultaneously. "
            f"Of {len(mo_scores):,} compounds, <strong>{len(mo_pareto)}</strong> are Pareto-optimal "
            f"across all 4 dimensions. "
            f"(For comparison, the old 2-objective front had only {len(pareto_2obj)} point(s).)"
            f"</div>",
            unsafe_allow_html=True,
        )

        # --- Main scatter: IC50 vs QED, colored by risk, sized by SA ---
        pareto_ids = {p["compound_id"] for p in mo_pareto}
        dominated = [s for s in mo_scores if s["compound_id"] not in pareto_ids]

        fig_pareto = go.Figure()

        # Dominated compounds
        if dominated:
            fig_pareto.add_trace(go.Scatter(
                x=[s["ic50_nm"] for s in dominated],
                y=[s["qed"] for s in dominated],
                mode="markers",
                marker=dict(
                    size=6,
                    color=[s["risk_score"] for s in dominated],
                    colorscale="RdYlGn_r",
                    cmin=0, cmax=0.6,
                    opacity=0.3,
                ),
                name="Dominated",
                text=[s["compound_name"] for s in dominated],
                customdata=[[s["risk_score"], s["sa_score"], s["phase"]] for s in dominated],
                hovertemplate=(
                    "%{text}<br>IC50: %{x:.2f} nM<br>QED: %{y:.3f}<br>"
                    "Risk: %{customdata[0]:.2f} | SA: %{customdata[1]:.1f} | Phase: %{customdata[2]}"
                    "<extra>Dominated</extra>"
                ),
            ))

        # Pareto-optimal compounds
        fig_pareto.add_trace(go.Scatter(
            x=[s["ic50_nm"] for s in mo_pareto],
            y=[s["qed"] for s in mo_pareto],
            mode="markers+text",
            marker=dict(
                size=[max(8, 20 - s["sa_score"] * 1.5) for s in mo_pareto],
                color=COLORS["primary"],
                symbol="diamond",
                line=dict(width=1.5, color="white"),
            ),
            name="Pareto-optimal",
            text=[s["compound_name"][:12] for s in mo_pareto],
            textposition="top center",
            textfont=dict(size=9),
            customdata=[[s["risk_score"], s["sa_score"], s["phase"], s["mw"]] for s in mo_pareto],
            hovertemplate=(
                "<b>%{text}</b><br>IC50: %{x:.2f} nM<br>QED: %{y:.3f}<br>"
                "Risk: %{customdata[0]:.2f} | SA: %{customdata[1]:.1f}<br>"
                "Phase: %{customdata[2]} | MW: %{customdata[3]:.0f}"
                "<extra>Pareto Front</extra>"
            ),
        ))

        # Highlight approved drugs on the Pareto front
        approved_pareto = [s for s in mo_pareto if s["phase"] >= 4]
        if approved_pareto:
            fig_pareto.add_trace(go.Scatter(
                x=[s["ic50_nm"] for s in approved_pareto],
                y=[s["qed"] for s in approved_pareto],
                mode="markers",
                marker=dict(size=16, color=COLORS["success"], symbol="star",
                            line=dict(width=2, color="white")),
                name="Approved Drug",
                text=[s["compound_name"] for s in approved_pareto],
                hovertemplate="<b>%{text}</b> (Approved)<br>IC50: %{x:.2f} nM<br>QED: %{y:.3f}<extra></extra>",
            ))

        fig_pareto.update_layout(
            **_BASE_LAYOUT,
            height=560,
            xaxis_title="IC50 (nM)  [log scale, lower = more potent]",
            yaxis_title="QED Drug-likeness  [higher = more drug-like]",
            xaxis_type="log",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="center", x=0.5, title=None),
        )
        fig_pareto.add_annotation(
            x=math.log10(max(0.01, min(s["ic50_nm"] for s in mo_scores))),
            y=max(s["qed"] for s in mo_scores),
            text="Ideal", showarrow=True,
            ax=40, ay=20, arrowhead=2, arrowcolor=COLORS["success"],
            font=dict(size=11, color=COLORS["success"]),
        )
        st.plotly_chart(fig_pareto, width="stretch")

        # --- Pareto table ---
        st.markdown(f"#### {len(mo_pareto)} Pareto-Optimal Compounds (4 objectives)")
        df_pareto = pd.DataFrame([
            {
                "#": i + 1,
                "Compound": s["compound_name"],
                "IC50 (nM)": s["ic50_nm"],
                "QED": s["qed"],
                "Risk": s["risk_score"],
                "SA Score": s["sa_score"],
                "MW": s["mw"],
                "LogP": s["logp"],
                "Phase": s["phase"],
            }
            for i, s in enumerate(mo_pareto)
        ])
        st.dataframe(df_pareto.style.format({
            "IC50 (nM)": lambda x: f"{x:.3f}" if x < 1 else (f"{x:.1f}" if x < 100 else f"{x:.0f}"),
            "QED": "{:.3f}",
            "Risk": "{:.2f}",
            "SA Score": "{:.1f}",
            "MW": "{:.0f}",
            "LogP": "{:.2f}",
        }), hide_index=True, width="stretch")

        # --- Objective trade-off radar for top 5 ---
        st.markdown("#### Objective Trade-offs (Top 5 Pareto Compounds)")
        categories = ["Potency\n(IC50)", "Drug-likeness\n(QED)", "MW\nOptimal", "Easy\nSynthesis", "Low\nLogP"]
        fig_radar = go.Figure()
        # Normalize within Pareto set for visible differentiation
        p_ic50s = [math.log10(max(0.01, s["ic50_nm"])) for s in mo_pareto[:5]]
        p_min_ic50, p_max_ic50 = min(p_ic50s), max(p_ic50s) if len(p_ic50s) > 1 else (min(p_ic50s), min(p_ic50s) + 1)
        ic50_range = max(0.01, p_max_ic50 - p_min_ic50)

        for s in mo_pareto[:5]:
            log_ic50 = math.log10(max(0.01, s["ic50_nm"]))
            potency = max(0, min(1, 1 - (log_ic50 - p_min_ic50) / ic50_range))
            qed = s["qed"]
            mw_score = max(0, 1 - abs(s["mw"] - 350) / 400)  # optimal ~350
            sa_score = max(0, 1 - (s["sa_score"] - 1) / 9)
            logp_score = max(0, 1 - abs(s["logp"] - 2.5) / 5)  # optimal ~2.5
            values = [potency, qed, mw_score, sa_score, logp_score]
            fig_radar.add_trace(go.Scatterpolar(
                r=values + [values[0]],
                theta=categories + [categories[0]],
                name=s["compound_name"][:15],
                fill="toself",
                opacity=0.3,
            ))
        fig_radar.update_layout(
            **_BASE_LAYOUT,
            height=420,
            polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
            showlegend=True,
            legend=dict(orientation="h", yanchor="bottom", y=-0.15, xanchor="center", x=0.5),
        )
        st.plotly_chart(fig_radar, width="stretch")


# ===========================================================================
# TAB 3 — Interaction Network
# ===========================================================================
with tab_network:
    if not display_recs:
        st.warning("No compounds to visualize.")
    else:
        st.markdown("### Compound-Target Interaction Network")
        st.markdown(
            f'<div class="info-box">'
            f"Each line connects a compound to <strong>{selected_target_name}</strong>. "
            f"Line thickness indicates binding strength (thicker = lower IC50 = stronger binding). "
            f"Node color indicates confidence level. Distance from center reflects the optimization weight "
            f"(closer = better candidate)."
            f"</div>",
            unsafe_allow_html=True,
        )

        net_n = len(display_recs)
        node_x, node_y, node_text, node_color, node_size = [], [], [], [], []
        edge_x, edge_y, edge_widths = [], [], []

        # Target at center
        node_x.append(0.0)
        node_y.append(0.0)
        node_text.append(f"<b>{selected_target_name}</b><br>(target protein)")
        node_color.append(COLORS["primary"])
        node_size.append(40)

        # Weight range for normalization
        weights = [r.weight for r in display_recs]
        w_min, w_max = min(weights), max(weights)
        w_range = w_max - w_min if w_max > w_min else 1.0

        for i, rec in enumerate(display_recs):
            angle = 2 * math.pi * i / net_n
            # Normalize radius: better compounds closer to center
            norm_w = (rec.weight - w_min) / w_range
            radius = 0.6 + norm_w * 1.4  # 0.6 to 2.0

            cx = radius * math.cos(angle)
            cy = radius * math.sin(angle)

            node_x.append(cx)
            node_y.append(cy)
            node_text.append(
                f"<b>{rec.compound_name}</b><br>"
                f"IC50: {rec.ic50_nm:.2f} nM<br>"
                f"Toxicity: {rec.toxicity:.3f}<br>"
                f"Confidence: {rec.confidence}<br>"
                f"Rank: #{rec.rank}"
            )
            node_color.append(CONFIDENCE_COLORS.get(rec.confidence, "#95a5a6"))
            node_size.append(max(10, 28 - rec.rank))

            # Edge width: inverse IC50 (thicker = more potent)
            ic50_log = math.log10(max(rec.ic50_nm, 0.01))
            edge_w = max(0.5, 4.0 - ic50_log)
            edge_x.extend([cx, 0.0, None])
            edge_y.extend([cy, 0.0, None])
            edge_widths.append(edge_w)

        fig_net = go.Figure()

        # Draw edges individually for variable width
        for i in range(net_n):
            fig_net.add_trace(go.Scatter(
                x=[edge_x[i * 3], edge_x[i * 3 + 1]],
                y=[edge_y[i * 3], edge_y[i * 3 + 1]],
                mode="lines",
                line=dict(width=edge_widths[i], color="#D0D5DD"),
                hoverinfo="none",
                showlegend=False,
            ))

        # Compound nodes
        fig_net.add_trace(go.Scatter(
            x=node_x[1:], y=node_y[1:],
            mode="markers+text",
            marker=dict(
                size=node_size[1:],
                color=node_color[1:],
                line=dict(width=1.5, color="white"),
            ),
            text=[t.split("<br>")[0].replace("<b>", "").replace("</b>", "") for t in node_text[1:]],
            textposition="top center",
            textfont=dict(size=9, color="#555"),
            hovertext=node_text[1:],
            hoverinfo="text",
            name="Compounds",
        ))

        # Target node (on top)
        fig_net.add_trace(go.Scatter(
            x=[node_x[0]], y=[node_y[0]],
            mode="markers+text",
            marker=dict(size=node_size[0], color=COLORS["primary"],
                        line=dict(width=2, color="white"), symbol="square"),
            text=[selected_target_name],
            textposition="bottom center",
            textfont=dict(size=13, color=COLORS["primary"]),
            hovertext=[node_text[0]],
            hoverinfo="text",
            name="Target",
        ))

        fig_net.update_layout(
            **_BASE_LAYOUT,
            height=620,
            showlegend=False,
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False,
                       scaleanchor="x", scaleratio=1),
        )
        st.plotly_chart(fig_net, width="stretch")

        # Legend
        st.markdown(
            f'<div class="legend-row">'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["high"]}"></span> HIGH confidence</span>'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["medium"]}"></span> MEDIUM confidence</span>'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["low"]}"></span> LOW confidence</span>'
            f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["primary"]}; border-radius:2px"></span> Target protein</span>'
            f"</div>",
            unsafe_allow_html=True,
        )


# ===========================================================================
# TAB 4 — Selectivity Profile
# ===========================================================================
with tab_selectivity:
    if not all_recs:
        st.warning("No compounds to analyze.")
    else:
        st.markdown("### Selectivity Profile")
        st.markdown(
            f'<div class="info-box">'
            f"Selectivity measures how specifically a compound binds to <strong>{selected_target_name}</strong> "
            f"versus other targets. A <strong>selectivity ratio > 10x</strong> means the compound is at least "
            f"10 times more potent against the primary target than the off-target. Higher is better."
            f"</div>",
            unsafe_allow_html=True,
        )

        # Data coverage warning
        multi_target_count = sum(
            1 for r in all_recs[:50]
            if sum(1 for _ in G.successors(r.compound_id)) > 1
        )
        if multi_target_count < len(all_recs[:50]) // 2:
            st.warning(
                f"Only **{multi_target_count}** of the top 50 compounds have multi-target data. "
                f"Most compounds were tested against {selected_target_name} only, "
                f"so selectivity cannot be assessed for the majority. "
                f"Select compounds with off-target data below for meaningful analysis."
            )

        # Compound selector — shared key for cross-tab use
        compound_options_sel = {r.compound_name: r.compound_id for r in all_recs[:50]}
        selected_compound_name = st.selectbox(
            "Select a compound to analyze",
            list(compound_options_sel.keys()),
            key="selectivity_compound",
            help="Choose from the top 50 compounds for this target.",
        )
        compound_id_sel = compound_options_sel[selected_compound_name]

        sel = compute_selectivity(G, compound_id_sel, selected_target_id)
        if "error" in sel:
            st.error(sel["error"])
        else:
            # KPI row
            s1, s2, s3 = st.columns(3)
            s1.metric(f"IC50 vs {selected_target_name}", f"{sel['primary_ic50_nm']:.2f} nM")
            s2.metric("Off-Targets Detected", sel["off_target_count"])
            selective_label = "Selective" if sel["is_selective"] else "Non-selective"
            s3.metric("Selectivity Status", selective_label,
                       delta="All ratios > 10x" if sel["is_selective"] else "Some ratios < 10x",
                       delta_color="normal" if sel["is_selective"] else "inverse")

            if sel["off_targets"]:
                chart_col, table_col = st.columns([3, 2])

                with chart_col:
                    df_sel = pd.DataFrame(sel["off_targets"])

                    # Color based on selectivity threshold
                    bar_colors = [
                        COLORS["high"] if r > 100 else (COLORS["medium"] if r > 10 else COLORS["low"])
                        for r in df_sel["selectivity_ratio"]
                    ]

                    fig_sel = go.Figure()
                    fig_sel.add_trace(go.Bar(
                        x=df_sel["target_name"],
                        y=df_sel["selectivity_ratio"],
                        marker_color=bar_colors,
                        text=[f"{r:.1f}x" for r in df_sel["selectivity_ratio"]],
                        textposition="outside",
                        hovertemplate=(
                            "<b>%{x}</b><br>"
                            "IC50: %{customdata:.1f} nM<br>"
                            "Selectivity: %{y:.1f}x<extra></extra>"
                        ),
                        customdata=df_sel["ic50_nm"],
                    ))
                    fig_sel.add_hline(
                        y=10, line_dash="dash", line_color=COLORS["success"], line_width=2,
                        annotation_text="10x threshold (selective)",
                        annotation_position="top right",
                        annotation_font_color=COLORS["success"],
                    )
                    fig_sel.update_layout(
                        **_BASE_LAYOUT,
                        height=400,
                        xaxis_title="Off-Target Protein",
                        yaxis_title="Selectivity Ratio (fold difference)",
                        yaxis_type="log" if max(df_sel["selectivity_ratio"]) > 100 else "linear",
                        showlegend=False,
                    )
                    st.plotly_chart(fig_sel, width="stretch")

                with table_col:
                    st.markdown("#### Off-Target Details")
                    df_display = df_sel[["target_name", "ic50_nm", "selectivity_ratio"]].copy()
                    df_display.columns = ["Off-Target", "IC50 (nM)", "Ratio"]
                    df_display["Risk"] = df_display["Ratio"].apply(
                        lambda r: "Low" if r > 100 else ("Moderate" if r > 10 else "High")
                    )
                    st.dataframe(df_display, hide_index=True, width="stretch")

                st.markdown(
                    f'<div class="legend-row">'
                    f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["high"]}"></span> > 100x (Low risk)</span>'
                    f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["medium"]}"></span> 10-100x (Moderate)</span>'
                    f'<span class="legend-item"><span class="legend-dot" style="background:{COLORS["low"]}"></span> < 10x (High risk)</span>'
                    f"</div>",
                    unsafe_allow_html=True,
                )
            else:
                st.info(
                    f"**{selected_compound_name}** was only tested against {selected_target_name}. "
                    f"No off-target binding data is available for selectivity analysis."
                )


# ===========================================================================
# TAB 5 — Drug-likeness (Molecular Properties)
# ===========================================================================
with tab_properties:
    if not all_recs:
        st.warning("No compounds.")
    else:
        st.markdown("### Drug-likeness Assessment")
        st.markdown(
            f'<div class="info-box">'
            f"Drug-likeness is estimated using <strong>Lipinski's Rule of Five</strong>, which predicts oral "
            f"bioavailability based on molecular weight, hydrogen bond donors/acceptors, and LogP. "
            f"Compounds violating 2+ rules are less likely to be orally bioavailable."
            f"</div>",
            unsafe_allow_html=True,
        )

        compound_options_prop = {r.compound_name: r.compound_id for r in all_recs[:50]}
        selected_comp_prop = st.selectbox(
            "Select a compound to analyze",
            list(compound_options_prop.keys()),
            key="properties_compound",
            help="Choose from the top 50 compounds for this target.",
        )
        comp_id_prop = compound_options_prop[selected_comp_prop]

        smiles = G.nodes[comp_id_prop].get("smiles", "")
        if not smiles:
            st.warning("No SMILES data available for this compound.")
        else:
            props = compute_properties(smiles)

            # KPI row with context
            p1, p2, p3, p4, p5 = st.columns(5)
            p1.metric("Est. MW", f"{props.estimated_mw:.0f}", delta="OK" if props.estimated_mw <= 500 else "High",
                       delta_color="normal" if props.estimated_mw <= 500 else "inverse")
            p2.metric("Est. LogP", f"{props.estimated_logp:.1f}", delta="OK" if props.estimated_logp <= 5 else "High",
                       delta_color="normal" if props.estimated_logp <= 5 else "inverse")
            p3.metric("HBD", props.estimated_hbd, delta="OK" if props.estimated_hbd <= 5 else "High",
                       delta_color="normal" if props.estimated_hbd <= 5 else "inverse")
            p4.metric("HBA", props.estimated_hba, delta="OK" if props.estimated_hba <= 10 else "High",
                       delta_color="normal" if props.estimated_hba <= 10 else "inverse")
            p5.metric("Violations", f"{props.lipinski_violations}/4",
                       delta="Drug-like" if props.lipinski_violations <= 1 else "Risky",
                       delta_color="normal" if props.lipinski_violations <= 1 else "inverse")

            st.code(f"SMILES: {smiles}", language=None)

            radar_col, detail_col = st.columns([3, 2])

            with radar_col:
                # Radar chart
                categories = ["MW / 500", "HBD / 5", "HBA / 10", "LogP / 5", "RotBonds / 10"]
                values = [
                    props.estimated_mw / 500,
                    props.estimated_hbd / 5,
                    props.estimated_hba / 10,
                    props.estimated_logp / 5 if props.estimated_logp > 0 else 0,
                    props.estimated_rotatable / 10,
                ]
                # Close the polygon
                categories_closed = categories + [categories[0]]
                values_closed = values + [values[0]]
                limit_closed = [1] * len(categories) + [1]

                fig_radar = go.Figure()
                fig_radar.add_trace(go.Scatterpolar(
                    r=values_closed,
                    theta=categories_closed,
                    fill="toself",
                    fillcolor="rgba(0,102,255,0.12)",
                    line=dict(color=COLORS["primary"], width=2),
                    name=selected_comp_prop,
                ))
                fig_radar.add_trace(go.Scatterpolar(
                    r=limit_closed,
                    theta=categories_closed,
                    fill="none",
                    line=dict(color=COLORS["danger"], width=2, dash="dash"),
                    name="Lipinski Limit",
                ))
                fig_radar.update_layout(
                    **_BASE_LAYOUT,
                    polar=dict(
                        radialaxis=dict(visible=True, range=[0, max(2, max(values) * 1.1)],
                                         showticklabels=True, tickfont_size=9),
                    ),
                    height=400,
                    showlegend=True,
                    legend=dict(orientation="h", yanchor="bottom", y=-0.15, xanchor="center", x=0.5),
                )
                st.plotly_chart(fig_radar, width="stretch")

            with detail_col:
                st.markdown("#### Property Details")
                detail_data = {
                    "Property": ["Molecular Weight", "H-Bond Donors", "H-Bond Acceptors", "LogP", "Rotatable Bonds"],
                    "Value": [f"{props.estimated_mw:.0f}", str(props.estimated_hbd),
                              str(props.estimated_hba), f"{props.estimated_logp:.1f}",
                              str(props.estimated_rotatable)],
                    "Limit": ["<= 500", "<= 5", "<= 10", "<= 5", "<= 10"],
                    "Status": [
                        "Pass" if props.estimated_mw <= 500 else "Fail",
                        "Pass" if props.estimated_hbd <= 5 else "Fail",
                        "Pass" if props.estimated_hba <= 10 else "Fail",
                        "Pass" if props.estimated_logp <= 5 else "Fail",
                        "Pass" if props.estimated_rotatable <= 10 else "Warn",
                    ],
                }
                df_detail = pd.DataFrame(detail_data)

                def color_status(val):
                    if val == "Pass":
                        return f"color: {COLORS['high']}; font-weight: 600"
                    elif val == "Fail":
                        return f"color: {COLORS['low']}; font-weight: 600"
                    return f"color: {COLORS['medium']}; font-weight: 600"

                styled_detail = df_detail.style.map(color_status, subset=["Status"])
                st.dataframe(styled_detail, hide_index=True, width="stretch")

                # Risk assessment
                st.markdown("#### Risk Assessment")
                risk_color = COLORS["high"] if props.risk_score < 0.2 else (
                    COLORS["medium"] if props.risk_score < 0.5 else COLORS["low"])
                risk_label = "Low Risk" if props.risk_score < 0.2 else (
                    "Moderate Risk" if props.risk_score < 0.5 else "High Risk")
                st.markdown(f"**Overall Risk Score:** <span style='color:{risk_color}; font-weight:700'>"
                            f"{props.risk_score:.2f} ({risk_label})</span>", unsafe_allow_html=True)

                if props.risk_flags:
                    for flag in props.risk_flags:
                        st.markdown(f"- {flag}")
                else:
                    st.markdown("_No structural risk flags detected._")


# ===========================================================================
# TAB 6 — Retrospective Validation
# ===========================================================================
with tab_validation:
    st.markdown("### Retrospective Validation: Does ChemPath Find Approved Drugs?")
    st.markdown(
        f'<div class="info-box">'
        f"<strong>Ground truth test:</strong> Our dataset contains {sum(1 for c in data['compounds'] if c.get('phase', 0) >= 4)} "
        f"FDA-approved drugs (Phase 4). If ChemPath's ranking is meaningful, these approved drugs "
        f"should appear in the top rankings more often than random chance. "
        f"The <strong>Enrichment Factor (EF)</strong> measures this: EF > 1.0 means better than random."
        f"</div>",
        unsafe_allow_html=True,
    )

    # Run validation — use top 10% of compounds as window
    recs_for_val = rank_compounds_for_target(G, selected_target_id, strategy=strategy)
    n_total = len(recs_for_val)
    val_top_n = max(10, n_total // 10)  # top 10% or at least 10
    val_result = validate_target(G, selected_target_id, top_n=val_top_n, strategy=strategy)

    if val_result is None:
        st.info(f"No approved drugs (Phase 4+) found among compounds tested against {selected_target_name}.")
    else:
        # Find percentile ranks for all approved drugs
        approved_ranks = []
        for rec in recs_for_val:
            phase = G.nodes.get(rec.compound_id, {}).get("phase", 0)
            if phase >= 4:
                percentile = round((1 - rec.rank / n_total) * 100, 1)
                approved_ranks.append({
                    "name": rec.compound_name,
                    "rank": rec.rank,
                    "percentile": percentile,
                    "ic50": rec.ic50_nm,
                    "phase": phase,
                })

        best_percentile = max(a["percentile"] for a in approved_ranks) if approved_ranks else 0

        # KPIs
        vc1, vc2, vc3, vc4, vc5 = st.columns(5)
        vc1.metric("Approved Drugs", len(approved_ranks))
        vc2.metric(f"In Top {val_top_n}", len(val_result.approved_in_top_n))
        vc3.metric("Enrichment Factor", f"{val_result.enrichment_factor:.1f}x")
        vc4.metric("Best Percentile", f"Top {100 - best_percentile:.0f}%")
        vc5.metric("Total Compounds", n_total)

        if val_result.enrichment_factor > 1:
            st.success(
                f"ChemPath is **{val_result.enrichment_factor:.1f}x better than random** at finding "
                f"approved drugs for {selected_target_name} in the top {val_top_n} (top {val_top_n*100//n_total}%)."
            )
        else:
            # Check if approved drugs are at least in top 50%
            in_top_half = sum(1 for a in approved_ranks if a["percentile"] >= 50)
            if in_top_half == len(approved_ranks):
                st.info(
                    f"All {len(approved_ranks)} approved drugs rank in the **top {100 - min(a['percentile'] for a in approved_ranks):.0f}%** "
                    f"of {n_total} compounds. They don't reach top {val_top_n} because {val_top_n} compounds have even lower IC50 values. "
                    f"This is expected: approved drugs are optimized for safety/PK, not just potency."
                )
            else:
                st.warning(
                    f"Approved drugs for {selected_target_name} have moderate IC50 values. "
                    f"They were approved for properties (safety, selectivity, PK) not captured in our IC50-only data."
                )

        # Approved drugs with percentile ranks
        st.markdown("#### Where Do Approved Drugs Rank?")
        df_approved = pd.DataFrame(approved_ranks)
        df_approved = df_approved.sort_values("rank")
        df_approved.columns = ["Drug Name", "Rank", "Percentile", "IC50 (nM)", "Phase"]
        st.dataframe(df_approved.style.format({
            "IC50 (nM)": lambda x: f"{x:.2f}" if x < 10 else f"{x:.0f}",
            "Percentile": lambda x: f"Top {100-x:.0f}%",
        }), hide_index=True, width="stretch")

        # Visual: rank position bar chart
        st.markdown("#### Rank Distribution")
        fig_ranks = go.Figure()
        # All compounds as background
        fig_ranks.add_trace(go.Bar(
            x=list(range(1, min(n_total + 1, 101))),
            y=[1] * min(n_total, 100),
            marker_color="#E8EAF0",
            name="All compounds",
            hoverinfo="skip",
        ))
        # Approved drugs as highlighted bars
        for a in approved_ranks:
            bar_pos = min(a["rank"], 100)
            fig_ranks.add_trace(go.Bar(
                x=[bar_pos],
                y=[1.5],
                marker_color=COLORS["success"] if a["rank"] <= val_top_n else COLORS["warning"],
                name=a["name"],
                text=a["name"],
                hovertemplate=f"<b>{a['name']}</b><br>Rank: #{a['rank']}/{n_total}<br>IC50: {a['ic50']:.1f} nM<extra></extra>",
            ))
        fig_ranks.update_layout(
            **_BASE_LAYOUT,
            height=200,
            showlegend=False,
            xaxis_title=f"Rank (1 = best, {n_total} = worst)",
            yaxis_visible=False,
            bargap=0,
        )
        fig_ranks.add_annotation(
            x=val_top_n, y=1.8,
            text=f"Top {val_top_n*100//n_total}% cutoff",
            showarrow=True, arrowhead=2, ay=-30,
            font=dict(size=10, color=COLORS["primary"]),
        )
        st.plotly_chart(fig_ranks, width="stretch")

    # Cross-target validation summary
    st.markdown("---")
    st.markdown("#### Cross-Target Validation Summary")
    st.caption("Enrichment factor (top 10%) across all targets with approved drugs")

    targets_with_data = []
    for tid, tname in [(t["chembl_id"], t["name"]) for t in data["targets"]]:
        t_recs = rank_compounds_for_target(G, tid, strategy=strategy)
        t_top_n = max(10, len(t_recs) // 10)  # consistent top 10%
        vr = validate_target(G, tid, top_n=t_top_n, strategy=strategy)
        if vr is not None:
            targets_with_data.append({
                "Target": vr.target_name,
                "Total Compounds": vr.total_compounds,
                "Approved Drugs": len(vr.approved_compounds),
                f"In Top 10%": len(vr.approved_in_top_n),
                "Enrichment Factor": vr.enrichment_factor,
                "Hit Rate": vr.hit_rate,
            })

    if targets_with_data:
        df_cross = pd.DataFrame(targets_with_data)
        avg_ef = df_cross["Enrichment Factor"].mean()
        def _color_ef(val):
            if val >= 5:
                return f"background-color: {COLORS['high']}33; color: {COLORS['high']}; font-weight:600"
            elif val >= 1:
                return f"background-color: {COLORS['medium']}33; color: {COLORS['medium']}; font-weight:600"
            else:
                return f"background-color: {COLORS['low']}33; color: {COLORS['low']}; font-weight:600"

        st.dataframe(df_cross.style.format({
            "Enrichment Factor": "{:.1f}x",
            "Hit Rate": "{:.0%}",
        }).map(_color_ef, subset=["Enrichment Factor"]),
            hide_index=True, width="stretch")
        st.markdown(f"**Average Enrichment Factor: {avg_ef:.1f}x** (>1.0 = better than random)")

        # Interpretation
        st.markdown("---")
        st.markdown("#### Why Some Approved Drugs Don't Rank High")
        st.markdown("""
**This is expected and reveals an important truth about drug discovery:**

1. **Approved drugs aren't always the most potent.** A drug with IC50=50nM but excellent safety, selectivity,
   and pharmacokinetics beats a compound with IC50=0.1nM but poor ADMET properties.

2. **Our data is IC50-only.** Real drug approval considers: selectivity panels, metabolic stability,
   CYP inhibition, hERG liability, oral bioavailability, in-vivo efficacy, and clinical trial outcomes.

3. **ChemPath ranks by what it has.** With only IC50 + estimated toxicity, it correctly finds the most potent
   compounds. Approved drugs that rank lower were approved for properties we can't measure from SMILES alone.

**What this means for users:** The top ChemPath rankings are valid starting points for hit identification,
but candidates should be validated through ADMET assays, selectivity screening, and dose-response studies
before advancing to lead optimization.
""")
    else:
        st.info("No targets in the current dataset have approved drugs for validation.")


# ===========================================================================
# TAB 7 — Methodology & Benchmarks
# ===========================================================================
with tab_methods:
    st.markdown("### How ChemPath Works")
    st.markdown(
        f'<div class="info-box">'
        f"ChemPath transforms drug screening into a <strong>shortest-path graph optimization</strong> problem. "
        f"Compounds and targets become graph nodes; bioactivity measurements (IC50) become directed edges "
        f"with weights derived from a pIC50 transformation plus a tunable toxicity penalty. "
        f"Ranking compounds then reduces to finding <strong>minimum-weight paths</strong> from compounds to targets."
        f"</div>",
        unsafe_allow_html=True,
    )

    # --- Pipeline diagram ---
    st.markdown("#### Optimization Pipeline")
    pipe_cols = st.columns(5)
    pipeline_steps = [
        ("1. Data Ingestion", "ChEMBL bioactivity data\n(IC50, targets, SMILES)"),
        ("2. Graph Construction", "Compounds + targets as nodes\nIC50 -> edge weights via pIC50"),
        ("3. Weight Optimization", "weight = (14 - pIC50)\n+ toxicity_penalty * tox_risk"),
        ("4. Ranking / Pareto", "Sort by weight (balanced)\nor multi-objective Pareto front"),
        ("5. Confidence Labeling", "4-dimensional scoring:\ndata source, phase, potency, tox"),
    ]
    for col, (title, desc) in zip(pipe_cols, pipeline_steps):
        col.markdown(f"**{title}**")
        col.caption(desc)

    st.markdown("---")

    # --- Key design choices ---
    st.markdown("#### Key Design Choices")
    design_left, design_right = st.columns(2)

    with design_left:
        st.markdown("##### Weight Function")
        st.markdown(
            "The IC50-to-weight transform uses `weight = 14 - pIC50` (clamped to [0.01, 10.0]), "
            "where `pIC50 = -log10(IC50 * 1e-9)`. This maps potent binders (low IC50) to low weights, "
            "making shortest-path = most potent compound.\n\n"
            "**Toxicity penalty** adds `tox_risk * penalty_factor` to each edge weight, "
            "creating a tunable efficacy-safety tradeoff controlled by the sidebar slider."
        )

    with design_right:
        st.markdown("##### Confidence System")
        st.markdown(
            "Each recommendation gets a confidence label (HIGH/MEDIUM/LOW) based on four dimensions:\n\n"
            "| Dimension | HIGH score | LOW score |\n"
            "|-----------|-----------|----------|\n"
            "| Data source | Experimental IC50 (+3) | Predicted (+0) |\n"
            "| Clinical phase | Approved drug (+3) | No trial data (+0) |\n"
            "| Potency | IC50 < 10 nM (+2) | IC50 > 100 nM (+0) |\n"
            "| Toxicity | Risk < 0.2 (+2) | Risk > 0.35 (+0) |"
        )

    st.markdown("---")

    # --- SOTA Comparison ---
    st.markdown("### Comparison with State-of-the-Art Methods")
    st.markdown(
        f'<div class="info-box">'
        f"This section compares ChemPath's graph-based approach against established and cutting-edge "
        f"drug discovery methods. Each method has distinct strengths; the right choice depends on "
        f"your screening stage, available data, and computational budget."
        f"</div>",
        unsafe_allow_html=True,
    )

    # Comparison data
    methods_data = [
        {
            "Method": "ChemPath (This Tool)",
            "Category": "Graph Optimization",
            "Input": "IC50 bioactivity + SMILES",
            "Output": "Ranked compounds + Pareto front",
            "DUD-E AUC": "N/A",
            "LIT-PCBA AUC": "N/A",
            "Speed": "< 1 sec",
            "Key Strength": "Multi-objective, interpretable, real-time",
            "Key Limitation": "Requires existing IC50 data",
        },
        {
            "Method": "AutoDock Vina",
            "Category": "Structure-based Docking",
            "Input": "3D protein + ligand structures",
            "Output": "Binding pose + score",
            "DUD-E AUC": "~0.68",
            "LIT-PCBA AUC": "0.55-0.65",
            "Speed": "~min/compound",
            "Key Strength": "Physics-based, no training data needed",
            "Key Limitation": "Needs 3D structures, high false-positive rate",
        },
        {
            "Method": "Glide (Schrodinger)",
            "Category": "Structure-based Docking",
            "Input": "3D protein + ligand structures",
            "Output": "Docking score + pose",
            "DUD-E AUC": "~0.70",
            "LIT-PCBA AUC": "0.55-0.65",
            "Speed": "~min/compound",
            "Key Strength": "Industry standard, well-validated",
            "Key Limitation": "Commercial license, needs crystal structures",
        },
        {
            "Method": "Chemprop (D-MPNN)",
            "Category": "Graph Neural Network",
            "Input": "SMILES strings",
            "Output": "Predicted bioactivity",
            "DUD-E AUC": "~0.85+",
            "LIT-PCBA AUC": "~0.65-0.75",
            "Speed": "ms/compound",
            "Key Strength": "SOTA on molecular benchmarks (HIV ROC-AUC 0.81)",
            "Key Limitation": "Needs labeled training data, black-box",
        },
        {
            "Method": "TxGNN",
            "Category": "Biomedical Knowledge Graph",
            "Input": "Disease-gene-drug KG",
            "Output": "Drug repurposing predictions",
            "DUD-E AUC": "N/A",
            "LIT-PCBA AUC": "N/A",
            "Speed": "Seconds",
            "Key Strength": "49.2% improvement zero-shot (Nature Med 2024)",
            "Key Limitation": "Requires large KG, no dose-response info",
        },
        {
            "Method": "RF-Score-VS / XGBoost",
            "Category": "Classical ML",
            "Input": "Molecular fingerprints",
            "Output": "Active/inactive classification",
            "DUD-E AUC": "0.80-0.92",
            "LIT-PCBA AUC": "0.59-0.75",
            "Speed": "ms/compound",
            "Key Strength": "Fast, interpretable, 55.6% hit rate at top 1%",
            "Key Limitation": "Fingerprint-dependent, limited chemical space",
        },
        {
            "Method": "Deffini (family DNN)",
            "Category": "Deep Neural Network",
            "Input": "Molecular descriptors",
            "Output": "Activity prediction",
            "DUD-E AUC": "0.92",
            "LIT-PCBA AUC": "Varies",
            "Speed": "ms/compound",
            "Key Strength": "Highest single-method AUC on DUD-E",
            "Key Limitation": "Family-specific, DUD-E bias inflates results",
        },
        {
            "Method": "LLM-based (DrugGPT/Gen)",
            "Category": "Large Language Model",
            "Input": "Text / SMILES / protein seq",
            "Output": "Generated molecules",
            "DUD-E AUC": "N/A",
            "LIT-PCBA AUC": "N/A",
            "Speed": "Seconds",
            "Key Strength": "Generative, broad chemical knowledge",
            "Key Limitation": "Hallucination risk, not clinically validated",
        },
    ]

    df_methods = pd.DataFrame(methods_data)

    # Styled comparison table
    def highlight_chempath(row):
        if row["Method"].startswith("ChemPath"):
            return [f"background-color: rgba(0,102,255,0.06); font-weight: 600"] * len(row)
        return [""] * len(row)

    styled_methods = df_methods.style.apply(highlight_chempath, axis=1)
    st.dataframe(styled_methods, hide_index=True, width="stretch", height=360)

    # --- Benchmark context ---
    st.markdown("---")
    st.markdown("#### Understanding the Benchmarks")
    bench_left, bench_right = st.columns(2)

    with bench_left:
        st.markdown(
            f'<div class="info-box">'
            f"<strong>DUD-E</strong> (Directory of Useful Decoys - Enhanced): 102 targets, "
            f"66,695 actives from ChEMBL with property-matched decoys. "
            f"<strong>Caveat:</strong> Known to have analog bias — ML models can achieve high AUC by "
            f"learning physicochemical differences between actives/decoys rather than true binding signals. "
            f"Results on DUD-E alone are insufficient for credible benchmarking."
            f"</div>",
            unsafe_allow_html=True,
        )

    with bench_right:
        st.markdown(
            f'<div class="info-box">'
            f"<strong>LIT-PCBA</strong>: 15 targets, 7,844 confirmed actives, 407,381 confirmed inactives "
            f"from dose-response PubChem assays. Uses asymmetric validation embedding (AVE) to reduce "
            f"bias. <strong>Considered the most rigorous current benchmark</strong> — methods that "
            f"perform well on DUD-E often struggle here. Realistic hit rates mimic actual screening."
            f"</div>",
            unsafe_allow_html=True,
        )

    st.caption(
        "Standard evaluation metrics: AUC-ROC (ranking quality), BEDROC alpha=20 (early enrichment), "
        "EF at 1%/5% (enrichment factor), AUC-PRC (for imbalanced data). "
        "AUC-ROC alone is insufficient — BEDROC and EF at 1% are more practically relevant because "
        "only top-ranked compounds proceed to experimental testing."
    )

    st.markdown("---")

    # --- Detailed advantage/disadvantage analysis ---
    st.markdown("#### Detailed Comparison")

    cmp_left, cmp_right = st.columns(2)

    with cmp_left:
        st.markdown("##### ChemPath Advantages")
        st.markdown("""
- **Multi-objective optimization**: Jointly optimizes efficacy and safety via Pareto front, unlike single-score methods (docking, ML classifiers)
- **Interpretable rankings**: Every weight and confidence label is traceable to concrete data (IC50, clinical phase, toxicity risk) — no black box
- **Real-time interactive**: Rankings update instantly when parameters change; docking takes minutes per compound
- **Sensitivity analysis**: Built-in robustness testing (IC50 noise, toxicity penalty, missing data) — most tools lack this
- **No 3D structures needed**: Works with IC50 + SMILES, while docking requires protein crystal structures
- **Selectivity profiling**: Cross-target comparison is native to the graph; docking methods screen one target at a time
- **Confidence calibration**: 4-dimensional confidence scoring (source, phase, potency, toxicity) provides actionable uncertainty estimates
""")

    with cmp_right:
        st.markdown("##### ChemPath Limitations")
        st.markdown("""
- **Requires existing bioactivity data**: Cannot screen novel compounds with no IC50 data (docking and ML can predict for unseen molecules)
- **No 3D binding mode**: Cannot predict binding poses or specific interactions (docking provides atomic-level detail)
- **Simplified toxicity model**: Estimates toxicity from SMILES-derived Lipinski properties, not from ADMET assays or in-vivo data
- **No dose-response modeling**: Uses single IC50 values, not full dose-response curves with Hill slopes
- **Approximate molecular properties**: SMILES-based MW/LogP/HBD/HBA estimates are less accurate than RDKit or experimental values
- **No PAINS filtering**: Does not flag pan-assay interference compounds that may give false positives
- **No generative capability**: Cannot design novel molecules; limited to ranking existing compounds
""")

    st.markdown("---")

    # --- Real-world screening coverage ---
    st.markdown("### Real-World Drug Screening Coverage")
    st.markdown(
        f'<div class="info-box">'
        f"This section maps ChemPath's current capabilities to the stages of a real drug screening workflow, "
        f"highlighting what is covered and what requires additional tools or experimental validation."
        f"</div>",
        unsafe_allow_html=True,
    )

    coverage_data = [
        {
            "Screening Stage": "Target selection",
            "Real-World Need": "Validate target with genetic/clinical evidence",
            "ChemPath Coverage": "Partial",
            "Details": "Provides curated anticancer targets from ChEMBL; does not validate target-disease association",
        },
        {
            "Screening Stage": "Hit identification",
            "Real-World Need": "Screen large library, identify active compounds",
            "ChemPath Coverage": "Full",
            "Details": "Ranks all compounds with IC50 data against target using graph optimization",
        },
        {
            "Screening Stage": "Hit-to-lead prioritization",
            "Real-World Need": "Multi-parameter optimization (potency, selectivity, ADMET)",
            "ChemPath Coverage": "Partial",
            "Details": "Potency + toxicity Pareto optimization and selectivity profiling; no ADMET prediction",
        },
        {
            "Screening Stage": "Dose-response confirmation",
            "Real-World Need": "Full dose-response curves, Hill slope, EC50/IC50",
            "ChemPath Coverage": "Limited",
            "Details": "Uses single IC50 values from ChEMBL; does not model full dose-response curves",
        },
        {
            "Screening Stage": "Selectivity screening",
            "Real-World Need": "Panel screening against related and unrelated targets",
            "ChemPath Coverage": "Full",
            "Details": "Cross-target selectivity ratios with automatic off-target detection",
        },
        {
            "Screening Stage": "Drug-likeness assessment",
            "Real-World Need": "Lipinski Ro5, Veber rules, PAINS filter, lead-likeness",
            "ChemPath Coverage": "Partial",
            "Details": "Lipinski Ro5 from SMILES; no PAINS filter, no Veber rules, approximate property estimation",
        },
        {
            "Screening Stage": "ADMET profiling",
            "Real-World Need": "Absorption, distribution, metabolism, excretion, toxicity",
            "ChemPath Coverage": "Minimal",
            "Details": "Toxicity estimated from molecular properties; no CYP inhibition, permeability, or metabolic stability",
        },
        {
            "Screening Stage": "Clinical translation",
            "Real-World Need": "Phase I/II/III trial data, PK/PD modeling",
            "ChemPath Coverage": "Metadata only",
            "Details": "Displays clinical phase from ChEMBL (used in confidence scoring); no PK/PD modeling",
        },
    ]

    df_coverage = pd.DataFrame(coverage_data)

    def color_coverage(val):
        colors = {
            "Full": f"background-color: rgba(0,196,140,0.15); color: {COLORS['high']}; font-weight: 600",
            "Partial": f"background-color: rgba(255,176,32,0.15); color: {COLORS['medium']}; font-weight: 600",
            "Limited": f"background-color: rgba(255,77,79,0.12); color: {COLORS['low']}; font-weight: 600",
            "Minimal": f"background-color: rgba(255,77,79,0.12); color: {COLORS['low']}; font-weight: 600",
            "Metadata only": f"background-color: rgba(140,140,140,0.1); color: {COLORS['text_muted']}; font-weight: 600",
        }
        return colors.get(val, "")

    styled_coverage = df_coverage.style.map(color_coverage, subset=["ChemPath Coverage"])
    st.dataframe(styled_coverage, hide_index=True, width="stretch", height=350)

    # --- Visual: coverage radar ---
    st.markdown("#### Coverage Heatmap")
    coverage_scores = {
        "Hit identification": 1.0,
        "Selectivity screening": 0.9,
        "Hit-to-lead": 0.6,
        "Drug-likeness": 0.5,
        "Target selection": 0.4,
        "Dose-response": 0.2,
        "ADMET": 0.15,
        "Clinical translation": 0.1,
    }
    cat_names = list(coverage_scores.keys())
    cat_values = list(coverage_scores.values())
    cat_names_closed = cat_names + [cat_names[0]]
    cat_values_closed = cat_values + [cat_values[0]]

    fig_coverage = go.Figure()
    fig_coverage.add_trace(go.Scatterpolar(
        r=cat_values_closed,
        theta=cat_names_closed,
        fill="toself",
        fillcolor="rgba(0,102,255,0.12)",
        line=dict(color=COLORS["primary"], width=2),
        name="ChemPath Coverage",
    ))
    fig_coverage.add_trace(go.Scatterpolar(
        r=[1.0] * len(cat_names_closed),
        theta=cat_names_closed,
        fill="none",
        line=dict(color=COLORS["text_muted"], width=1, dash="dash"),
        name="Ideal (full coverage)",
    ))
    fig_coverage.update_layout(
        **_BASE_LAYOUT,
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 1.1], tickvals=[0.25, 0.5, 0.75, 1.0],
                             ticktext=["25%", "50%", "75%", "100%"], tickfont_size=9),
        ),
        height=450,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=-0.15, xanchor="center", x=0.5),
    )
    st.plotly_chart(fig_coverage, width="stretch")

    st.markdown("---")

    # --- When to use ChemPath vs alternatives ---
    st.markdown("#### When to Use Each Method")
    when_data = [
        {
            "Scenario": "I have IC50 data and want to prioritize compounds quickly",
            "Best Method": "ChemPath",
            "Why": "Instant multi-objective ranking with confidence labels",
        },
        {
            "Scenario": "I need to screen novel compounds with no bioactivity data",
            "Best Method": "Docking (Vina/Glide) or Chemprop",
            "Why": "Can predict binding from structure alone",
        },
        {
            "Scenario": "I want to understand how a drug binds at the atomic level",
            "Best Method": "Molecular Docking + MD Simulation",
            "Why": "Provides 3D binding poses and interaction maps",
        },
        {
            "Scenario": "I need to repurpose existing drugs for a rare disease",
            "Best Method": "TxGNN or Knowledge Graph methods",
            "Why": "Zero-shot capability for diseases with no prior screening data",
        },
        {
            "Scenario": "I want a quick filter before expensive experiments",
            "Best Method": "ChemPath + Lipinski filter",
            "Why": "Fast triage with interpretable confidence and risk flags",
        },
        {
            "Scenario": "I need to design novel molecules from scratch",
            "Best Method": "Generative models (VAE, Diffusion, LLM)",
            "Why": "Can propose new molecular structures with desired properties",
        },
        {
            "Scenario": "I want robust rankings that survive parameter changes",
            "Best Method": "ChemPath (sensitivity analysis)",
            "Why": "Built-in sensitivity analysis tests ranking stability under perturbation",
        },
    ]
    df_when = pd.DataFrame(when_data)
    st.dataframe(df_when, hide_index=True, width="stretch", height=300)
