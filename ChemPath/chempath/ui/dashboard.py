"""
ChemPath — Multi-Hop Drug Discovery Dashboard v2.1

Sidebar-driven navigation with integrated views.
Run with: uv run streamlit run chempath/ui/dashboard.py
"""

import streamlit as st
import plotly.graph_objects as go
import pandas as pd
import math
from pathlib import Path

from chempath.data.chembl_client import load_saved_data, ANTICANCER_TARGETS
from chempath.data.biological_network import DISEASES
from chempath.graph.network import build_multihop_graph, get_multihop_summary
from chempath.graph.pathfinding import (
    find_repurposing_candidates,
    find_mechanism_of_action,
    identify_targets_for_disease,
    find_combination_therapy,
)

DATA_PATH = Path(__file__).parent.parent.parent / "data" / "chembl_real.json"

# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------
C = {
    "blue": "#2563EB", "blue_bg": "#EFF6FF",
    "green": "#059669", "green_bg": "#ECFDF5",
    "amber": "#D97706", "amber_bg": "#FFFBEB",
    "red": "#DC2626", "red_bg": "#FEF2F2",
    "purple": "#7C3AED", "purple_bg": "#F5F3FF",
    "gray": "#6B7280", "gray_bg": "#F9FAFB",
    "border": "#E5E7EB", "text": "#111827",
    # Node type colors
    "compound": "#2563EB", "target": "#059669",
    "pathway": "#D97706", "disease": "#DC2626",
}

_LAYOUT = dict(
    template="plotly_white",
    font=dict(family="Inter, -apple-system, sans-serif", size=13, color=C["text"]),
    margin=dict(l=40, r=24, t=40, b=40),
    hoverlabel=dict(bgcolor="white", font_size=12, bordercolor=C["border"]),
    paper_bgcolor="white", plot_bgcolor="white",
)

NODE_COLORS = {"compound": C["compound"], "target": C["target"],
               "pathway": C["pathway"], "disease": C["disease"]}


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="ChemPath",
    page_icon="https://img.icons8.com/fluency/48/test-tube.png",
    layout="wide", initial_sidebar_state="expanded",
)

st.markdown("""<style>
.block-container{padding-top:1rem;max-width:1100px}
section[data-testid="stSidebar"]>div:first-child{padding-top:.6rem}
div[data-testid="stMetric"]{background:#F9FAFB;border:1px solid #E5E7EB;border-radius:10px;padding:12px 16px}
div[data-testid="stMetric"] label{color:#6B7280!important;font-size:.72rem!important;text-transform:uppercase;letter-spacing:.04em}
div[data-testid="stMetric"] div[data-testid="stMetricValue"]{font-size:1.3rem!important;font-weight:700!important;color:#111827!important}
.path-card{background:#F9FAFB;border:1px solid #E5E7EB;border-radius:10px;padding:14px 18px;margin:6px 0;font-size:.87rem;line-height:1.6}
.path-card-top{border-left:4px solid #059669}
.tag{display:inline-block;padding:2px 8px;border-radius:4px;font-size:.75rem;font-weight:600;margin-right:4px}
.tag-compound{background:#EFF6FF;color:#2563EB}
.tag-target{background:#ECFDF5;color:#059669}
.tag-pathway{background:#FFFBEB;color:#D97706}
.tag-disease{background:#FEF2F2;color:#DC2626}
.legend-row{display:flex;gap:16px;margin:6px 0 12px 0;flex-wrap:wrap}
.legend-item{display:flex;align-items:center;gap:5px;font-size:.78rem;color:#6B7280}
.legend-dot{width:9px;height:9px;border-radius:50%;display:inline-block}
.section-desc{color:#6B7280;font-size:.88rem;margin-bottom:12px;line-height:1.5}
</style>""", unsafe_allow_html=True)


# ---------------------------------------------------------------------------
# Data + Graph (cached)
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_data():
    return load_saved_data(DATA_PATH)

@st.cache_resource(show_spinner=False)
def build_graph():
    return build_multihop_graph(load_data(), verbose=False)

G = build_graph()
summary = get_multihop_summary(G)
n = summary["nodes_by_type"]
e = summary["edges_by_type"]


# ---------------------------------------------------------------------------
# Lookups
# ---------------------------------------------------------------------------
disease_map = {d["name"]: d["id"] for d in DISEASES}
target_map = {t["name"]: t["chembl_id"] for t in ANTICANCER_TARGETS}

approved_drugs = sorted(
    [(nid, d.get("name", nid)) for nid, d in G.nodes(data=True)
     if d.get("node_type") == "compound" and d.get("phase", 0) >= 4],
    key=lambda x: x[1],
)
drug_options = {name: nid for nid, name in approved_drugs}
# Add a sample of non-approved too
_all_compounds = sorted(
    [(nid, d.get("name", nid)) for nid, d in G.nodes(data=True)
     if d.get("node_type") == "compound"],
    key=lambda x: x[1],
)
for nid, name in _all_compounds[:30]:
    if name not in drug_options:
        drug_options[name] = nid


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def path_html(pr, top=False):
    """Render a path result as a styled card."""
    parts = []
    for nm, nt in zip(pr.path_names, pr.path_types):
        parts.append(f'<span class="tag tag-{nt}">{nm}</span>')
    chain = ' &rarr; '.join(parts)
    cls = "path-card path-card-top" if top else "path-card"
    dims = (f'P(eff)={pr.probs_by_dim["efficacy"]:.3f} &middot; '
            f'P(safe)={pr.probs_by_dim["safety"]:.3f} &middot; '
            f'P(evid)={pr.probs_by_dim["evidence"]:.3f}')
    return (f'<div class="{cls}">{chain}<br>'
            f'<span style="color:#6B7280;font-size:.8rem">'
            f'P = {pr.total_probability:.4f} &middot; {pr.n_hops} hops &middot; {dims}</span></div>')


def legend():
    return (
        '<div class="legend-row">'
        f'<span class="legend-item"><span class="legend-dot" style="background:{C["compound"]}"></span>Compound</span>'
        f'<span class="legend-item"><span class="legend-dot" style="background:{C["target"]}"></span>Target</span>'
        f'<span class="legend-item"><span class="legend-dot" style="background:{C["pathway"]}"></span>Pathway</span>'
        f'<span class="legend-item"><span class="legend-dot" style="background:{C["disease"]}"></span>Disease</span>'
        '</div>'
    )


def prob_color(p):
    if p > 0.5: return C["green"]
    if p > 0.2: return C["amber"]
    return C["red"]


# ---------------------------------------------------------------------------
# Sidebar — full navigation
# ---------------------------------------------------------------------------
with st.sidebar:
    st.markdown("## ChemPath")
    st.caption("Multi-Hop Drug Discovery  |  v2.1")
    st.markdown("---")

    page = st.radio(
        "Navigate",
        ["Drug Repurposing", "Mechanism of Action",
         "Target Identification", "Combination Therapy",
         "Network & Methodology"],
        label_visibility="collapsed",
    )

    st.markdown("---")

    # Optimization dimension — always visible
    st.markdown("##### Optimize for")
    weight_dim = st.radio(
        "dim", ["Efficacy", "Safety", "Evidence"],
        horizontal=True, label_visibility="collapsed",
    )
    weight_key = {"Efficacy": "w_efficacy", "Safety": "w_safety", "Evidence": "w_evidence"}[weight_dim]

    st.markdown("---")

    # Context-sensitive controls
    if page == "Drug Repurposing":
        st.markdown("##### Select Compound")
        sel_drug_name = st.selectbox(
            "drug", list(drug_options.keys()), label_visibility="collapsed",
            index=list(drug_options.keys()).index("GEFITINIB") if "GEFITINIB" in drug_options else 0,
        )
        sel_drug_id = drug_options[sel_drug_name]

    elif page == "Mechanism of Action":
        st.markdown("##### Compound")
        moa_drug_name = st.selectbox(
            "moa_drug", list(drug_options.keys()), label_visibility="collapsed",
            index=list(drug_options.keys()).index("IMATINIB") if "IMATINIB" in drug_options else 0,
        )
        moa_drug_id = drug_options[moa_drug_name]
        st.markdown("##### Disease")
        moa_disease_name = st.selectbox(
            "moa_disease", list(disease_map.keys()), label_visibility="collapsed",
            index=list(disease_map.keys()).index("CML") if "CML" in disease_map else 0,
        )
        moa_disease_id = disease_map[moa_disease_name]

    elif page == "Target Identification":
        st.markdown("##### Select Disease")
        tid_disease_name = st.selectbox(
            "tid_disease", list(disease_map.keys()), label_visibility="collapsed",
            index=list(disease_map.keys()).index("Melanoma") if "Melanoma" in disease_map else 0,
        )
        tid_disease_id = disease_map[tid_disease_name]

    elif page == "Combination Therapy":
        st.markdown("##### Select Targets")
        combo_targets = st.multiselect(
            "combo_targets", list(target_map.keys()), label_visibility="collapsed",
            default=["EGFR", "HER2"] if "EGFR" in target_map else list(target_map.keys())[:2],
        )
        combo_target_ids = [target_map[t] for t in combo_targets]

    st.markdown("---")
    st.markdown("##### Network")
    st.caption(
        f"{n['compound']:,} compounds  ·  {n['target']} targets\n\n"
        f"{n['pathway']} pathways  ·  {n['disease']} diseases\n\n"
        f"{summary['total_edges']:,} edges"
    )


# ============================================================================
# PAGE: Drug Repurposing
# ============================================================================
if page == "Drug Repurposing":
    st.markdown(f"# Drug Repurposing: **{sel_drug_name}**")
    st.markdown(
        '<p class="section-desc">'
        'Find new disease indications via multi-hop paths through the biological network. '
        'Each path represents a candidate biological mechanism.'
        '</p>', unsafe_allow_html=True)
    st.markdown(legend(), unsafe_allow_html=True)

    with st.spinner("Searching paths..."):
        results = find_repurposing_candidates(G, sel_drug_id, k_paths=5, weight_key=weight_key)

    if not results:
        st.warning("No paths found. This compound may lack target connectivity in the network.")
    else:
        # KPI row
        k1, k2, k3 = st.columns(3)
        k1.metric("Diseases Reachable", len(results))
        k2.metric("Best P(mechanism)", f"{results[0].best_probability:.4f}")
        k3.metric("Best Match", results[0].disease_name)

        # Bar chart + table side by side
        col_chart, col_table = st.columns([3, 2])

        with col_chart:
            diseases = [r.disease_name for r in results]
            probs = [r.best_probability for r in results]
            fig = go.Figure()
            fig.add_trace(go.Bar(
                y=diseases[::-1], x=probs[::-1], orientation="h",
                marker_color=[prob_color(p) for p in probs[::-1]],
                text=[f"{p:.3f}" for p in probs[::-1]], textposition="outside",
                hovertemplate="<b>%{y}</b><br>P = %{x:.4f}<extra></extra>",
            ))
            fig.update_layout(**_LAYOUT, height=max(280, len(diseases) * 38),
                              xaxis_title="P(drug reaches disease)",
                              xaxis_range=[0, max(probs) * 1.35] if probs else [0, 1])
            st.plotly_chart(fig, use_container_width=True)

        with col_table:
            st.markdown("#### Probability Breakdown")
            rows = []
            for r in results:
                bp = r.paths[0]
                rows.append({
                    "Disease": r.disease_name,
                    "P(mech)": f"{r.best_probability:.4f}",
                    "Hops": bp.n_hops,
                    "P(eff)": f"{bp.probs_by_dim['efficacy']:.3f}",
                    "P(safe)": f"{bp.probs_by_dim['safety']:.3f}",
                    "P(evid)": f"{bp.probs_by_dim['evidence']:.3f}",
                })
            st.dataframe(pd.DataFrame(rows), hide_index=True, use_container_width=True)

        # Detailed paths
        st.markdown("---")
        st.markdown("#### Candidate Mechanisms")
        for i, result in enumerate(results[:5]):
            with st.expander(
                f"{sel_drug_name} → {result.disease_name}  |  "
                f"P = {result.best_probability:.4f}  |  {len(result.paths)} path(s)",
                expanded=(i == 0),
            ):
                for j, p in enumerate(result.paths[:3]):
                    st.markdown(path_html(p, top=(j == 0)), unsafe_allow_html=True)

                # Edge evidence
                best = result.paths[0]
                for edge in best.edges_detail:
                    ev = edge.get("biological_evidence", "")
                    ic50 = edge.get("ic50_nm")
                    extra = f" (IC50 = {ic50:.1f} nM)" if ic50 else ""
                    if ev:
                        st.caption(f"**{edge['from']}** → **{edge['to']}**{extra}: {ev}")


# ============================================================================
# PAGE: Mechanism of Action
# ============================================================================
elif page == "Mechanism of Action":
    st.markdown(f"# Mechanism: **{moa_drug_name}** → **{moa_disease_name}**")
    st.markdown(
        '<p class="section-desc">'
        'Trace the most probable biological mechanism connecting a drug to a disease. '
        'Each path passes through specific targets and signaling pathways.'
        '</p>', unsafe_allow_html=True)
    st.markdown(legend(), unsafe_allow_html=True)

    with st.spinner("Finding mechanisms..."):
        moa = find_mechanism_of_action(G, moa_drug_id, moa_disease_id, k_paths=10)

    if moa is None or not moa.paths:
        st.warning(f"No path found from {moa_drug_name} to {moa_disease_name}.")
    else:
        pm = moa.primary_mechanism

        # Primary mechanism highlight
        st.markdown("#### Primary Mechanism")
        st.markdown(path_html(pm, top=True), unsafe_allow_html=True)

        # Sankey diagram
        if pm.n_hops >= 2:
            labels = pm.path_names
            node_colors = [NODE_COLORS.get(t, "#999") for t in pm.path_types]
            fig_s = go.Figure(go.Sankey(
                node=dict(pad=20, thickness=25, label=labels, color=node_colors,
                          line=dict(color="white", width=1.5)),
                link=dict(
                    source=list(range(len(labels) - 1)),
                    target=list(range(1, len(labels))),
                    value=[e["p_efficacy"] for e in pm.edges_detail],
                    label=[f'{e["edge_type"].replace("_"," ").title()}: P={e["p_efficacy"]:.3f}'
                           for e in pm.edges_detail],
                    color=[c + "40" for c in node_colors[:-1]],
                ),
            ))
            fig_s.update_layout(**_LAYOUT, height=260)
            st.plotly_chart(fig_s, use_container_width=True)

        # All candidate mechanisms
        if len(moa.paths) > 1:
            st.markdown("#### Alternative Mechanisms")
            for p in moa.paths[1:5]:
                st.markdown(path_html(p), unsafe_allow_html=True)

        # Evidence detail
        with st.expander("Biological evidence for each link"):
            for edge in pm.edges_detail:
                ev = edge.get("biological_evidence", "")
                ic50 = edge.get("ic50_nm")
                ic50_str = f"  |  IC50 = {ic50:.2f} nM" if ic50 else ""
                st.markdown(
                    f"**{edge['from']}** → **{edge['to']}** "
                    f"({edge['edge_type'].replace('_', ' ')}){ic50_str}"
                )
                if ev:
                    st.caption(f"_{ev}_")

        # Probability dimension comparison
        st.markdown("---")
        st.markdown("#### Probability Dimensions")
        dim_cols = st.columns(3)
        dim_cols[0].metric("P(efficacy)", f"{pm.probs_by_dim['efficacy']:.4f}")
        dim_cols[1].metric("P(safety)", f"{pm.probs_by_dim['safety']:.4f}")
        dim_cols[2].metric("P(evidence)", f"{pm.probs_by_dim['evidence']:.4f}")


# ============================================================================
# PAGE: Target Identification
# ============================================================================
elif page == "Target Identification":
    st.markdown(f"# Targets for: **{tid_disease_name}**")
    st.markdown(
        '<p class="section-desc">'
        'Discover which protein targets are most strongly connected to the disease '
        'by tracing backwards through pathways. Higher probability = stronger connection.'
        '</p>', unsafe_allow_html=True)
    st.markdown(legend(), unsafe_allow_html=True)

    with st.spinner("Identifying targets..."):
        tr = identify_targets_for_disease(G, tid_disease_id, weight_key=weight_key)

    if tr is None or not tr.targets:
        st.warning(f"No targets found for {tid_disease_name}.")
    else:
        # Horizontal bar chart
        tnames = [t["target_name"] for t in tr.targets]
        tprobs = [t["probability"] for t in tr.targets]

        fig_t = go.Figure()
        fig_t.add_trace(go.Bar(
            y=tnames[::-1], x=tprobs[::-1], orientation="h",
            marker_color=[prob_color(p) for p in tprobs[::-1]],
            text=[f"P={p:.3f}" for p in tprobs[::-1]], textposition="outside",
            hovertemplate="<b>%{y}</b><br>P = %{x:.4f}<extra></extra>",
        ))
        fig_t.update_layout(**_LAYOUT, height=max(300, len(tnames) * 45),
                            xaxis_title="P(target drives disease)",
                            xaxis_range=[0, max(tprobs) * 1.35] if tprobs else [0, 1])
        st.plotly_chart(fig_t, use_container_width=True)

        # Mechanism detail
        st.markdown("#### Mechanisms")
        for t in tr.targets:
            bp = t["best_path"]
            st.markdown(path_html(bp, top=(t == tr.targets[0])), unsafe_allow_html=True)

        # Radar: compare top targets across dims
        if len(tr.targets) >= 2:
            st.markdown("---")
            st.markdown("#### Probability Dimension Comparison")
            cats = ["Efficacy", "Safety", "Evidence"]
            fig_r = go.Figure()
            for t in tr.targets[:6]:
                bp = t["best_path"]
                v = [bp.probs_by_dim["efficacy"], bp.probs_by_dim["safety"], bp.probs_by_dim["evidence"]]
                fig_r.add_trace(go.Scatterpolar(
                    r=v + [v[0]], theta=cats + [cats[0]],
                    fill="toself", opacity=0.2, name=t["target_name"],
                ))
            fig_r.update_layout(**_LAYOUT, height=380,
                                polar=dict(radialaxis=dict(visible=True, range=[0, 1])),
                                legend=dict(orientation="h", y=-0.1, xanchor="center", x=0.5))
            st.plotly_chart(fig_r, use_container_width=True)


# ============================================================================
# PAGE: Combination Therapy
# ============================================================================
elif page == "Combination Therapy":
    st.markdown("# Combination Therapy")
    st.markdown(
        '<p class="section-desc">'
        'Find the minimum set of compounds covering multiple targets simultaneously. '
        'This is a set cover problem — a genuine graph optimization.'
        '</p>', unsafe_allow_html=True)

    if len(combo_targets) < 2:
        st.info("Select at least 2 targets in the sidebar.")
    else:
        st.markdown(f"**Targets:** {', '.join(combo_targets)}")

        with st.spinner("Computing combinations..."):
            combo = find_combination_therapy(G, combo_target_ids, weight_key=weight_key)

        if combo is None or not combo.best_combination:
            st.warning("No compounds found covering the selected targets.")
        else:
            # Best combination
            st.markdown(f"#### Optimal: {len(combo.best_combination)} compound(s)")
            for i, comp in enumerate(combo.best_combination):
                covered = ", ".join(comp["target_names_covered"])
                st.markdown(
                    f'<div class="path-card path-card-top">'
                    f'<span class="tag tag-compound">{comp["compound_name"]}</span> '
                    f'covers <strong>{covered}</strong> '
                    f'(avg cost = {comp["avg_cost"]:.3f})</div>',
                    unsafe_allow_html=True,
                )

            # Coverage matrix
            st.markdown("#### Top Candidates")
            matrix = []
            for comp in combo.compounds[:20]:
                row = {"Compound": comp["compound_name"], "Covered": comp["n_covered"]}
                for tname in combo_targets:
                    tid = target_map[tname]
                    row[tname] = "Yes" if tid in comp["targets_covered"] else "-"
                row["Avg Cost"] = round(comp["avg_cost"], 3)
                matrix.append(row)
            df_m = pd.DataFrame(matrix)

            def _hl(val):
                if val == "Yes":
                    return f"background-color:{C['green_bg']};color:{C['green']};font-weight:600"
                return ""

            st.dataframe(
                df_m.style.map(_hl, subset=combo_targets).format({"Avg Cost": "{:.3f}"}),
                hide_index=True, use_container_width=True,
            )

            # Coverage distribution
            cov_counts = {}
            for comp in combo.compounds:
                nc = comp["n_covered"]
                cov_counts[nc] = cov_counts.get(nc, 0) + 1
            xs = sorted(cov_counts.keys())
            fig_c = go.Figure()
            fig_c.add_trace(go.Bar(
                x=[f"{x} target{'s' if x > 1 else ''}" for x in xs],
                y=[cov_counts[x] for x in xs],
                marker_color=[C["green"] if x == len(combo_targets) else C["blue"] for x in xs],
                text=[cov_counts[x] for x in xs], textposition="outside",
            ))
            fig_c.update_layout(**_LAYOUT, height=320,
                                xaxis_title="Targets Covered", yaxis_title="Compounds")
            st.plotly_chart(fig_c, use_container_width=True)


# ============================================================================
# PAGE: Network & Methodology
# ============================================================================
elif page == "Network & Methodology":
    st.markdown("# Network & Methodology")

    # --- Network Backbone Visualization ---
    st.markdown("#### Signaling Network Backbone")
    st.markdown(
        '<p class="section-desc">'
        'Target → Pathway → Disease connections that form multi-hop routes. '
        'Edge thickness reflects probability strength. Dashed lines = PPI cross-links.'
        '</p>', unsafe_allow_html=True)
    st.markdown(legend(), unsafe_allow_html=True)

    targets_n = [(nid, d) for nid, d in G.nodes(data=True) if d.get("node_type") == "target"]
    pathways_n = [(nid, d) for nid, d in G.nodes(data=True) if d.get("node_type") == "pathway"]
    diseases_n = [(nid, d) for nid, d in G.nodes(data=True) if d.get("node_type") == "disease"]

    pos = {}
    for i, (nid, _) in enumerate(targets_n):
        pos[nid] = (0, -i * 1.2)
    for i, (nid, _) in enumerate(pathways_n):
        pos[nid] = (2, -i * 1.2)
    for i, (nid, _) in enumerate(diseases_n):
        pos[nid] = (4, -i * 1.2)

    fig_net = go.Figure()
    for u, v, d in G.edges(data=True):
        etype = d.get("edge_type", "")
        if etype not in ("target_pathway", "pathway_disease", "ppi"):
            continue
        if u not in pos or v not in pos:
            continue
        p_eff = d.get("p_efficacy", 0.5)
        fig_net.add_trace(go.Scatter(
            x=[pos[u][0], pos[v][0]], y=[pos[u][1], pos[v][1]],
            mode="lines",
            line=dict(width=max(0.5, p_eff * 4),
                      color="#C4B5FD" if etype == "ppi" else "#D0D5DD",
                      dash="dot" if etype == "ppi" else "solid"),
            hovertext=(f"{G.nodes[u].get('name', u)} → {G.nodes[v].get('name', v)}<br>"
                       f"P(eff)={d.get('p_efficacy', 0):.2f}"),
            hoverinfo="text", showlegend=False,
        ))

    for lname, lnodes, color in [
        ("Targets", targets_n, C["target"]),
        ("Pathways", pathways_n, C["pathway"]),
        ("Diseases", diseases_n, C["disease"]),
    ]:
        xs = [pos[nid][0] for nid, _ in lnodes if nid in pos]
        ys = [pos[nid][1] for nid, _ in lnodes if nid in pos]
        names = [d.get("name", nid) for nid, d in lnodes if nid in pos]
        tp = "middle right" if lname != "Diseases" else "middle left"
        fig_net.add_trace(go.Scatter(
            x=xs, y=ys, mode="markers+text",
            marker=dict(size=20, color=color, line=dict(width=2, color="white")),
            text=names, textposition=tp, textfont=dict(size=10, color=color),
            name=lname, hovertemplate="<b>%{text}</b><extra>" + lname + "</extra>",
        ))

    fig_net.update_layout(
        **_LAYOUT, height=600, showlegend=True,
        legend=dict(orientation="h", y=1.02, xanchor="center", x=0.5),
        xaxis=dict(visible=False), yaxis=dict(visible=False),
    )
    st.plotly_chart(fig_net, use_container_width=True)

    # --- Methodology ---
    st.markdown("---")
    st.markdown("#### Probability Framework")
    st.markdown("""
Every edge carries a **3-dimensional probability vector**:

| Dimension | Meaning | Transform |
|-----------|---------|-----------|
| P(efficacy) | Therapeutic contribution | IC50 → sigmoid; pathway centrality |
| P(safety) | No adverse effects | 1 - toxicity; PPI type penalty |
| P(evidence) | Link is experimentally real | Clinical phase; data source |

In log-space: `w = -log(P)`, so **shortest path = maximum probability path**.

Each dimension is independently additive along a path:
`P(path) = P(edge1) × P(edge2) × ... × P(edgeN)`

Multi-objective Pareto optimization finds paths that are optimal across all 3 dimensions simultaneously.
""")

    st.markdown("#### Graph Structure")
    st.markdown(f"""
```
Layer 0: Compounds ({n['compound']:,} from ChEMBL)
  ↓  P(binding) from IC50 sigmoid
Layer 1: Targets ({n['target']} anti-cancer kinases)
  ↓  P(membership) from KEGG/Reactome
  ↔  PPI: P(interaction) from STRING
Layer 2: Pathways ({n['pathway']} signaling cascades)
  ↓  P(association) from DisGeNET/KEGG
Layer 3: Diseases ({n['disease']} cancer types)
```
""")

    st.markdown("#### Validation")
    st.markdown("""
| Known Pair | Path Found | P | Correct? |
|-----------|-----------|---|----------|
| Imatinib → CML | Imatinib → ABL1 → BCR-ABL → CML | 0.587 | Yes |
| Gefitinib → NSCLC | Gefitinib → EGFR → ErbB → NSCLC | 0.485 | Yes |
| BRAF → Melanoma | BRAF → MAPK → Melanoma | 0.640 | Yes |
| Gefitinib → HCC | Gefitinib → EGFR → MET → HGF → HCC | 0.277 | Repurposing hypothesis |
""")

    st.markdown("---")
    st.caption("ChemPath v2.1 — Data: ChEMBL | Network: KEGG, STRING, DisGeNET (curated)")
