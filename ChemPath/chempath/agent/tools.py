"""
Claude Tool Use definitions for ChemPath.

Each tool maps a structured function call to our graph optimization pipeline.
Tools are defined in Anthropic's tool_use format (MCP-compatible).
"""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from chempath.data.chembl_client import load_saved_data, ANTICANCER_TARGETS
from chempath.data.loader import get_target_by_name
from chempath.graph.builder import build_drug_target_graph, get_graph_summary
from chempath.graph.optimizer import (
    rank_compounds_for_target,
    compute_pareto_front,
    find_knee_point,
    compute_selectivity,
    compare_compounds_across_targets,
    DrugRecommendation,
)
from chempath.chemistry.properties import compute_properties
from chempath.graph.analysis import (
    toxicity_penalty_sensitivity,
    ic50_perturbation_sensitivity,
    missing_data_sensitivity,
)

DATA_PATH = Path(__file__).parent.parent.parent / "data" / "chembl_real.json"

# --- Tool Definitions (Anthropic format) ---

TOOL_DEFINITIONS = [
    {
        "name": "list_targets",
        "description": (
            "List all available drug targets in the database. "
            "Use this when the user asks what targets are available or wants to explore the data."
        ),
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": [],
        },
    },
    {
        "name": "screen_compounds",
        "description": (
            "Screen and rank compounds for a specific protein target. "
            "Returns ranked compounds with IC50, toxicity, confidence labels, and optimization weights. "
            "Use this when the user asks to find drugs for a target, rank inhibitors, or screen compounds."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {
                    "type": "string",
                    "description": "Target protein name (e.g., 'EGFR', 'BRAF', 'ALK') or ChEMBL ID (e.g., 'CHEMBL203')",
                },
                "strategy": {
                    "type": "string",
                    "enum": ["balanced", "efficacy", "safety"],
                    "description": "Ranking strategy: 'balanced' (default, considers both IC50 and toxicity), 'efficacy' (IC50 only), 'safety' (toxicity first)",
                },
                "top_n": {
                    "type": "integer",
                    "description": "Number of top compounds to return (default: 10)",
                },
                "max_ic50": {
                    "type": "number",
                    "description": "Maximum IC50 in nM to filter results (optional)",
                },
            },
            "required": ["target"],
        },
    },
    {
        "name": "compute_pareto",
        "description": (
            "Compute the Pareto-optimal set of compounds for a target, balancing efficacy (IC50) and safety (toxicity). "
            "Also identifies the 'knee point' — the compound with the best balance. "
            "Use this when the user asks about trade-offs, Pareto front, or multi-objective optimization."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {
                    "type": "string",
                    "description": "Target protein name or ChEMBL ID",
                },
            },
            "required": ["target"],
        },
    },
    {
        "name": "run_sensitivity",
        "description": (
            "Run sensitivity analysis to test how robust the compound rankings are. "
            "Tests three dimensions: toxicity penalty weight, IC50 measurement noise, and missing data. "
            "Use this when the user asks about robustness, reliability, or sensitivity of recommendations."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {
                    "type": "string",
                    "description": "Target protein name or ChEMBL ID",
                },
                "analysis_type": {
                    "type": "string",
                    "enum": ["all", "toxicity_penalty", "ic50_noise", "missing_data"],
                    "description": "Which sensitivity analysis to run (default: 'all')",
                },
            },
            "required": ["target"],
        },
    },
    {
        "name": "get_compound_info",
        "description": (
            "Get detailed information about a specific compound, including all targets it binds to. "
            "Use this when the user asks about a specific drug or compound."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "compound": {
                    "type": "string",
                    "description": "Compound name (e.g., 'Sirolimus') or ChEMBL ID",
                },
            },
            "required": ["compound"],
        },
    },
    {
        "name": "compare_selectivity",
        "description": (
            "Analyze selectivity of a compound for its primary target vs off-targets. "
            "Shows selectivity ratios (higher = more selective). "
            "Use this when the user asks about selectivity, off-target effects, or specificity."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "compound": {
                    "type": "string",
                    "description": "Compound name or ChEMBL ID",
                },
                "primary_target": {
                    "type": "string",
                    "description": "The intended target name or ChEMBL ID",
                },
            },
            "required": ["compound", "primary_target"],
        },
    },
    {
        "name": "head_to_head",
        "description": (
            "Compare multiple compounds across multiple targets side-by-side. "
            "Shows an IC50 matrix and molecular properties. "
            "Use this when the user wants to compare drugs or asks 'which is better'."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "target": {
                    "type": "string",
                    "description": "Primary target to find top compounds for",
                },
                "top_n": {
                    "type": "integer",
                    "description": "Number of top compounds to compare (default: 5)",
                },
            },
            "required": ["target"],
        },
    },
    {
        "name": "graph_summary",
        "description": (
            "Get a summary of the drug-target interaction graph: number of compounds, targets, and edges. "
            "Use this when the user asks about the database size or graph structure."
        ),
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": [],
        },
    },
]


# --- Tool Execution ---

class ChemPathToolExecutor:
    """Executes tool calls against the ChemPath pipeline."""

    def __init__(self, data_path: str | Path = DATA_PATH):
        self.data_path = Path(data_path)
        self._data = None
        self._graph = None

    @property
    def data(self):
        if self._data is None:
            self._data = load_saved_data(self.data_path)
        return self._data

    @property
    def graph(self):
        if self._graph is None:
            self._graph = build_drug_target_graph(self.data, toxicity_penalty=1.0, verbose=False)
        return self._graph

    def _resolve_target(self, target_str: str) -> tuple[str, str]:
        """Resolve a target name or ChEMBL ID to (chembl_id, name)."""
        # Direct ChEMBL ID
        if target_str.startswith("CHEMBL"):
            for t in self.data["targets"]:
                if t["chembl_id"] == target_str:
                    return t["chembl_id"], t["name"]
            return target_str, target_str

        # Name lookup (case-insensitive, partial match)
        target_lower = target_str.lower()
        for t in self.data["targets"]:
            if t["name"].lower() == target_lower:
                return t["chembl_id"], t["name"]
        # Partial match
        for t in self.data["targets"]:
            if target_lower in t["name"].lower() or target_lower in t.get("description", "").lower():
                return t["chembl_id"], t["name"]
        # Also check ANTICANCER_TARGETS for descriptions
        for t in ANTICANCER_TARGETS:
            if target_lower in t["name"].lower() or target_lower in t["description"].lower():
                return t["chembl_id"], t["name"]

        return None, target_str

    def _rec_to_dict(self, rec: DrugRecommendation) -> dict:
        """Convert a DrugRecommendation to a clean dict for LLM consumption."""
        return {
            "rank": rec.rank,
            "compound_id": rec.compound_id,
            "compound_name": rec.compound_name,
            "ic50_nm": rec.ic50_nm,
            "toxicity": rec.toxicity,
            "weight": round(rec.weight, 4),
            "source": rec.source,
            "confidence": rec.confidence,
            "confidence_reasons": rec.confidence_reasons,
        }

    def execute(self, tool_name: str, tool_input: dict) -> str:
        """Execute a tool call and return a string result."""
        dispatch = {
            "list_targets": self._list_targets,
            "screen_compounds": self._screen_compounds,
            "compute_pareto": self._compute_pareto,
            "run_sensitivity": self._run_sensitivity,
            "get_compound_info": self._get_compound_info,
            "compare_selectivity": self._compare_selectivity,
            "head_to_head": self._head_to_head,
            "graph_summary": self._graph_summary,
        }
        handler = dispatch.get(tool_name)
        if handler is None:
            return f"Unknown tool: {tool_name}"
        return handler(tool_input)

    def _list_targets(self, _input: dict) -> str:
        targets = self.data["targets"]
        lines = [f"Available targets ({len(targets)}):"]
        for t in targets:
            count = sum(1 for a in self.data["bioactivities"] if a["target"] == t["chembl_id"])
            lines.append(f"  - {t['name']} ({t['chembl_id']}): {count} compounds tested")
        return "\n".join(lines)

    def _screen_compounds(self, input: dict) -> str:
        target_str = input["target"]
        strategy = input.get("strategy", "balanced")
        top_n = input.get("top_n", 10)
        max_ic50 = input.get("max_ic50")

        target_id, target_name = self._resolve_target(target_str)
        if target_id is None:
            return f"Target '{target_str}' not found. Use list_targets to see available targets."

        recs = rank_compounds_for_target(
            self.graph, target_id, strategy=strategy
        )
        if not recs:
            return f"No compounds found for target {target_name} ({target_id})."

        if max_ic50 is not None:
            recs = [r for r in recs if r.ic50_nm <= max_ic50]
            if not recs:
                return f"No compounds found with IC50 <= {max_ic50}nM for {target_name}."

        top = recs[:top_n]
        results = [self._rec_to_dict(r) for r in top]

        lines = [
            f"Top {len(top)} compounds for {target_name} ({target_id}), strategy={strategy}:",
            f"Total candidates: {len(recs)}",
            "",
        ]
        for r in results:
            lines.append(
                f"  #{r['rank']} [{r['confidence']}] {r['compound_name']} — "
                f"IC50={r['ic50_nm']:.2f}nM, Tox={r['toxicity']:.2f}, "
                f"Weight={r['weight']:.4f}, Source={r['source']}"
            )
            for reason in r["confidence_reasons"]:
                lines.append(f"      - {reason}")

        return "\n".join(lines)

    def _compute_pareto(self, input: dict) -> str:
        target_str = input["target"]
        target_id, target_name = self._resolve_target(target_str)
        if target_id is None:
            return f"Target '{target_str}' not found."

        recs = rank_compounds_for_target(
            self.graph, target_id, self.data.get("toxicity"), strategy="balanced"
        )
        if not recs:
            return f"No compounds found for {target_name}."

        pareto = compute_pareto_front(recs)
        knee = find_knee_point(pareto)

        lines = [
            f"Pareto-optimal compounds for {target_name} (efficacy vs toxicity):",
            f"  {len(pareto)} Pareto-optimal from {len(recs)} total candidates",
            "",
        ]
        for r in pareto[:10]:
            marker = " <<< KNEE POINT" if knee and r.compound_id == knee.compound_id else ""
            lines.append(
                f"  #{r.rank} [{r.confidence}] {r.compound_name} — "
                f"IC50={r.ic50_nm:.2f}nM, Tox={r.toxicity:.2f}{marker}"
            )

        if knee:
            lines.append(f"\nRecommended (knee point): {knee.compound_name}")
            lines.append(f"  Best balance of efficacy (IC50={knee.ic50_nm:.2f}nM) and safety (Tox={knee.toxicity:.2f})")

        return "\n".join(lines)

    def _run_sensitivity(self, input: dict) -> str:
        target_str = input["target"]
        analysis_type = input.get("analysis_type", "all")
        target_id, target_name = self._resolve_target(target_str)
        if target_id is None:
            return f"Target '{target_str}' not found."

        lines = [f"Sensitivity analysis for {target_name}:"]
        analyses = []

        if analysis_type in ("all", "toxicity_penalty"):
            result = toxicity_penalty_sensitivity(self.data, target_id)
            analyses.append(("Toxicity Penalty Weight", result))

        if analysis_type in ("all", "ic50_noise"):
            result = ic50_perturbation_sensitivity(self.data, target_id, n_trials=3)
            analyses.append(("IC50 Measurement Noise", result))

        if analysis_type in ("all", "missing_data"):
            result = missing_data_sensitivity(self.data, target_id, n_trials=3)
            analyses.append(("Missing Data", result))

        for name, result in analyses:
            lines.append(f"\n  {name} ({result.parameter}):")
            for val in result.values_tested:
                ranking = result.rankings_per_value.get(val, [])
                lines.append(f"    {val}: {' > '.join(ranking) if ranking else '(no data)'}")
            lines.append(f"    Stable compounds: {', '.join(result.stable_compounds) or 'None'}")
            lines.append(f"    Volatile compounds: {', '.join(result.volatile_compounds) or 'None'}")

        return "\n".join(lines)

    def _get_compound_info(self, input: dict) -> str:
        compound_str = input["compound"]

        # Search by name or ID
        found = None
        for c in self.data["compounds"]:
            if (c["chembl_id"].lower() == compound_str.lower() or
                    c["name"].lower() == compound_str.lower()):
                found = c
                break
        # Partial name match
        if found is None:
            for c in self.data["compounds"]:
                if compound_str.lower() in c["name"].lower():
                    found = c
                    break

        if found is None:
            return f"Compound '{compound_str}' not found in the database."

        cid = found["chembl_id"]
        activities = [a for a in self.data["bioactivities"] if a["compound"] == cid]

        lines = [
            f"Compound: {found['name']} ({cid})",
            f"  SMILES: {found['smiles'][:80]}{'...' if len(found['smiles']) > 80 else ''}",
            f"  Clinical Phase: {found.get('phase', 'Unknown')}",
            f"  Tested against {len(activities)} target(s):",
        ]
        for a in activities:
            target_name = a.get("target", "")
            for t in self.data["targets"]:
                if t["chembl_id"] == a["target"]:
                    target_name = f"{t['name']} ({a['target']})"
                    break
            lines.append(f"    - {target_name}: IC50={a['value']}nM")

        return "\n".join(lines)

    def _compare_selectivity(self, input: dict) -> str:
        compound_str = input["compound"]
        primary_str = input["primary_target"]

        # Resolve compound
        found = None
        for c in self.data["compounds"]:
            if (c["chembl_id"].lower() == compound_str.lower() or
                    c["name"].lower() == compound_str.lower()):
                found = c
                break
        if found is None:
            for c in self.data["compounds"]:
                if compound_str.lower() in c["name"].lower():
                    found = c
                    break
        if found is None:
            return f"Compound '{compound_str}' not found."

        target_id, target_name = self._resolve_target(primary_str)
        if target_id is None:
            return f"Target '{primary_str}' not found."

        result = compute_selectivity(self.graph, found["chembl_id"], target_id)
        if "error" in result:
            return result["error"]

        lines = [
            f"Selectivity profile for {result['compound_name']}:",
            f"  Primary target: {target_name} — IC50={result['primary_ic50_nm']:.2f}nM",
            f"  Off-targets tested: {result['off_target_count']}",
            f"  Selective (>10x ratio for all off-targets): {'Yes' if result['is_selective'] else 'No'}",
        ]
        if result["off_targets"]:
            lines.append("\n  Off-target activity:")
            for ot in result["off_targets"]:
                selectivity = f"{ot['selectivity_ratio']:.1f}x"
                risk = "LOW RISK" if ot["selectivity_ratio"] > 100 else (
                    "MODERATE" if ot["selectivity_ratio"] > 10 else "HIGH RISK"
                )
                lines.append(
                    f"    - {ot['target_name']}: IC50={ot['ic50_nm']:.2f}nM, "
                    f"Ratio={selectivity} [{risk}]"
                )
        else:
            lines.append("  No off-target data available (tested against only 1 target)")

        # Add molecular properties
        props = compute_properties(found["smiles"])
        lines.append(f"\n  Molecular properties:")
        lines.append(f"    Est. MW: {props.estimated_mw:.0f}")
        lines.append(f"    Lipinski violations: {props.lipinski_violations}/4")
        lines.append(f"    Risk score: {props.risk_score:.2f}")
        if props.risk_flags:
            for flag in props.risk_flags:
                lines.append(f"    Warning: {flag}")

        return "\n".join(lines)

    def _head_to_head(self, input: dict) -> str:
        target_str = input["target"]
        top_n = input.get("top_n", 5)

        target_id, target_name = self._resolve_target(target_str)
        if target_id is None:
            return f"Target '{target_str}' not found."

        recs = rank_compounds_for_target(
            self.graph, target_id, self.data.get("toxicity"), strategy="balanced"
        )
        if not recs:
            return f"No compounds for {target_name}."

        top = recs[:top_n]
        compound_ids = [r.compound_id for r in top]
        target_ids = [t["chembl_id"] for t in self.data["targets"]]

        matrix = compare_compounds_across_targets(self.graph, compound_ids, target_ids)

        lines = [
            f"Head-to-head comparison: Top {len(top)} {target_name} inhibitors",
            f"{'=' * 80}",
        ]

        # Header row
        target_names = [t["name"] for t in self.data["targets"]]
        header = f"  {'Compound':20s}"
        for tn in target_names:
            header += f" {tn:>8s}"
        header += "  Risk  Conf"
        lines.append(header)
        lines.append("  " + "-" * (len(header) - 2))

        for row, rec in zip(matrix, top):
            props = compute_properties(
                next(c["smiles"] for c in self.data["compounds"] if c["chembl_id"] == row["compound_id"])
            )
            line = f"  {row['compound_name'][:20]:20s}"
            for tn in target_names:
                val = row["targets"].get(tn)
                if val is not None:
                    if val < 1:
                        line += f" {val:>7.3f}*"
                    elif val < 100:
                        line += f" {val:>8.1f}"
                    else:
                        line += f" {val:>8.0f}"
                else:
                    line += f" {'—':>8s}"
            line += f"  {props.risk_score:.2f}  {rec.confidence}"
            lines.append(line)

        lines.append(f"{'=' * 80}")
        lines.append("  Values are IC50 in nM (lower = more potent). * = sub-nanomolar")
        lines.append("  Risk = molecular property risk score (0=safe, 1=risky)")
        lines.append(f"  '—' = not tested against that target")

        return "\n".join(lines)

    def _graph_summary(self, _input: dict) -> str:
        summary = get_graph_summary(self.graph)

        # Per-target compound counts
        target_stats = []
        for t in self.data["targets"]:
            count = sum(1 for a in self.data["bioactivities"] if a["target"] == t["chembl_id"])
            target_stats.append(f"    {t['name']}: {count}")

        # Clinical phase stats
        phases = {}
        for c in self.data["compounds"]:
            p = c.get("phase", 0)
            phases[p] = phases.get(p, 0) + 1
        phase_str = ", ".join(f"Phase {k}: {v}" for k, v in sorted(phases.items()) if k > 0)

        # IC50 range
        ic50s = [a["value"] for a in self.data["bioactivities"] if a.get("value", 0) > 0]
        ic50_min = min(ic50s) if ic50s else 0
        ic50_max = max(ic50s) if ic50s else 0

        lines = [
            f"Drug-Target Interaction Graph:",
            f"  Compounds: {summary['compounds']}",
            f"  Targets: {summary['targets']}",
            f"  Total edges: {summary['total_edges']} ({summary['experimental_edges']} experimental, {summary['predicted_edges']} predicted)",
            f"  IC50 range: {ic50_min:.3f} — {ic50_max:.0f} nM",
            f"  Clinical phase data: {phase_str or 'None'}",
            f"  Compounds per target:",
        ]
        lines.extend(target_stats)
        lines.append(f"  Data source: ChEMBL REST API")

        return "\n".join(lines)
