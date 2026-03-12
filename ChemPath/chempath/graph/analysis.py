"""
Sensitivity analysis — demonstrates awareness of model limitations.

Three analyses:
1. Toxicity penalty sensitivity — ranking stability vs penalty weight
2. IC50 perturbation — ranking stability vs measurement noise
3. Missing data impact — ranking stability vs data coverage
"""

import copy
import random
from collections import Counter
from dataclasses import dataclass

from chempath.graph.builder import build_drug_target_graph
from chempath.graph.optimizer import rank_compounds_for_target


@dataclass
class SensitivityResult:
    parameter: str
    values_tested: list[float]
    rankings_per_value: dict
    stable_compounds: list[str]
    volatile_compounds: list[str]


def toxicity_penalty_sensitivity(
    data: dict, target_id: str, penalties: list[float] | None = None, top_n: int = 3
) -> SensitivityResult:
    """Test how rankings change as toxicity penalty weight varies."""
    if penalties is None:
        penalties = [0.0, 0.5, 1.0, 2.0, 5.0]

    rankings = {}
    all_top_n = []

    for penalty in penalties:
        G = build_drug_target_graph(data, toxicity_penalty=penalty, verbose=False)
        recs = rank_compounds_for_target(G, target_id, data.get("toxicity"), strategy="balanced")
        top = [r.compound_name for r in recs[:top_n]]
        rankings[penalty] = top
        all_top_n.append(set(top))

    stable = list(set.intersection(*all_top_n)) if all_top_n else []
    all_seen = set.union(*all_top_n) if all_top_n else set()
    volatile = [c for c in all_seen if c not in stable]

    return SensitivityResult(
        parameter="toxicity_penalty",
        values_tested=penalties,
        rankings_per_value=rankings,
        stable_compounds=stable,
        volatile_compounds=volatile,
    )


def ic50_perturbation_sensitivity(
    data: dict, target_id: str,
    noise_levels: list[float] | None = None, n_trials: int = 5, top_n: int = 3,
) -> SensitivityResult:
    """Test ranking stability under IC50 measurement noise (±X%)."""
    if noise_levels is None:
        noise_levels = [0.0, 0.1, 0.2, 0.5, 1.0]

    rankings = {}
    all_top_n = []

    for noise in noise_levels:
        trial_tops = []
        for _ in range(n_trials):
            perturbed_data = copy.deepcopy(data)
            for activity in perturbed_data["bioactivities"]:
                if noise > 0:
                    factor = 1.0 + random.uniform(-noise, noise)
                    activity["value"] = max(0.01, activity["value"] * factor)

            G = build_drug_target_graph(perturbed_data, verbose=False)
            recs = rank_compounds_for_target(G, target_id, data.get("toxicity"), strategy="balanced")
            trial_tops.append(tuple(r.compound_name for r in recs[:top_n]))

        most_common = Counter(trial_tops).most_common(1)[0][0]
        rankings[noise] = list(most_common)
        all_top_n.append(set(most_common))

    stable = list(set.intersection(*all_top_n)) if all_top_n else []
    all_seen = set.union(*all_top_n) if all_top_n else set()
    volatile = [c for c in all_seen if c not in stable]

    return SensitivityResult(
        parameter="ic50_noise_fraction",
        values_tested=noise_levels,
        rankings_per_value=rankings,
        stable_compounds=stable,
        volatile_compounds=volatile,
    )


def missing_data_sensitivity(
    data: dict, target_id: str,
    removal_fractions: list[float] | None = None, n_trials: int = 5, top_n: int = 3,
) -> SensitivityResult:
    """Test what happens when bioactivity edges are randomly removed."""
    if removal_fractions is None:
        removal_fractions = [0.0, 0.1, 0.2, 0.3, 0.5]

    rankings = {}
    all_top_n = []

    for frac in removal_fractions:
        trial_tops = []
        for _ in range(n_trials):
            perturbed_data = copy.deepcopy(data)
            activities = perturbed_data["bioactivities"]
            n_remove = int(len(activities) * frac)
            if n_remove > 0:
                to_remove = set(random.sample(range(len(activities)), min(n_remove, len(activities))))
                perturbed_data["bioactivities"] = [a for i, a in enumerate(activities) if i not in to_remove]

            G = build_drug_target_graph(perturbed_data, verbose=False)
            recs = rank_compounds_for_target(G, target_id, data.get("toxicity"), strategy="balanced")
            trial_tops.append(tuple(r.compound_name for r in recs[:top_n]))

        if trial_tops:
            most_common = Counter(trial_tops).most_common(1)[0][0]
            rankings[frac] = list(most_common)
            all_top_n.append(set(most_common))

    stable = list(set.intersection(*all_top_n)) if all_top_n else []
    all_seen = set.union(*all_top_n) if all_top_n else set()
    volatile = [c for c in all_seen if c not in stable]

    return SensitivityResult(
        parameter="data_removal_fraction",
        values_tested=removal_fractions,
        rankings_per_value=rankings,
        stable_compounds=stable,
        volatile_compounds=volatile,
    )


def format_sensitivity_report(result: SensitivityResult) -> str:
    """Format a sensitivity analysis result as a readable report."""
    lines = [
        f"\n{'=' * 70}",
        f" Sensitivity Analysis: {result.parameter}",
        f"{'=' * 70}",
    ]

    for val in result.values_tested:
        ranking = result.rankings_per_value.get(val, [])
        ranking_str = " > ".join(ranking) if ranking else "(no data)"
        lines.append(f"  {result.parameter}={val:<5} : {ranking_str}")

    lines.append(f"\n  Stable (robust) compounds: {', '.join(result.stable_compounds) or 'None'}")
    lines.append(f"  Volatile (sensitive) compounds: {', '.join(result.volatile_compounds) or 'None'}")

    if result.stable_compounds:
        lines.append(f"\n  INTERPRETATION: {', '.join(result.stable_compounds)} "
                      f"consistently rank in the top across all tested {result.parameter} values.")
        lines.append(f"  These recommendations are ROBUST to parameter changes.")
    if result.volatile_compounds:
        lines.append(f"  {', '.join(result.volatile_compounds)} are SENSITIVE to {result.parameter}.")
        lines.append(f"  Recommendations involving these compounds should be treated with CAUTION.")

    lines.append(f"{'=' * 70}")
    return "\n".join(lines)


def run_full_sensitivity_analysis(data: dict, target_id: str) -> list[SensitivityResult]:
    """Run all three sensitivity analyses."""
    results = []

    print("\n  Running toxicity penalty sensitivity...")
    results.append(toxicity_penalty_sensitivity(data, target_id))

    print("  Running IC50 perturbation sensitivity...")
    results.append(ic50_perturbation_sensitivity(data, target_id))

    print("  Running missing data sensitivity...")
    results.append(missing_data_sensitivity(data, target_id))

    return results
