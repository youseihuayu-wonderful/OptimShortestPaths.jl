"""
Fetch real IC50 data from ChEMBL for anti-cancer targets.
Saves result to data/chembl_real.json for offline use.
"""

from pathlib import Path
from chempath.data.chembl_client import (
    fetch_real_data, enrich_compounds_with_phase, save_data, ANTICANCER_TARGETS,
)

DATA_DIR = Path(__file__).parent.parent / "data"


def main():
    print("=" * 70)
    print(" ChemPath — Fetch Real ChEMBL Data")
    print("=" * 70)

    data = fetch_real_data(
        targets=ANTICANCER_TARGETS,
        max_per_target=500,
        use_cache=True,
        verbose=True,
    )

    # Enrich with clinical phase
    print("\nEnriching compounds with clinical phase data...")
    enrich_compounds_with_phase(data, use_cache=True, verbose=True)

    out_path = DATA_DIR / "chembl_real.json"
    save_data(data, out_path)

    # Print summary per target
    print("\nPer-target breakdown:")
    for t in ANTICANCER_TARGETS:
        tid = t["chembl_id"]
        count = sum(1 for a in data["bioactivities"] if a["target"] == tid)
        print(f"  {t['name']:8s} ({tid}): {count} activities")

    print("\nDone. Run 'uv run python scripts/demo_real.py' to test with real data.")


if __name__ == "__main__":
    main()
