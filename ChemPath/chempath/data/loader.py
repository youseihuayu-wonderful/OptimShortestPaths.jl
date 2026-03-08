"""
Unified data loading interface.
Currently wraps mock data; will add ChEMBL REST API client later.
"""


def get_compound_by_name(data: dict, name: str) -> dict | None:
    """Look up a compound by name (case-insensitive)."""
    for c in data["compounds"]:
        if c["name"].lower() == name.lower():
            return c
    return None


def get_target_by_name(data: dict, name: str) -> dict | None:
    """Look up a target by name (case-insensitive)."""
    for t in data["targets"]:
        if t["name"].lower() == name.lower():
            return t
    return None


def get_bioactivities_for_target(data: dict, target_chembl_id: str) -> list[dict]:
    """Get all bioactivity records for a given target."""
    return [b for b in data["bioactivities"] if b["target"] == target_chembl_id]


def get_bioactivities_for_compound(data: dict, compound_chembl_id: str) -> list[dict]:
    """Get all bioactivity records for a given compound."""
    return [b for b in data["bioactivities"] if b["compound"] == compound_chembl_id]
