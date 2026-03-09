"""Backward compatibility — use curated_data instead."""
from chempath.data.curated_data import CURATED_DATA as MOCK_DATA, load_curated_data as load_mock_data

__all__ = ["MOCK_DATA", "load_mock_data"]
