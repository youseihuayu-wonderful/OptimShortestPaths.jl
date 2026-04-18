"""
Streamlit Cloud entry point for ChemPath Drug Discovery Dashboard.

This thin wrapper sets up paths and patches DATA_PATH so the dashboard
works when deployed from the repo root on Streamlit Cloud.
"""
import sys
from pathlib import Path

# Repo root
ROOT = Path(__file__).resolve().parent

# Add ChemPath to sys.path
sys.path.insert(0, str(ROOT / "ChemPath"))

# Patch DATA_PATH before dashboard.py tries to use it:
# dashboard.py uses Path(__file__).parent.parent.parent / "data" / "chembl_real.json"
# but __file__ inside exec would point to this file, not dashboard.py.
# So we patch the module-level constant by setting it in the exec namespace.

import chempath.data.chembl_client  # noqa: E402
import chempath.ui  # noqa: E402

# Read and exec the dashboard code with a patched __file__
dashboard_path = ROOT / "ChemPath" / "chempath" / "ui" / "dashboard.py"
dashboard_code = dashboard_path.read_text()

# Replace the DATA_PATH line with the correct absolute path
data_path = ROOT / "ChemPath" / "data" / "chembl_real.json"
dashboard_code = dashboard_code.replace(
    'DATA_PATH = Path(__file__).parent.parent.parent / "data" / "chembl_real.json"',
    f'DATA_PATH = Path("{data_path}")',
)

exec(compile(dashboard_code, str(dashboard_path), "exec"))
