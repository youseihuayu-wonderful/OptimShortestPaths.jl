"""
Streamlit Cloud entry point.

Adds ChemPath to sys.path and runs the dashboard.
This file exists solely for Streamlit Cloud deployment —
the actual app code lives in ChemPath/chempath/ui/dashboard.py.
"""
import sys
from pathlib import Path

# Add ChemPath to path so 'from chempath...' imports work
sys.path.insert(0, str(Path(__file__).parent / "ChemPath"))

# Streamlit Cloud looks for this file and runs it directly.
# We exec the real dashboard so all st.* calls happen in this process.
exec(open(Path(__file__).parent / "ChemPath" / "chempath" / "ui" / "dashboard.py").read())
