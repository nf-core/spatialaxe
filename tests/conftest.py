"""Pytest configuration for nf-xenium-processing tests."""

from __future__ import annotations

import sys
from pathlib import Path

# Add bin/ to path so skill_*.py modules can be imported
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bin"))
