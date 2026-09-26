"""
elmer_env.py - Cross-platform environment and path resolver for Elmer FEM tools.
Supports Linux, macOS, and Windows. Respects ELMER_SOLVER and ELMERGRID_BIN env vars,
with PATH lookups via shutil.which and OS-specific fallbacks.
"""

import os
import shutil
from pathlib import Path


def get_elmer_solver() -> str:
    """Resolve the path to the ElmerSolver executable."""
    env_bin = os.environ.get("ELMER_SOLVER")
    if env_bin and (shutil.which(env_bin) or os.path.isfile(env_bin)):
        return env_bin

    which_bin = shutil.which("ElmerSolver")
    if which_bin:
        return which_bin

    candidates = [
        r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe",
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerSolver.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerSolver.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerSolver.exe",
        "/usr/local/bin/ElmerSolver",
        "/usr/bin/ElmerSolver",
        "/opt/elmerfem/bin/ElmerSolver",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c

    return "ElmerSolver"


def get_elmergrid_bin() -> str:
    """Resolve the path to the ElmerGrid executable."""
    env_bin = os.environ.get("ELMERGRID_BIN")
    if env_bin and (shutil.which(env_bin) or os.path.isfile(env_bin)):
        return env_bin

    which_bin = shutil.which("ElmerGrid")
    if which_bin:
        return which_bin

    candidates = [
        r"C:\Users\bresc\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe",
        os.path.expanduser(r"~\ElmerFEM\ElmerFEM-nogui-nompi-Windows-AMD64\bin\ElmerGrid.exe"),
        r"C:\Program Files\Elmer 9.0-Release\bin\ElmerGrid.exe",
        r"C:\Program Files (x86)\Elmer\bin\ElmerGrid.exe",
        "/usr/local/bin/ElmerGrid",
        "/usr/bin/ElmerGrid",
        "/opt/elmerfem/bin/ElmerGrid",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c

    return "ElmerGrid"
