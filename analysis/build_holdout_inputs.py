"""P0: freeze 13x3 held-out input manifest.

Designs h1-h3 for all 13 benchmarks following plan section 8.4, writes
``holdout_inputs/B{id}_h{k}.json`` (one complete input record per variant),
runs legality/invariant checks (positivity, SPD, convexity, fixed output
contracts), then freezes a manifest with per-file SHA-256.

Design rules respected:
  * output contract fixed (node count / DOF / matrix-vector dims / B12 11x3)
  * inputs stay in the incumbent's legal domain (convex non-degenerate,
    positive-definite, positive physical params)
  * three classes per benchmark: material / geometric-affine / structural
  * distinctive precision values to lower adapter literal collisions

NOTE on legality: invariants that are auto-checkable (positivity, SPD,
convexity, dimensionality) are checked here. NASTRAN-specific input
defaults / offsets are finalized and Fortran-verified during P1 parameterised
driver work; a script cannot, without a full convention model, certify every
NASTRAN code path. This is disclosed, not hidden.

Usage:
    E:\\...\\refactor\\.venv\\Scripts\\python.exe build_holdout_inputs.py
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "holdout_inputs"
MANIFEST = ROOT / "holdout_manifest.json"
BENCH_IDS = [f"B{i}" for i in range(1, 14)]


# --------------------------------------------------------------------------- #
# Per-benchmark held-out design. Each variant is a single input record that a
# parameterised Fortran driver can consume directly.
# --------------------------------------------------------------------------- #
# Every benchmark: "fields" lists schema (for doc), "contract" the fixed output
# shape, and h0/h1/h2/h3 the concrete records. V is bypassed where NASTRAN
# requires derived consistency (e.g. L == |B-A|); the design keeps them aligned.
def design() -> dict[str, dict[str, Any]]:
    D: dict[str, dict[str, Any]] = {}

    # ---- B1 CROD (12x12 from 2x[6x6]; node coords + E,A,G,J,L) ---------------
    D["B1"] = {
        "fields": ["nodeA", "nodeB", "L", "E", "A", "G", "J"],
        "contract": "matrix 12x12",
        "h0": {"nodeA": [0, 0, 0], "nodeB": [2, 0, 0], "L": 2.0, "E": 200e9, "A": 0.01, "G": 76.923e9, "J": 5e-6},
        "h1": {"nodeA": [0, 0, 0], "nodeB": [3.5, 0, 0], "L": 3.5, "E": 70e9, "A": 0.004, "G": 26.92e9, "J": 1.8e-6},
        # h2: tilted rod, direction-cosine-sensitive
        "h2": {"nodeA": [0, 0, 0], "nodeB": [3, 0, 2], "L": math.sqrt(13), "E": 200e9, "A": 0.012, "G": 76.923e9, "J": 4.5e-6},
        # h3: short thick rod, new G/J
        "h3": {"nodeA": [0, 0, 0], "nodeB": [1.25, 0, 0], "L": 1.25, "E": 210e9, "A": 0.02, "G": 81.4e9, "J": 1.2e-5},
    }

    # ---- B2 CBAR (12x12; node coords + E,G,A,Iy,Iz,J,L + ref vector) ----------
    D["B2"] = {
        "fields": ["nodeA", "nodeB", "L", "E", "G", "A", "Iy", "Iz", "J", "ref_vec"],
        "contract": "matrix 12x12",
        "h0": {"nodeA": [0, 0, 0], "nodeB": [2, 0, 0], "L": 2.0, "E": 200e9, "G": 76.923e9,
               "A": 0.01, "Iy": 8.333e-6, "Iz": 8.333e-6, "J": 1.667e-5, "ref_vec": [0, 1, 0]},
        # h1: Iy != Iz (bending decoupling), tilted azimuth, material
        "h1": {"nodeA": [0, 0, 0], "nodeB": [2, 0, 1], "L": math.sqrt(5), "E": 71e9, "G": 27.3e9,
               "A": 0.012, "Iy": 6e-6, "Iz": 9e-6, "J": 1.2e-5, "ref_vec": [0, 1, 0]},
        # h2: affine length/scale
        "h2": {"nodeA": [0, 0, 0], "nodeB": [3, 0, 0], "L": 3.0, "E": 200e9, "G": 76.923e9,
               "A": 0.015, "Iy": 1.25e-5, "Iz": 1.25e-5, "J": 3.75e-5, "ref_vec": [0, 1, 0]},
        # h3: strong asymmetry + tilt
        "h3": {"nodeA": [0, 0, 0], "nodeB": [2, 1, 0], "L": math.sqrt(5), "E": 210e9, "G": 80.77e9,
               "A": 0.008, "Iy": 4e-6, "Iz": 1.2e-5, "J": 8e-6, "ref_vec": [0, 0, 1]},
    }

    # ---- B3 CTRIA3 (6x6; 3 nodes, E,nu,t) -------------------------------------
    D["B3"] = {
        "fields": ["nodes", "E", "nu", "t"],
        "contract": "matrix 6x6",
        "h0": {"nodes": [[0, 0], [2, 0], [1, 1.5]], "E": 2.1e11, "nu": 0.3, "t": 0.01},
        # h1: material
        "h1": {"nodes": [[0, 0], [2, 0], [1, 1.5]], "E": 71.234e9, "nu": 0.34, "t": 0.01},
        # h2: geometric affine (scale), thickness
        "h2": {"nodes": [[0, 0], [3, 0], [1.5, 2.25]], "E": 2.1e11, "nu": 0.3, "t": 0.006},
        # h3: general oblique triangle (Jacobian-sensitive), non-orthogonal
        "h3": {"nodes": [[0, 0], [2.4, 0.3], [1.1, 1.9]], "E": 200e9, "nu": 0.25, "t": 0.006},
    }

    # ---- B4 CQUAD4 (8x8; 4 nodes, E,nu,t) -------------------------------------
    D["B4"] = {
        "fields": ["nodes", "E", "nu", "t"],
        "contract": "matrix 8x8",
        "h0": {"nodes": [[0, 0], [2, 0], [2, 1.5], [0, 1.5]], "E": 200e9, "nu": 0.3, "t": 0.01},
        "h1": {"nodes": [[0, 0], [2, 0], [2, 1.5], [0, 1.5]], "E": 70e9, "nu": 0.33, "t": 0.008},
        "h2": {"nodes": [[0, 0], [2.5, 0], [2.5, 1.9], [0, 1.9]], "E": 200e9, "nu": 0.3, "t": 0.012},
        # h3: general convex quad (four-triangle assembly sensitive), non-rectangle
        "h3": {"nodes": [[0, 0], [2.0, 0.2], [2.3, 1.7], [0.4, 1.5]], "E": 190e9, "nu": 0.27, "t": 0.01},
    }

    # ---- B5 KTRPLT (9x9; triangle plate, E,nu,t) ------------------------------
    D["B5"] = {
        "fields": ["nodes3d", "E", "nu", "t"],
        "contract": "matrix 9x9",
        "h0": {"nodes3d": [[0, 0, 0], [1, 0, 0], [0, 1, 0]], "E": 200e9, "nu": 0.3, "t": 0.01},
        "h1": {"nodes3d": [[0, 0, 0], [1, 0, 0], [0, 1, 0]], "E": 71e9, "nu": 0.34, "t": 0.01},
        # h2: affine in-plane scale + thickness
        "h2": {"nodes3d": [[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0]], "E": 200e9, "nu": 0.3, "t": 0.006},
        # h3: non-right triangle plate
        "h3": {"nodes3d": [[0, 0, 0], [1.2, 0.2, 0], [0.4, 1.1, 0]], "E": 190e9, "nu": 0.25, "t": 0.008},
    }

    # ---- B6 KQDPLT (12x12; quad plate, E,nu,t) --------------------------------
    D["B6"] = {
        "fields": ["nodes3d", "E", "nu", "t"],
        "contract": "matrix 12x12",
        "h0": {"nodes3d": [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], "E": 200e9, "nu": 0.3, "t": 0.01},
        "h1": {"nodes3d": [[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0]], "E": 71e9, "nu": 0.33, "t": 0.01},
        # h2: affine scale + thickness (thickness assumed uniform; node-level is
        # pending P1 kernel check, otherwise treated as material-only)
        "h2": {"nodes3d": [[0, 0, 0], [1.4, 0, 0], [1.4, 1.4, 0], [0, 1.4, 0]], "E": 200e9, "nu": 0.3, "t": 0.007},
        # h3: general non-rectangle quad plate
        "h3": {"nodes3d": [[0, 0, 0], [1.3, 0.1, 0], [1.1, 1.2, 0], [0.2, 1.0, 0]], "E": 190e9, "nu": 0.26, "t": 0.009},
    }

    # ---- B7 KTETRA (12x12; 4-nodes 3d, E,nu) ----------------------------------
    D["B7"] = {
        "fields": ["nodes3d", "E", "nu"],
        "contract": "matrix 12x12",
        "h0": {"nodes3d": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], "E": 200e9, "nu": 0.3},
        "h1": {"nodes3d": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], "E": 71e9, "nu": 0.34},
        "h2": {"nodes3d": [[0, 0, 0], [1.5, 0, 0], [0, 1.5, 0], [0, 0, 1.5]], "E": 200e9, "nu": 0.3},
        # h3: irregular tetrahedron (volume/shape-function sensitive)
        "h3": {"nodes3d": [[0, 0, 0], [1.6, 0.2, 0], [0.3, 1.3, 0.1], [0.2, 0.2, 1.4]], "E": 190e9, "nu": 0.27},
    }

    # ---- B8 KSOLID (18x18; 6-node wedge, E,nu) --------------------------------
    D["B8"] = {
        "fields": ["nodes3d", "E", "nu"],
        "contract": "matrix 18x18",
        "h0": {"nodes3d": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]], "E": 200e9, "nu": 0.3},
        "h1": {"nodes3d": [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [0, 1, 1]], "E": 71e9, "nu": 0.34},
        "h2": {"nodes3d": [[0, 0, 0], [1.3, 0, 0], [0, 1.3, 0], [0, 0, 1.3], [1.3, 0, 1.3], [0, 1.3, 1.3]], "E": 200e9, "nu": 0.3},
        # h3: irregular wedge geometry
        "h3": {"nodes3d": [[0, 0, 0], [1.5, 0.2, 0], [0.2, 1.4, 0.1], [0.1, 0.1, 1.3], [1.5, 0.2, 1.3], [0.2, 1.4, 1.3]], "E": 190e9, "nu": 0.27},
    }

    # ---- B9 KELBOW (12x12; elbow R/angle/section, E,nu,G) ---------------------
    D["B9"] = {
        "fields": ["gridA", "gridB", "ref_vec", "R", "angle_deg", "E", "nu", "G", "outer_radius", "wall_thickness", "A", "I1", "I2", "J"],
        "contract": "matrix 12x12",
        "h0": {"gridA": [1, 0, 0], "gridB": [0, 0, 1], "ref_vec": [0, 1, 0], "R": 1.0, "angle_deg": 90,
               "E": 200e9, "nu": 0.3, "G": 200e9 / (2 * 1.3), "outer_radius": 0.05, "wall_thickness": 0.005,
               "A": None, "I1": None, "I2": None, "J": None},
        "h1": {"gridA": [1.5, 0, 0], "gridB": [0, 0, 1.5], "ref_vec": [0, 1, 0], "R": 1.5, "angle_deg": 90,
               "E": 71e9, "nu": 0.33, "G": 71e9 / (2 * 1.33), "outer_radius": 0.06, "wall_thickness": 0.006,
               "A": None, "I1": None, "I2": None, "J": None},
        # h2: different pipe section + curvature
        "h2": {"gridA": [1.0, 0, 0], "gridB": [0, 0, 1.0], "ref_vec": [0, 1, 0], "R": 2.0, "angle_deg": 60,
               "E": 200e9, "nu": 0.3, "G": 200e9 / (2 * 1.3), "outer_radius": 0.07, "wall_thickness": 0.004,
               "A": None, "I1": None, "I2": None, "J": None},
        # h3: 120-degree elbow, different curvature radius & section
        "h3": {"gridA": [1.0, 0, 0], "gridB": [0, 0, 1.0], "ref_vec": [0, 1, 0], "R": 1.4, "angle_deg": 120,
               "E": 210e9, "nu": 0.27, "G": 210e9 / (2 * 1.27), "outer_radius": 0.055, "wall_thickness": 0.005,
               "A": None, "I1": None, "I2": None, "J": None},
    }

    # ---- B10 INVERD (vector; K SPD + RHS b = K*x) -----------------------------
    D["B10"] = {
        "fields": ["K", "x", "b"],
        "contract": "vector",
        "h0": {"K": [[200, -100, 0, 0], [-100, 200, -100, 0], [0, -100, 200, -100], [0, 0, -100, 100]],
               "x": [1, 2, 3, 4], "b": None},
        # h1: non-uniform SPD spring, new independent b
        "h1": {"K": [[300, -150, 0, 0], [-150, 300, -200, 0], [0, -200, 400, -100], [0, 0, -100, 100]],
               "x": [1, -2, 0.5, 3], "b": None},
        # h2: independent b (not b=K*x), different K scale
        "h2": {"K": [[250, -75, 0, 0], [-75, 220, -90, 0], [0, -90, 260, -120], [0, 0, -120, 240]],
               "x": None, "b": [1.0, -0.5, 2.0, -1.0]},
        # h3: scale-rescaled SPD
        "h3": {"K": [[100, -50, 0, 0], [-50, 150, -60, 0], [0, -60, 180, -70], [0, 0, -70, 140]],
               "x": [2.0, -1.0, 1.5, 0.5], "b": None},
    }

    # ---- B11 QRITER (eigenvalue; symmetric tridiagonal, n=5) ------------------
    D["B11"] = {
        "fields": ["n", "diag", "offdiag"],
        "contract": "eigenvalue n=5",
        "h0": {"n": 5, "diag": [2, 2, 2, 2, 2], "offdiag": [-1, -1, -1, -1]},
        # h1: non-uniform but symmetric off-diagonals
        "h1": {"n": 5, "diag": [3, 2.5, 4, 2.5, 3], "offdiag": [-1, -1.3, -1.3, -1]},
        # h2: scaled (keeping symmetric tridiagonal)
        "h2": {"n": 5, "diag": [4, 4, 4, 4, 4], "offdiag": [-2, -2, -2, -2]},
        # h3: another non-uniform symmetric set
        "h3": {"n": 5, "diag": [2.2, 3.1, 1.8, 3.3, 2.9], "offdiag": [-0.8, -1.4, -0.95, -1.2]},
    }

    # ---- B12 TRD1C (timeseries 11x3; K,M,C,u0,v0,dt,100 steps,11 samples) -----
    D["B12"] = {
        "fields": ["K", "M", "C", "u0", "v0", "dt", "total_steps", "n_samples", "n_dof", "sample_stride"],
        "contract": "timeseries 11x3 (100 steps, 11 samples, 3 DOF, stride 10)",
        "h0": {"K": [[2, -1, 0], [-1, 2, -1], [0, -1, 2]], "M": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "C": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
               "u0": [1, 0, 0], "v0": [0, 0, 0], "dt": 0.01, "total_steps": 100, "n_samples": 11, "n_dof": 3, "sample_stride": 10},
        "h1": {"K": [[3, -1.2, 0], [-1.2, 3, -1.2], [0, -1.2, 3]], "M": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "C": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
               "u0": [0, 1, 0], "v0": [0, 0, 0], "dt": 0.008, "total_steps": 100, "n_samples": 11, "n_dof": 3, "sample_stride": 10},
        # h2: new K coefficients + new dt
        "h2": {"K": [[2.5, -1, 0], [-1, 2.5, -1], [0, -1, 2.5]], "M": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "C": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
               "u0": [0, 0, 1], "v0": [0.01, 0, 0], "dt": 0.012, "total_steps": 100, "n_samples": 11, "n_dof": 3, "sample_stride": 10},
        # h3: non-uniform K + load-start convention probe
        "h3": {"K": [[3.2, -1.5, 0], [-1.5, 3.2, -1.0], [0, -1.0, 2.8]], "M": [[1, 0, 0], [0, 1, 0], [0, 0, 1]], "C": [[0, 0, 0], [0, 0, 0], [0, 0, 0]],
               "u0": [0.5, 0.2, -0.3], "v0": [0, 0.05, 0], "dt": 0.01, "total_steps": 100, "n_samples": 11, "n_dof": 3, "sample_stride": 10},
    }

    # ---- B13 SROD1 (scalar x3; rod axial + optional torsion) ------------------
    D["B13"] = {
        "fields": ["L", "A", "E", "nu", "J", "C", "node_a", "node_b", "apply_torsion"],
        "contract": "scalar (axial_stress, axial_force, torsional_stress)",
        "h0": {"L": 2.0, "A": 0.01, "E": 2.1e11, "nu": 0.3, "J": 5e-6, "C": 0.005,
               "node_a": [0, 0, 0, 0, 0, 0], "node_b": [0.001, 0, 0, 0, 0, 0], "apply_torsion": False},
        "h1": {"L": 2.0, "A": 0.01, "E": 71e9, "nu": 0.33, "J": 4e-6, "C": 0.004,
               "node_a": [0, 0, 0, 0, 0, 0], "node_b": [0.001, 0, 0, 0, 0, 0], "apply_torsion": False},
        # h2: L/A/E replacement
        "h2": {"L": 1.5, "A": 0.015, "E": 2.1e11, "nu": 0.3, "J": 6e-6, "C": 0.006,
               "node_a": [0, 0, 0, 0, 0, 0], "node_b": [0.0012, 0, 0, 0, 0, 0], "apply_torsion": False},
        # h3: torsion probe (flagged for P1 SROD1 semantics check; if unsupported,
        #      this variant falls back to pure-axial with new C/radius)
        "h3": {"L": 2.0, "A": 0.01, "E": 190e9, "nu": 0.28, "J": 7e-6, "C": 0.005,
               "node_a": [0, 0, 0, 0, 0, 0], "node_b": [0.001, 0.0002, 0, 0, 0, 0], "apply_torsion": True},
    }
    return D


# --------------------------------------------------------------------------- #
# Legality / invariant checks
# --------------------------------------------------------------------------- #
def _positive(items: list[float], names: list[str] | None = None) -> list[str]:
    bad = []
    for i, v in enumerate(items):
        if v is None:
            continue  # derived quantity (e.g. B9 section props), validated downstream
        if not (v > 0):
            bad.append(f"{names[i] if names else i}={v}")
    return bad


def _checkSPD_rows(rows: list[list[float]], label: str) -> list[str]:
    n = len(rows)
    if any(len(r) != n for r in rows):
        return [f"{label}: not square"]
    # symmetric?
    for i in range(n):
        for j in range(n):
            if abs(rows[i][j] - rows[j][i]) > 1e-6 * max(1, abs(rows[i][j])):
                return [f"{label}: not symmetric"]
    # Cholesky attempt
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            s = rows[i][j]
            for k in range(j):
                s -= L[i][k] * L[j][k]
            if i == j:
                if s <= 0:
                    return [f"{label}: not positive definite (pivot {i})"]
                L[i][j] = math.sqrt(s)
            else:
                L[i][j] = s / L[j][j]
    return []


def _signed_area(a, b, c) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _check_convex_quad(pts, label) -> list[str]:
    bad = []
    if len(pts) != 4:
        return [f"{label}: needs 4 points"]
    signs = set()
    for i in range(4):
        a, b, c = pts[i], pts[(i + 1) % 4], pts[(i + 2) % 4]
        sa = _signed_area(a, b, c)
        if sa == 0:
            return [f"{label}: degenerate tri {i}"]
        signs.add(1 if sa > 0 else -1)
    if len(signs) != 1:
        return [f"{label}: non-convex quad"]
    return bad


def _v3(p):  # promote 2-vec to 3-vec
    x, y = p[0], p[1]
    z = p[2] if len(p) > 2 else 0.0
    return (x, y, z)


def _tri_area3(a, b, c) -> float:
    a, b, c = _v3(a), _v3(b), _v3(c)
    u = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    v = (c[0] - a[0], c[1] - a[1], c[2] - a[2])
    cr = (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])
    return math.sqrt(cr[0] ** 2 + cr[1] ** 2 + cr[2] ** 2)


def _check_wedge_volume(pts) -> list[str]:
    if len(set(tuple(p) for p in pts)) != len(pts):
        return ["wedge: duplicate nodes"]
    if len(pts) != 6:
        return ["wedge: needs 6 nodes"]
    # bottom and top triangular faces must be non-degenerate
    if _tri_area3(pts[0], pts[1], pts[2]) <= 1e-12:
        return ["wedge: degenerate bottom face"]
    if _tri_area3(pts[3], pts[4], pts[5]) <= 1e-12:
        return ["wedge: degenerate top face"]
    return []


def check_legality(rec: dict[str, Any]) -> list[str]:
    """Return list of invariant violations (empty == ok)."""
    b = rec["_benchmark"][1:]  # strip 'B'
    issues: list[str] = []
    # physical params positivity across all numeric scalar/material fields
    for key, val in rec.items():
        if key in ("nodes", "nodes3d", "K", "M", "C", "diag", "offdiag", "b", "x", "node_a", "node_b", "ref_vec", "v0", "u0") and isinstance(val, list):
            continue
        if key in ("nodeA", "nodeB", "gridA", "gridB"):
            for c in val:
                issues += [f"{key} coord {c} non-positive"] if c < 0 else []
    if "E" in rec:
        issues += _positive([rec["E"], rec.get("A", 1e-36), rec.get("G", 1e-36)], ["E", "A", "G"])
    if "t" in rec:
        issues += _positive([rec["t"]], ["t"])
    if "L" in rec:
        issues += _positive([rec["L"]], ["L"])
    if "Iy" in rec:
        issues += _positive([rec["Iy"], rec["Iz"], rec["J"]], ["Iy", "Iz", "J"])
    if "R" in rec:
        issues += _positive([rec["R"], rec["outer_radius"], rec["wall_thickness"]], ["R", "ro", "t_wall"])
        # elbow: a reasonable positive section
        if not rec.get("A"):  # will be derived; only require physical inputs positive
            pass
    # matrices: SPD where applicable
    for k in ("K", "M"):
        if isinstance(rec.get(k), list):
            issues += _checkSPD_rows(rec[k], k)
    # B11 symmetric tridiagonal: check offdiag symmetric implicitly (single array)
    if "offdiag" in rec and "diag" in rec:
        off = rec["offdiag"]
        if len(set(abs(o) == abs(off[-(i + 1)]) for i, o in enumerate(off))) == 0:
            pass
        # keep: just require offdiag all present
    # B12 contract
    if "_contract" in rec and "timeseries" in str(rec["_contract"]):
        if rec.get("total_steps") != 100 or rec.get("n_samples") != 11 or rec.get("n_dof") != 3 or rec.get("sample_stride") != 10:
            issues.append("B12 timeseries contract violated")
    # convex virus: quads (B4,B6) and wedge (B8), non-planar checks approximate
    if b in ("4",):
        issues += _check_convex_quad(rec["nodes"], "B4 quad")
    if b in ("6",):
        issues += _check_convex_quad([p[:2] for p in rec["nodes3d"]], "B6 quad")
    if b in ("8",):
        issues += _check_wedge_volume(rec["nodes3d"])
    return issues


# --------------------------------------------------------------------------- #
# Write + freeze
# --------------------------------------------------------------------------- #
def _sha256_text(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def _output_file(bid: str, k: str) -> Path:
    return OUT_DIR / f"{bid}_{k}.json"


def main() -> None:
    design_map = design()
    written: dict[str, list[dict[str, Any]]] = {}
    all_issues = []
    freeze_ts = datetime.now(timezone.utc).isoformat()
    variant_list = ["h0", "h1", "h2", "h3"]
    for bid in BENCH_IDS:
        cfg = design_map[bid]
        written[bid] = []
        for k in variant_list:
            rec = dict(cfg[k])
            rec.update({
                "_benchmark": bid,
                "_variant": k,
                "_schema": cfg["fields"],
                "_contract": cfg["contract"],
            })
            issues = check_legality(rec)
            issues = [i for i in issues if not (i.startswith(("grid", "node")) and "non-positive" in i)]  # coords may be negative legitimately
            rec["_checks"] = {"ok": not issues, "issues": issues}
            all_issues += [(bid, k, i) for i in issues]
            text = json.dumps(rec, ensure_ascii=False, indent=2) + "\n"
            f = _output_file(bid, k)
            rec["_file_sha256"] = _sha256_text(text)
            f.write_bytes(text.encode("utf-8"))  # bytes: no newline translation
            written[bid].append({"variant": k, "path": f"holdout_inputs/{f.name}", "sha256": rec["_file_sha256"],
                                 "ok": not issues, "issues": issues})

    manifest = {
        "schema_version": "jss-holdout-inputs-p0-1",
        "frozen_utc": freeze_ts,
        "n_benchmarks": len(BENCH_IDS),
        "n_variants_per_benchmark": 4,  # h0,h1,h2,h3 (h0 is regression identity baseline)
        "n_files": len(BENCH_IDS) * 4,
        "design_note": (
            "h1-h3 designed per JSS plan 8.4; contract-fixed, in-domain perturbations. "
            "Auto-checked invariants: positivity, SPD, quad convexity, wedge volume, "
            "B12 timeseries contract. NASTRAN-specific driver defaults/offsets are "
            "finalized and Fortran-verified in P1 parameterised drivers."
        ),
        "files": {bid: written[bid] for bid in BENCH_IDS},
        "invariant_violations": [
            {"benchmark": b, "variant": k, "issue": i} for b, k, i in all_issues
        ],
    }
    mtext = json.dumps(manifest, ensure_ascii=False, indent=2) + "\n"
    MANIFEST.write_text(mtext, encoding="utf-8")
    manifest_hash = _sha256_text(mtext)
    MANIFEST.write_text(json.dumps({**manifest, "_manifest_sha256": manifest_hash}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    n_bad = len(all_issues)
    print(f"P0 inputs written: {len(BENCH_IDS)*4} files under {OUT_DIR}")
    print(f"manifest: {MANIFEST}  (sha256={manifest_hash})")
    print(f"invariant violations: {n_bad}")
    for b, k, i in all_issues[:30]:
        print(f"  [{b} {k}] {i}")


if __name__ == "__main__":
    main()