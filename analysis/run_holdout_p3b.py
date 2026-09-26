"""P3b: offline held-out revalidation of retention-tested candidates (h1-h3).

Zero API cost - local MSVC compile + output comparison only.

For each of the 63 ``retention_tested`` candidates frozen in
``jss_holdout_adapter_scan_v2.json`` we adapt the h0 input to each held-out
variant, recompile with the frozen MSVC+Eigen toolchain, and compare output to
the frozen oracle ``_result`` with canonical-scoring tolerance semantics.

Adaptation tiering (honest boundary of *value* adaptation):
  clean    - only value-unique scalar params change (E,A,G,J,Iy,Iz,t,nu,C).
             Pure token-value substitution.
  geometry - node/geometry coordinates also change (B3/B4/B7 h2-h3, B1, B2).
             Param substitution + best-effort Vector3d node-coordinate rebind.
  limited  - structural/array/split-section input (B2 Iy/Iz split, B11
             diag/offdiag, B13 h2/h3 length+torsion+2D displacement). Not
             representable by literal substitution; recorded without compile.

Outcome per (candidate, variant):
  pass / fail / compile_error / run_error / parse_error / no_change /
  no_target_literal / geometry_unlocatable / adapter_limited
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util as _ilu
import json
import os
import re
import shutil
import sys
import tempfile
import types
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "analysis"
HOLDIN = ROOT / "holdout_inputs"
HOLDOUT = ROOT / "holdout_oracles"
BASELINE = PAPER / "jss_v4_holdout_baseline.json"
SCAN = PAPER / "jss_holdout_adapter_scan_v2.json"
OUT_JSON = PAPER / "holdout_results.json"

sys.path.insert(0, str(ROOT))

# Offline safety (same as the canonical rescore): never construct an API client.
if "openai" not in sys.modules:
    offline_openai = types.ModuleType("openai")

    class _OfflineOpenAI:
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            raise RuntimeError("LLM/API disabled in offline P3b revalidation")

    offline_openai.OpenAI = _OfflineOpenAI
    sys.modules["openai"] = offline_openai

from run_multilang_v2 import compile_and_run_cpp  # noqa: E402
from run_experiments_v3 import (  # noqa: E402
    load_benchmark_config, parse_output, verify_output,
)

_spec = _ilu.spec_from_file_location("holdout_adapter_v2", PAPER / "holdout_adapter_v2.py")
_adap = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_adap)
signed_spans = _adap.signed_spans
vnorm = _adap._vnorm

VARIANTS = ["h1", "h2", "h3"]

# Value-unique scalar roles we substitute by raw token value.
SAFE_FIELDS = {
    "B1":  ["E", "A", "G", "J"],
    "B2":  ["E", "G", "A", "J"],       # Iy/Iz split excluded -> limited when they differ
    "B3":  ["E", "nu", "t"],
    "B4":  ["E", "nu", "t"],
    "B7":  ["E", "nu"],
    "B13": ["E", "A", "J", "C"],
    "B11": [],
}
GEOMETRY_CAPABLE = {"B1", "B2", "B3", "B4", "B7"}

# Per-benchmark geometry: which input field(s) hold node coordinates.
_GEOM_SCALAR_KEY = {"B13": "L"}
_NODE_LIST = {
    "B3": {"key": "nodes", "dim": 2},
    "B4": {"key": "nodes", "dim": 2},
    "B7": {"key": "nodes3d", "dim": 3},
}
_NODE_PAIR = {"B1", "B2"}


def _geometry_changed(benchmark: str, h0: dict, hk: dict) -> bool:
    if benchmark in _NODE_LIST:
        k = _NODE_LIST[benchmark]["key"]
        return h0.get(k) != hk.get(k)
    if benchmark in _NODE_PAIR:
        return (h0.get("nodeA") != hk.get("nodeA")) or (h0.get("nodeB") != hk.get("nodeB"))
    if benchmark == "B13":
        return (h0.get("node_b") != hk.get("node_b")) or (h0.get("L") != hk.get("L"))
    return False


def classify(benchmark: str, h0: dict, hk: dict) -> tuple[str, dict[str, tuple]]:
    """Return (tier, safe_changes)."""
    safe_ch = compute_changes(h0, hk, SAFE_FIELDS[benchmark])
    geo = _geometry_changed(benchmark, h0, hk)
    split = _split_section_changed(benchmark, h0, hk)
    # Any non-safe, non-flagged field that differs means value-only substitution
    # cannot fully represent the perturbation (arrays, diag/offdiag, booleans, ..).
    non_safe_changed = _non_safe_changed(benchmark, h0, hk)
    if not geo and not split and not non_safe_changed:
        return ("clean", safe_ch)
    if benchmark == "B13" or benchmark == "B2" or benchmark == "B11":
        # B13: length/torsion/displacement geometry; B2: Iy/Iz split section;
        # B11: array-built matrix - none are literal-substitutable faithfully.
        return ("limited", {})
    # B1, B3, B4, B7: node coordinates are literal Vector3d -> best effort.
    return ("geometry", safe_ch)


def _non_safe_changed(benchmark: str, h0: dict, hk: dict) -> bool:
    safe = set(SAFE_FIELDS[benchmark])
    for k in ("nodeA", "nodeB", "node_a", "node_b", "nodes", "nodes3d",
              "ref_vec", "n", "apply_torsion", "diag", "offdiag"):
        if k in h0 or k in hk:
            if h0.get(k) != hk.get(k):
                return True
    for k in ("L", "Iy", "Iz", "nu"):
        if k in SAFE_FIELDS[benchmark]:
            continue
        if h0.get(k) != hk.get(k):
            return True
    return False


def _split_section_changed(benchmark: str, h0: dict, hk: dict) -> bool:
    if benchmark != "B2":
        return False
    return (h0.get("Iy") != hk.get("Iy")) or (h0.get("Iz") != hk.get("Iz"))


def literal(value: float) -> str:
    s = repr(float(value))
    if "e" not in s and "E" not in s and "." not in s:
        s += ".0"
    return s


def compute_changes(h0: dict, hk: dict, fields: list[str]) -> dict[str, tuple]:
    ch = {}
    for f in fields:
        a, b = h0.get(f), hk.get(f)
        if isinstance(a, (int, float)) and isinstance(b, (int, float)):
            if not np.isclose(float(a), float(b), rtol=0, atol=0):
                ch[f] = (a, b)
    return ch


def rebind_values(source: str, changes: dict[str, tuple]) -> tuple[str, int]:
    spans = signed_spans(source)
    repl = []
    for val_old, val_new in changes.values():
        old_norm = vnorm(f"{val_old:g}")
        new_lit = literal(val_new)
        for v, s, e in spans:
            if vnorm(v) == old_norm:
                repl.append((s, e, new_lit))
    if not repl:
        return source, 0
    out = list(source)
    for s, e, new_lit in sorted(repl, key=lambda x: -x[0]):
        out[s:e] = list(new_lit)
    return "".join(out), len(repl)


# ------------------------------ geometry rebinding -------------------------- #
_VEC3D = re.compile(
    r"(Vector3d|Vector3)\s+([A-Za-z_]\w*)\s*\(\s*"
    r"([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)"
)
_VAR_INDEX = {f"node{i}": i - 1 for i in range(1, 5)} | {f"n{i}": i - 1 for i in range(1, 5)}
_VAR_PAIR = {
    "nodea": "nodeA", "nodeb": "nodeB", "node1": "nodeA", "node2": "nodeB",
    "n1": "nodeA", "n2": "nodeB", "node_a": "node_a", "node_b": "node_b",
}
# Scalar node-coordinate form: "double x2 = 2.0" (coord letter + node index).
_SCALAR_NODE = re.compile(r"\b([xyz])([1-4])\s*=\s*([-+0-9.eE]+)")


def rebind_nodes(source: str, benchmark: str, inp_hk: dict) -> tuple[str, int, list[str]]:
    note: list[str] = []
    repl = []
    if benchmark in _NODE_LIST:
        cfg = _NODE_LIST[benchmark]
        nodes = inp_hk.get(cfg["key"]) or []
        for m in _VEC3D.finditer(source):
            idx = _VAR_INDEX.get(m.group(2).lower())
            if idx is None or idx >= len(nodes):
                continue
            c = nodes[idx]
            x, y = float(c[0]), float(c[1])
            z = float(c[2]) if cfg["dim"] == 3 else 0.0
            new = f"{literal(x)}, {literal(y)}, {literal(z)}"
            repl.append((m.start(3), m.end(5), new))
        note = [f"node-list rebound ({cfg['key']})"]
        # Scalar coordinate form: x1,x2,.. / y1,y2,.. / z1,z2,..
        col = {"x": 0, "y": 1, "z": 2}
        for m in _SCALAR_NODE.finditer(source):
            coord, idx = m.group(1), int(m.group(2)) - 1
            if idx >= len(nodes):
                continue
            c = nodes[idx]
            val = float(c[col[coord]]) if col[coord] < len(c) else 0.0
            repl.append((m.start(3), m.end(3), literal(val)))
            if "scalar-coord" not in " ".join(note):
                note.append("scalar-coord rebound")
    elif benchmark in _NODE_PAIR:
        for m in _VEC3D.finditer(source):
            role = _VAR_PAIR.get(m.group(2).lower())
            if role is None:
                continue
            v = inp_hk.get(role)
            if not isinstance(v, list) or len(v) < 3:
                continue
            new = f"{literal(v[0])}, {literal(v[1])}, {literal(v[2])}"
            repl.append((m.start(3), m.end(5), new))
        note = ["node-pair rebound (nodeA/nodeB)"]
    if not repl:
        return source, 0, ["no Vector3d node initialiser matched"]
    out = list(source)
    for s, e, new in sorted(repl, key=lambda x: -x[0]):
        out[s:e] = list(new)
    return "".join(out), len(repl), note


# ------------------------------ oracles ------------------------------------ #
def oracle_reference(config: dict, result: dict) -> np.ndarray:
    vt = config["verification_type"]
    if vt == "matrix_symmetry":
        return np.array(result["matrix"], dtype=float)
    if vt == "scalar_compare":
        return np.array([
            float(result.get("axial_stress", np.nan)),
            float(result.get("axial_force", np.nan)),
            float(result.get("torsional_stress", np.nan)),
        ], dtype=float)
    if vt == "eigenvalue_compare":
        return np.sort(np.array(result["eigenvalues"], dtype=float))
    if vt == "vector_compare":
        return np.array(result["vector"] if "vector" in result else result.get("solution"), dtype=float)
    if vt == "timeseries_compare":
        return np.array(result["timeseries"], dtype=float)
    return np.array(result.get("matrix"), dtype=float)


# ------------------------------ run one ------------------------------------ #
def run_one(benchmark, candidate_src_path, inp_h0, inp_hk, variant, oracle_result,
            config, work_root, index):
    tier, safe_ch = classify(benchmark, inp_h0, inp_hk)
    rec = {
        "index": index, "benchmark": benchmark, "variant": variant,
        "candidate_path": str(candidate_src_path.relative_to(ROOT).as_posix()),
        "tier": tier, "outcome": None, "rebind_tokens": 0, "rebind_notes": None,
        "compile": False, "run": False, "equivalence": False,
        "stdout_sha256": None, "error_stage": None, "error": None,
        "max_rel_diff": None, "max_abs_diff": None, "verification": None,
    }
    try:
        src = candidate_src_path.read_text(encoding="utf-8", errors="replace")
        if tier == "clean":
            if not safe_ch:
                rec["outcome"] = "no_change"
                return rec
            rebound, n = rebind_values(src, safe_ch)
            rec["rebind_tokens"] = n
            if n == 0:
                rec["outcome"] = "no_target_literal"
                rec["error"] = f"no token matched h0 params {list(safe_ch)}"
                return rec
        elif tier == "geometry":
            total = 0
            rebound = src
            if safe_ch:
                rebound, n = rebind_values(src, safe_ch)
                total += n
            rebound, g, notes = rebind_nodes(rebound, benchmark, inp_hk)
            total += g
            rec["rebind_tokens"] = total
            rec["rebind_notes"] = notes
            if g == 0:
                rec["outcome"] = "geometry_unlocatable"
                rec["error"] = "; ".join(notes)
                return rec
        else:  # limited
            rec["outcome"] = "adapter_limited"
            rec["error"] = "structural/array/split-section input not value-adaptable"
            return rec

        work = work_root / f"h{index:04d}_{variant}"
        work.mkdir(parents=True, exist_ok=True)
        run = compile_and_run_cpp(rebound, work)
        if not run.get("success"):
            stage = run.get("stage") or "compile"
            rec["compile"] = stage != "compile"
            rec["error_stage"] = stage
            rec["error"] = str(run.get("error", ""))[:2000]
            rec["outcome"] = "compile_error" if stage == "compile" else "run_error"
            return rec
        rec["compile"] = True
        rec["run"] = True
        stdout = run.get("stdout", "")
        rec["stdout_sha256"] = hashlib.sha256(stdout.encode("utf-8", errors="replace")).hexdigest()
        try:
            computed = parse_output(stdout, config)
        except Exception as exc:
            rec["error_stage"] = "parse"
            rec["error"] = f"{type(exc).__name__}: {exc}"[:2000]
            rec["outcome"] = "parse_error"
            return rec
        if computed is None:
            rec["error_stage"] = "parse"
            rec["error"] = "Could not parse output"
            rec["outcome"] = "parse_error"
            return rec
        try:
            ref = oracle_reference(config, oracle_result)
            verif = verify_output(computed, ref, config)
            rec["verification"] = json.loads(json.dumps(
                verif, default=lambda o: o.item() if hasattr(o, "item") else repr(o)))
            rec["equivalence"] = bool(verif.get("passed"))
            rec["max_rel_diff"] = verif.get("max_rel_diff")
            rec["max_abs_diff"] = verif.get("max_abs_diff")
            rec["outcome"] = "pass" if rec["equivalence"] else "fail"
            if not rec["equivalence"]:
                rec["error_stage"] = "equivalence"
                rec["error"] = verif.get("reason") or verif.get("error") or "numerical mismatch"
        except Exception as exc:
            rec["error_stage"] = "verify"
            rec["error"] = f"{type(exc).__name__}: {exc}"[:2000]
            rec["outcome"] = "verify_error"
        finally:
            shutil.rmtree(work, ignore_errors=True)
    except Exception as exc:
        rec["error_stage"] = "orchestration"
        rec["error"] = f"{type(exc).__name__}: {exc}"[:2000]
        rec["outcome"] = "run_error"
    return rec


def main(only: str | None = None, verbose: bool = False):
    scan = json.load(open(SCAN, encoding="utf-8"))
    ret = []
    for grp in scan["records"]:
        for c in grp.get("details", []):
            if c.get("final_status") == "retention_tested":
                ret.append(c)
    if only:
        ret = [c for c in ret if c["benchmark"] == only]

    results = []
    work_root = Path(tempfile.mkdtemp(prefix="p3b_"))
    idx = 0
    for c in ret:
        bench = c["benchmark"]
        config = load_benchmark_config(bench)
        cand_abs = ROOT / Path(c["candidate_path"].replace("/", os.sep))
        inp = {}
        oracle = {}
        inp["h0"] = json.load(open(HOLDIN / f"{bench}_h0.json", encoding="utf-8"))
        for v in VARIANTS:
            inp[v] = json.load(open(HOLDIN / f"{bench}_{v}.json", encoding="utf-8"))
            oracle[v] = json.load(open(HOLDOUT / f"{bench}_{v}.json", encoding="utf-8"))
        for v in VARIANTS:
            idx += 1
            rec = run_one(bench, cand_abs, inp["h0"], inp[v], v, oracle[v]["_result"],
                          config, work_root, idx)
            rec["canonical_key"] = c["canonical_key"]
            results.append(rec)
            if verbose:
                print(f"[{bench} {v}] {c['candidate_path']} -> {rec['outcome']}")
    shutil.rmtree(work_root, ignore_errors=True)

    merged = {}
    if OUT_JSON.exists():
        try:
            prev = json.load(open(OUT_JSON, encoding="utf-8"))
            for r in prev.get("results", []):
                merged[(r["benchmark"], r["candidate_path"], r["variant"])] = r
        except Exception:
            merged = {}
    for r in results:
        merged[(r["benchmark"], r["candidate_path"], r["variant"])] = r
    merged_list = list(merged.values())

    doc = {
        "schema_version": "jss-holdout-p3b-1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_recorded": len(merged_list),
        "note": (
            "P3b offline revalidation of 63 retention_tested candidates on h1-h3. "
            "clean = value-unique param substitution; geometry = param + best-effort "
            "Vector3d node rebinding; limited = structural/array/split-section input "
            "recorded without compile. Zero API calls."
        ),
        "results": merged_list,
    }
    OUT_JSON.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")
    print("saved", OUT_JSON)
    print("outcomes:", Counter(r["outcome"] for r in merged_list))
    per_b = defaultdict(Counter)
    for r in merged_list:
        per_b[r["benchmark"]][r["outcome"]] += 1
    for b in sorted(per_b):
        print(b, dict(per_b[b]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    ap.add_argument("--verbose", action="store_true")
    a = ap.parse_args()
    main(a.only, a.verbose)
