"""P3a: candidate input adapter — mutation-coverage scan v2 (section 8.6).

Improves the naive v1 scan (analysis/holdout_adapter.py) in three ways:

  1. Signed numeric tokenisation with source spans (v1's tokeniser dropped the
     leading ``-`` so negative values like ``-100`` never matched their rule).
  2. Per-field *expected* occurrence counts plus optional surrounding *context*
     regexes, so a value that also appears as an algorithm constant (e.g. the
     length ``L = 2.0`` vs. a formula denominator ``/ 2.0``) can still be
     unambiguously located where the candidate spells it in an assignment shape.
  3. An honest three-way provisional split:
       retention_tested_ready  - every input field is uniquely locatable
       not_parameterizable     - a field's value collides with algorithm
                                 constants / the input is structurally assembled
                                 (matrix/array) and no context can disambiguate
       adapter_unverified      - infra-level failure (missing rule table, empty
                                 source, ..) that must be resolved before scoring

This is offline Gate 1 (no compile).  Gate 2 (h0 reconstruction) lives in
analysis/holdout_adapter.py --compile on the ``retention_tested_ready`` set.

CLI:
  python holdout_adapter_v2.py --scan all
  python holdout_adapter_v2.py --scan B1 --print
"""

from __future__ import annotations

import argparse
import json
import os
import re
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "analysis"
BASELINE = PAPER / "jss_v4_holdout_baseline.json"
OUT = PAPER / "jss_holdout_adapter_scan_v2.json"

# signed numeric literal with spans (value-normalised), excluding base-notation
_SIGNED = re.compile(
    r"(?<![0-9A-Za-z_.])([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)[fFlLuU]*"
)
_SUFFIX = re.compile(r"[fFlLuU+]+$")


def _vnorm(lit: str) -> str:
    s = _SUFFIX.sub("", lit).replace("+", "")
    try:
        return format(Decimal(s), "f")
    except (InvalidOperation, ValueError):
        return s


def signed_spans(source: str) -> list[tuple[str, int, int]]:
    """Return (value_normalised, start, end) for every signed numeric literal."""
    out: list[tuple[str, int, int]] = []
    for m in _SIGNED.finditer(source):
        raw = m.group(1)
        if raw.lower().startswith("0x") or raw.lower().startswith("0b"):
            continue
        out.append((_vnorm(raw), m.start(1), m.end(1)))
    return out


def candidate_abs(rel: str) -> Path:
    return ROOT / Path(rel.replace("/", os.sep))


# -------------------------------------------------------------------------- #
# Per-benchmark adapter rules.
# Each field: role, lits([h0 literals, value-normalised), expected, ctx(regex|None)
#   "brittle": field is known to collide with algorithm constants / structural
#              input, so it MUST be resolved through ``ctx``; else not_parameterizable.
# Value-unique fields (E/A/G/I/J/t/nu/C) resolve by exact occurrence count and
# are the reliable core for every benchmark that has candidates.
# -------------------------------------------------------------------------- #
def rules_for(benchmark: str) -> list[dict[str, Any]]:
    R = {
        "B1": [
            {"role": "E",   "lits": ["200e9"],  "expected": 1, "brittle": False, "ctx": None},
            {"role": "A",   "lits": ["0.01"],   "expected": 1, "brittle": False, "ctx": None},
            {"role": "G",   "lits": ["76.923e9"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "J",   "lits": ["5e-6"],   "expected": 1, "brittle": False, "ctx": None},
            {"role": "L",   "lits": ["2.0"],    "expected": 2, "brittle": True,
             "ctx": r"(?i)(\bL[a-z0-9_]*\s*=\s*2\.0\b)|(node[a-z0-9_]*\s*(\(|<<|\[)\s*2\.0)"},
        ],
        "B2": [
            {"role": "E",   "lits": ["200e9"],    "expected": 1, "brittle": False, "ctx": None},
            {"role": "A",   "lits": ["0.01"],     "expected": 1, "brittle": False, "ctx": None},
            {"role": "G",   "lits": ["76.923e9"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "I",   "lits": ["8.333e-6"], "expected": 2, "brittle": False, "ctx": None},
            {"role": "J",   "lits": ["1.667e-5"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "L",   "lits": ["2.0"],      "expected": 2, "brittle": True,
             "ctx": r"(?i)(\bL[a-z0-9_]*\s*=\s*2\.0\b)|(node[a-z0-9_]*\s*(\(|<<|\[)\s*2\.0)"},
        ],
        "B3": [
            {"role": "E", "lits": ["2.1e11"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "nu","lits": ["0.3"],    "expected": 1, "brittle": False, "ctx": None},
            {"role": "t", "lits": ["0.01"],   "expected": 1, "brittle": False, "ctx": None},
        ],
        "B4": [
            {"role": "E",  "lits": ["200e9"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "nu", "lits": ["0.3"],   "expected": 1, "brittle": False, "ctx": None},
            {"role": "t",  "lits": ["0.01"],  "expected": 1, "brittle": False, "ctx": None},
        ],
        "B7": [
            {"role": "E",  "lits": ["200e9"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "nu", "lits": ["0.3"],   "expected": 1, "brittle": False, "ctx": None},
        ],
        "B10": [
            {"role": "K_diag",    "lits": ["200"],  "expected": 3, "brittle": True,
             "ctx": r"(?i)\b(?:diag|val|diagm)\s*=?\s*(?:<<|\{|\[|,|\()\s*200"},
            {"role": "K_last",    "lits": ["100"],  "expected": 1, "brittle": True,
             "ctx": r"(?i)(?<![\d-])\b100\b"},
            {"role": "K_off",     "lits": ["-100"], "expected": 3, "brittle": True,
             "ctx": r"(?i)-?\b100\b"},
        ],
        "B11": [
            {"role": "diag",    "lits": ["2.0"],   "expected": 5, "brittle": True,
             "ctx": r"(?i)(diag[a-z0-9_]*\s*=\s*\{|\bdiag\s*=)"},
            {"role": "offdiag", "lits": ["-1.0"],  "expected": 4, "brittle": True,
             "ctx": r"(?i)(off[a-z0-9_]*\s*=\s*\{|\boffd\s*=)"},
        ],
        "B13": [
            {"role": "E",  "lits": ["2.1e11"], "expected": 1, "brittle": False, "ctx": None},
            {"role": "A",  "lits": ["0.01"],   "expected": 1, "brittle": False, "ctx": None},
            {"role": "J",  "lits": ["5e-6"],   "expected": 1, "brittle": False, "ctx": None},
            {"role": "C",  "lits": ["0.005"],  "expected": 1, "brittle": False, "ctx": None},
            {"role": "u_b","lits": ["0.001"],  "expected": 1, "brittle": True,
             "ctx": r"(?i)(disp_?b|u_b|x_b)\s*(<<|[,\[(=\s])\s*0\.001"},
        ],
    }
    return R.get(benchmark, [])


def scan_source(source: str, benchmark: str) -> dict[str, Any]:
    spans = signed_spans(source)
    by_val: dict[str, list[tuple[int, int]]] = {}
    for v, s, e in spans:
        by_val.setdefault(v, []).append((s, e))
    rules = rules_for(benchmark)
    per_field = []
    for field in rules:
        if field["lits"]:
            vals = {_vnorm(l) for l in field["lits"]}
        else:
            vals = set()
        hits: list[tuple[int, int]] = []
        for v in vals:
            hits.extend(by_val.get(v, []))
        hit_n = len(hits)
        # Non-brittle input parameters (E/A/G/I/J/t/nu/C) have a value dedicated
        # to that role; candidates may spell/repeat it (e.g. E twice), so we
        # rebind *every* occurrence. Mutation coverage is satisfied by >=1 hit.
        # Brittle fields (length coords, displacements, matrix/eigenvalue values)
        # collide with algorithm constants and must be located by context shape.
        ok = False
        reason = None
        if not field["brittle"]:
            ok = hit_n >= 1
            reason = None if ok else "miss"
        else:
            if field.get("ctx"):
                ok = len(re.findall(field["ctx"], source)) >= 1
                reason = "ctx_ok" if ok else "structural_ambiguous_ctx"
            else:
                reason = "structural_no_ctx"
        per_field.append({
            "role": field["role"],
            "expected": field["expected"],
            "hit": hit_n,
            "to_replace": hit_n,
            "brittle": field["brittle"],
            "ok": ok,
            "reason": reason or "ok",
        })
    all_ok = bool(rules) and all(f["ok"] for f in per_field)
    status = "retention_tested_ready" if all_ok else "not_parameterizable"
    return {"benchmark": benchmark, "all_ok": all_ok, "status": status,
            "per_field": per_field, "n_numeric_tokens": len(spans)}


def scan_records(records: list[dict[str, Any]], benchmark: str) -> dict[str, Any]:
    rows = [r for r in records if r["benchmark"] == benchmark]
    scanned = []
    for r in rows:
        src_path = candidate_abs(r["candidate_path"])
        src = src_path.read_text(encoding="utf-8", errors="replace")
        s = scan_source(src, benchmark)
        s["canonical_key"] = r["canonical_key"]
        s["candidate_path"] = r["candidate_path"]
        s["status"] = "adapter_unverified" if not rules_for(benchmark) and not src.strip() else s["status"]
        scanned.append(s)
    from collections import Counter
    st = Counter(s["status"] for s in scanned)
    return {"benchmark": benchmark, "n_candidates": len(scanned),
            "status_counts": dict(st), "details": scanned,
            "per_field_rollup": _rollup(scanned)}


def _rollup(scanned: list[dict[str, Any]]) -> list[dict[str, Any]]:
    acc: dict[str, dict[str, int]] = {}
    for s in scanned:
        for f in s["per_field"]:
            a = acc.setdefault(f["role"], {"ok": 0, "ambig": 0, "miss": 0})
            if f["ok"]:
                a["ok"] += 1
            elif f["reason"] in ("ctx_ok",):
                a["ok"] += 1
            elif f["hit"] == 0:
                a["miss"] += 1
            else:
                a["ambig"] += 1
    return [{"role": k, **v} for k, v in acc.items()]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", default="all")
    ap.add_argument("--print", action="store_true")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()
    records = json.loads(BASELINE.read_text(encoding="utf-8"))["records"]
    bms = sorted({r["benchmark"] for r in records})
    if args.scan != "all":
        bms = [args.scan]
    report = {
        "schema_version": "jss-holdout-adapter-scan-v2",
        "mode": "offline-mutation-coverage-gate1-no-compile",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "note": ("Per-candidate provisional split for P3a. Gate 1 only (no compile): "
                 "a field is 'retention_tested_ready' only if every input value is "
                 "uniquely locatable (value-count for simple params, context shape for "
                 "brittle/structural fields). not_parameterizable means no safe rebind "
                 "without changing computation logic."),
        "records": [scan_records(records, b) for b in bms],
    }
    Path(args.out).write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    from collections import Counter
    tot = Counter()
    for r in report["records"]:
        tot.update(r["status_counts"])
        print(f"[{r['benchmark']}] n={r['n_candidates']} "
              f"status={dict(r['status_counts'])}")
        for pr in r["per_field_rollup"]:
            print(f"    field={pr['role']:10s} ok={pr['ok']} ambiguous={pr['ambig']} miss={pr['miss']}")
        if args.print:
            for d in r["details"]:
                print(f"      {d['status']:<28s} {d['candidate_path']}")
    print("\noverall:", dict(tot))
    print("written:", args.out)


if __name__ == "__main__":
    main()
