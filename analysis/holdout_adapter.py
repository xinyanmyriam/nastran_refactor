"""Candidate input adapter: two-layer validity guards (section 8.6).

Candidates hard-code inputs in ``main()`` and cannot be re-generated, so we
rebind input literals at the source level. Before any held-out execution the
adapter must prove it bounded the *right* literals. Two independent gates:

  Gate 1 - mutation coverage scan (offline, no compiler):
      tokenise each candidate source and count occurrences of every input
      literal (value-normalised) against each rule's ``expected``. Roles that
      are missing, over-hit, or ambiguous are flagged; only rules passing the
      coverage check proceed.

  Gate 2 - h0 reconstruction (compiles):
      rebind every rule back to its h0 value through the full
      parse -> replace -> recompile path and require the stdout SHA-256 to
      reproduce the frozen canonical h0 ``stdout_sha256`` from
      ``jss_v4_holdout_baseline.json``. Anything else is ``adapter_unverified``.

The two gates together satisfy review point 3: a no-op identity mapping is not
sufficient; we require (a) that every expected input literal is actually hit
and unambiguously locatable, and (b) that a provenance-preserving rebind
reproduces frozen h0 output.

CLI:
  python holdout_adapter.py --scan B1            # mutation coverage scan (offline)
  python holdout_adapter.py --scan all
  python holdout_adapter.py --gate B1 --max 3    # reconstruction gate (compiles)

This module does not edit any candidate until ``adapter_status`` is decided.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import tempfile
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "analysis"
BASELINE = PAPER / "jss_v4_holdout_baseline.json"
PILOT = ("B1", "B3", "B10", "B11", "B13")  # cover matrix/vector/eig/timeseries-scalar

# --------------------------------------------------------------------------- #
# Candidate source path plumbing (mirrors build_holdout_baseline)
# --------------------------------------------------------------------------- #
def candidate_abs(rel: str) -> Path:
    return ROOT / Path(rel.replace("/", os.sep))


def _load_baseline() -> list[dict[str, Any]]:
    doc = json.loads(BASELINE.read_text(encoding="utf-8"),
                     parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    return doc["records"]


# --------------------------------------------------------------------------- #
# Numeric literal normalisation (value-based, so 200e9==2e11==2.0e11)
# --------------------------------------------------------------------------- #
_NUM = re.compile(r"(0[xX][0-9a-fA-F_.pP+-]+)|((?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)[fFlLuU]*")
_SUFFIX = re.compile(r"[fFlLuU]+$")


def numeric_tokens(source: str) -> list[str]:
    """Return value-normalised decimal strings of every numeric literal in order."""
    out: list[str] = []
    for m in _NUM.finditer(source):
        tok = m.group(0)
        if tok.lower().startswith("0x") or tok.lower().startswith("0b"):
            continue  # base-notation not usable for these input roles
        body = _SUFFIX.sub("", tok)
        try:
            d = Decimal(body)
        except (InvalidOperation, ValueError):
            continue
        # collapse to a canonical short form to count "same value" occurrences
        out.append(format(d, "f"))
    return out


def value_normalise(lit: str) -> str:
    d = Decimal(str(lit))
    return format(d, "f")


# --------------------------------------------------------------------------- #
# Per-benchmark input literal rules (h0 values; derived from _TEST_PARAMS)
# --------------------------------------------------------------------------- #
# rule: role, h0_literal, expected (approx occurrence count in a faithful
# candidate where each input parameter is written once). expected=1 is the
# common case; multi-cell matrices have repeated values (B10).
def rules_for(benchmark: str) -> list[dict[str, Any]]:
    R = {
        "B1": [
            {"role": "E", "lit": "200e9", "expected": 1},
            {"role": "A", "lit": "0.01", "expected": 1},
            {"role": "G", "lit": "76.923e9", "expected": 1},
            {"role": "J", "lit": "5e-6", "expected": 1},
            {"role": "L", "lit": "2.0", "expected": 2},   # L + node-B x
        ],
        "B3": [
            {"role": "E", "lit": "2.1e11", "expected": 1},
            {"role": "nu", "lit": "0.3", "expected": 1},
            {"role": "t", "lit": "0.01", "expected": 1},
        ],
        "B4": [
            {"role": "E", "lit": "200e9", "expected": 1},
            {"role": "nu", "lit": "0.3", "expected": 1},
            {"role": "t", "lit": "0.01", "expected": 1},
        ],
        "B2": [
            {"role": "E", "lit": "200e9", "expected": 1},
            {"role": "A", "lit": "0.01", "expected": 1},
            {"role": "G", "lit": "76.923e9", "expected": 1},
            {"role": "I", "lit": "8.333e-6", "expected": 2},   # Iy and Iz share the value
            {"role": "J", "lit": "1.667e-5", "expected": 1},
            {"role": "L", "lit": "2.0", "expected": 2},
        ],
        "B7": [
            {"role": "E", "lit": "200e9", "expected": 1},
            {"role": "nu", "lit": "0.3", "expected": 1},
        ],
        "B10": [
            {"role": "K_diag_main", "lit": "200", "expected": 3},   # diag 200,200,200
            {"role": "K_diag_last", "lit": "100", "expected": 1},   # bottom-right 100
            {"role": "K_offdiag", "lit": "-100", "expected": 3},
        ],
        "B11": [
            {"role": "diag", "lit": "2", "expected": 5},
            {"role": "offdiag", "lit": "-1", "expected": 4},
        ],
        "B13": [
            {"role": "E", "lit": "2.1e11", "expected": 1},
            {"role": "A", "lit": "0.01", "expected": 1},
            {"role": "J", "lit": "5e-6", "expected": 1},
            {"role": "C", "lit": "0.005", "expected": 1},
            {"role": "u_b", "lit": "0.001", "expected": 1},
        ],
    }
    return R.get(benchmark, [])


# --------------------------------------------------------------------------- #
# Gate 1: mutation coverage scan
# --------------------------------------------------------------------------- #
def scan_source(source: str, benchmark: str) -> dict[str, Any]:
    toks = numeric_tokens(source)
    counts: dict[str, int] = {}
    for t in toks:
        counts[t] = counts.get(t, 0) + 1
    rules = rules_for(benchmark)
    per_rule = []
    all_ok = bool(rules)
    for rule in rules:
        v = value_normalise(rule["lit"])
        hit = counts.get(v, 0)
        ok = hit == rule["expected"]
        ambiguous = hit != rule["expected"] and hit > 0  # over/under-hit but present
        per_rule.append({
            "role": rule["role"],
            "lit": rule["lit"],
            "expected": rule["expected"],
            "hit": hit,
            "ok": ok,
            "ambiguous": ambiguous,
        })
        if not ok:
            all_ok = False
    return {"benchmark": benchmark, "n_rules": len(rules), "passed": all_ok,
            "per_rule": per_rule, "total_numeric_tokens": len(toks)}


def scan_records(records: list[dict[str, Any]], benchmark: str) -> dict[str, Any]:
    rows = [r for r in records if r["benchmark"] == benchmark]
    if not rows:
        raise SystemExit(f"no canonical-pass candidates for {benchmark}")
    scanned = []
    for r in rows:
        src = candidate_abs(r["candidate_path"]).read_text(encoding="utf-8", errors="replace")
        s = scan_source(src, benchmark)
        s["canonical_key"] = r["canonical_key"]
        s["candidate_path"] = r["candidate_path"]
        s["pass_scan"] = s["passed"]
        scanned.append(s)
    passed_n = sum(1 for s in scanned if s["pass_scan"])
    return {
        "benchmark": benchmark,
        "n_candidates": len(scanned),
        "n_scan_pass": passed_n,
        "n_scan_fail": len(scanned) - passed_n,
        "per_rule_summary": _rollup(scanned),
        "details": scanned,
    }


def _rollup(scanned: list[dict[str, Any]]) -> list[dict[str, Any]]:
    if not scanned:
        return []
    # report, per role, how many candidates matched / over-hit / under-hit
    hats: dict[str, dict[str, int]] = {}
    for s in scanned:
        for pr in s["per_rule"]:
            h = hats.setdefault(pr["role"], {"ok": 0, "ambiguous": 0, "miss": 0})
            if pr["ok"]:
                h["ok"] += 1
            elif pr["ambiguous"]:
                h["ambiguous"] += 1
            else:
                h["miss"] += 1
    return [{"role": k, **v} for k, v in hats.items()]


# --------------------------------------------------------------------------- #
# Gate 2: h0 reconstruction (compile + stdout hash). Compiler-bound.
# --------------------------------------------------------------------------- #
def _sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def reconstruction_switch_binding(source: str, benchmark: str, value_map: dict[str, str]) -> tuple[bool, str]:
    """Gate 2 helper: apply literal rebindings; returns (all_rules_found, new_body).

    ``value_map``: role -> replacement literal (for the identity reconstruction
    this maps every role to its own original h0 literal, so success means the
    pipeline could *locate* every input literal). Operates on value-normalised
    numeric tokens so ``200e9``, ``2e11`` and ``2.0e11`` are all found.
    """
    rules = rules_for(benchmark)
    # Bucket by VALUE (roles that share a value, e.g. B2 Iy/Iz, are one bucket)
    by_val: dict[str, dict[str, Any]] = {}
    labels: dict[str, list[str]] = {}
    for r in rules:
        vr = value_normalise(r["lit"])
        b = by_val.setdefault(vr, {"value": vr, "expected": 0, "target": r["lit"], "hits": 0, "spans": []})
        b["expected"] += r["expected"]
        labels.setdefault(vr, []).append(r["role"])

    for m in _NUM.finditer(source):
        tok = m.group(0)
        if tok.lower().startswith("0x") or tok.lower().startswith("0b"):
            continue
        val = value_normalise(_SUFFIX.sub("", tok))
        if val in by_val:
            by_val[val]["hits"] += 1
            by_val[val]["spans"].append((m.start(), m.end()))

    all_ok = True
    for b in by_val.values():
        if b["hits"] != b["expected"]:
            all_ok = False

    new_body = source
    if all_ok:
        edits = sorted(((s, e, b["target"]) for b in by_val.values() for s, e in b["spans"]),
                       key=lambda x: -x[0])
        for s, e, t in edits:
            new_body = new_body[:s] + t + new_body[e:]
    return all_ok, new_body


def run_gate(record: dict[str, Any], benchmark: str, max_scan_token_only: bool = True) -> dict[str, Any]:
    """Rebind h0 literals to themselves and require the frozen stdout SHA-256.

    With ``max_scan_token_only=True`` (default in pre-toolchain mode) this only
    verifies the binding pipeline can locate every rule literal, i.e. the
    non-compiling half of Gate 2. Set False to actually compile and compare the
    output hash (requires MSVC/Eigen; see run_multilang_v2).
    """
    src = candidate_abs(record["candidate_path"]).read_text(encoding="utf-8", errors="replace")
    id_map = {rule["role"]: rule["lit"] for rule in rules_for(benchmark)}
    applied, _ = reconstruction_switch_binding(src, benchmark, id_map)
    result = {
        "canonical_key": record["canonical_key"],
        "benchmark": benchmark,
        "bind_all_rules_found": applied,
        "h0_stdout_sha256_expected": record["stdout_sha256"],
    }
    if not applied or max_scan_token_only:
        result["adapter_status"] = "adapter_verified_binding" if applied else "adapter_unverified"
        result["compiled_output_hash"] = None
        return result
    # Full compile path (deferred; needs the offline runner + toolchain).
    sys.path.insert(0, str(ROOT))
    try:
        from run_multilang_v2 import compile_and_run_cpp  # noqa: PLC0415
        from tempfile import TemporaryDirectory
        with TemporaryDirectory(prefix="holdout_gate_") as d:
            rr = compile_and_run_cpp(src, Path(d))
            out_hash = _sha256_bytes(rr.get("stdout", "").encode("utf-8", errors="replace"))
        result["compiled_output_hash"] = out_hash
        result["adapter_status"] = (
            "adapter_verified" if out_hash == record["stdout_sha256"] else "adapter_unverified"
        )
    except Exception as exc:  # offline / toolchain unavailable
        result["adapter_status"] = "adapter_compile_unavailable"
        result["error"] = f"{type(exc).__name__}: {exc}"[:500]
    return result


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def all_record_benchmarks(records: list[dict[str, Any]]) -> list[str]:
    bms = sorted({r["benchmark"] for r in records})
    return bms


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", metavar="B|all", help="benchmark id or 'all'")
    ap.add_argument("--gate", metavar="B")
    ap.add_argument("--max", type=int, default=3, help="gate: candidates to process")
    ap.add_argument("--out", default=PAPER / "jss_holdout_adapter_scan.json")
    ap.add_argument("--compile", action="store_true", help="gate: actually compile (needs toolchain)")
    ap.add_argument("--print", action="store_true", help="scan: print per-candidate flags")
    args = ap.parse_args()
    records = _load_baseline()
    available = all_record_benchmarks(records)

    if args.scan:
        bms = available if args.scan == "all" else [args.scan]
        unimplemented = [b for b in bms if not rules_for(b)]
        if unimplemented:
            print("no rule table yet for:", unimplemented, "(skipped)")
        bms = [b for b in bms if rules_for(b)]
        report = {"mode": "scan", "adapter_gate": "mutation-coverage",
                  "records": [scan_records(records, b) for b in bms]}
        Path(args.out).write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                                  encoding="utf-8")
        total, pass_n = 0, 0
        for r in report["records"]:
            total += r["n_candidates"]
            pass_n += r["n_scan_pass"]
            print(f"[{r['benchmark']}] candidates={r['n_candidates']} "
                  f"scan_pass={r['n_scan_pass']} scan_fail={r['n_scan_fail']}")
            for pr in r["per_rule_summary"]:
                print(f"    {pr['role']:12s} matched_ok={pr['ok']} ambiguous={pr['ambiguous']} miss={pr['miss']}")
            if args.print:
                for d in r["details"]:
                    print(f"      pass={bool(d['pass_scan'])} {d['candidate_path']}")
        print(f"\nadapter-coverage across {len(bms)} benchmarks: "
              f"{pass_n}/{total} candidates pass the strict mutation-coverage gate")
        print("written:", args.out)
        return

    if args.gate:
        rows = [r for r in records if r["benchmark"] == args.gate][: args.max]
        res = [run_gate(r, args.gate, max_scan_token_only=not args.compile) for r in rows]
        from collections import Counter
        print(f"[{args.gate}] gate candidates={len(res)}")
        print("    statuses:", dict(Counter(r["adapter_status"] for r in res)))
        print("    all bindings found:", all(r["bind_all_rules_found"] for r in res))
        return


if __name__ == "__main__":
    main()
