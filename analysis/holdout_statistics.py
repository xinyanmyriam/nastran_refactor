"""Held-out retention statistics: explicit denominators + benchmark clustering.

Implements the two review-critical fixes for the held-out plan (sections 8.8):

  * EXplicit denominators (review point 2): a held-out result is only
    measurable when the adapter could be validated. We report three coherent
    quantities and never let ``adapter_unverified`` / ``not_parameterizable``
    silently disappear from the denominator.

      - Conservative retention: denominator = all canonical-h0-pass candidates;
        unverified / not_parameterizable count as *not* retained.
      - Conditional retention: denominator = only candidates whose adapter was
        fully validated (``retention_tested``).
      - Identification bounds: lower = unverifiable all-fail; upper =
        unverifiable all-pass.

  * Benchmark clustering (review point 4): candidates are nested in 13
    benchmarks and the same candidate's h1-h3 are repeated measures, so pooled
    candidates are not independent. The main metric is the 13-benchmark
    macro-average; uncertainty is a benchmark-level bootstrap (resamples the
    13 benchmark blocks with replacement, each block kept intact).

Input data model (``holdout_results.json``)::

    {
      "candidates": [
        {
          "cohort": "...", "method": "...", "benchmark": "B1",
          "experiment_id": "...", "repeat": 1,
          "adapter_status": "retention_tested" | "adapter_unverified" | "not_parameterizable",
          "adapter_note": "...",
          "inputs": {
            "h1": {"pass": true, ...}, "h2": ..., "h3": ...
          }
        }, ...
      ]
    }

Usage:
    python holdout_statistics.py --selftest     # verify on synthetic data
    python holdout_statistics.py --report       # read holdout_results.json if present
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Callable

ROOT = Path(__file__).resolve().parents[1]
IN_JSON = ROOT / "analysis" / "holdout_results.json"
N_INPUTS = 3
RNG_SEED = 20260921

OK = "retention_tested"
UNVERIFIED = "adapter_unverified"
NOT_PARAM = "not_parameterizable"
VALID_STATUS = {OK, UNVERIFIED, NOT_PARAM}


# --------------------------------------------------------------------------- #
# Per-candidate classification
# --------------------------------------------------------------------------- #
def full_retention(candidate: dict[str, Any], n_inputs: int = N_INPUTS) -> bool | None:
    """True if the candidate is adapter-validated AND passes all held-out inputs.

    Returns None when the candidate is not measurable (adapter not validated),
    signalling it must be handled by explicit denominator rules rather than
    silently dropped.
    """
    status = candidate.get("adapter_status")
    if status not in VALID_STATUS:
        return None
    if status != OK:
        return None
    inputs = candidate.get("inputs") or {}
    if len([k for k in ("h1", "h2", "h3") if k in inputs]) < n_inputs:
        return None  # missing inputs -> not measurable
    return all(bool(inputs[k].get("pass")) for k in ("h1", "h2", "h3"))


# --------------------------------------------------------------------------- #
# Denominator metrics (review point 2)
# --------------------------------------------------------------------------- #
def _metric_block(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    total = len(candidates)
    n_ok = sum(1 for c in candidates if c.get("adapter_status") == OK)
    n_unverified = sum(1 for c in candidates if c.get("adapter_status") == UNVERIFIED)
    n_notparam = sum(1 for c in candidates if c.get("adapter_status") == NOT_PARAM)

    retained = sum(1 for c in candidates if full_retention(c) is True)
    # measurable = adapter-validated; conditional denominator
    if n_ok == 0:
        conditional = None
    else:
        conditional = retained / n_ok

    # Conservative: unverifiable count as failure
    conservative = retained / total if total else None

    # Upper bound: unverifiable count as pass
    upper_retained = retained + n_unverified + n_notparam
    upper = upper_retained / total if total else None

    return {
        "n": total,
        "n_ok": n_ok,
        "n_unverified": n_unverified,
        "n_not_parameterizable": n_notparam,
        "n_retained": retained,
        "conservative": conservative,
        "conditional": conditional,
        "bounds": {"lower": conservative, "upper": upper},
        "adapter_coverage": (n_ok / total) if total else None,
    }


def aggregate(candidates: list[dict[str, Any]],
              group_by: Callable[[dict[str, Any]], Any],
              include_missing_groups: list[Any] | None = None) -> dict[str, Any]:
    groups: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for c in candidates:
        groups[group_by(c)].append(c)
    out: dict[Any, dict[str, Any]] = {}
    for key in sorted(set(groups) | set(include_missing_groups or [])):
        out[str(key)] = _metric_block(groups.get(key, []))
    return out


def macro_average(blocks: dict[str, Any], metric: str) -> dict[str, Any]:
    """Mean of per-benchmark conservative/conditional over non-empty benchmarks."""
    values = [b.get(metric) for b in blocks.values() if b.get("n", 0) > 0]
    present = [b for b in blocks.values() if b.get("n", 0) > 0]
    present_metric = [b.get(metric) for b in present]
    return {
        "metric": metric,
        "n_benchmarks": len(present),
        "macro_mean": (sum(present_metric) / len(present_metric)) if present_metric else None,
        "per_benchmark": {k: b.get(metric) for k, b in blocks.items() if b.get("n", 0) > 0},
        # 13-benchmark domain: report how many of the 13 are empty (zero-pass cells)
        "empty_cells_among_13": 13 - len(present),
        "all_values": values,
    }


# --------------------------------------------------------------------------- #
# Benchmark-level bootstrap (review point 4)
# --------------------------------------------------------------------------- #
def benchmark_bootstrap(candidates: list[dict[str, Any]],
                        metric: str = "conservative",
                        n_boot: int = 2000,
                        seed: int = RNG_SEED) -> dict[str, Any]:
    """Resample the 13 benchmark *blocks* (not candidates) with replacement.

    Each block keeps all its candidates intact; the macro-average is recomputed
    over the resampled set of (non-empty) benchmark blocks. Returns a
    percentile 95% CI around the observed macro-average.
    """
    rng = random.Random(seed)
    by_bm: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for c in candidates:
        by_bm[c["benchmark"]].append(c)
    block_names = sorted(by_bm)

    def macro(rows: list[dict[str, Any]]) -> float | None:
        blocks = aggregate(rows, lambda c: c["benchmark"])
        return macro_average(blocks, metric)["macro_mean"]

    observed = macro(candidates)
    if observed is None or not block_names:
        return {"metric": metric, "observed": None, "ci95": None,
                "note": "no data or only empty benchmark cells"}

    boots: list[float] = []
    for _ in range(n_boot):
        sampled: list[dict[str, Any]] = []
        for _ in range(len(block_names)):
            name = rng.choice(block_names)
            sampled.extend(by_bm[name])
        m = macro(sampled)
        if m is not None:
            boots.append(m)
    boots.sort()
    lo = boots[int(0.025 * len(boots))]
    hi = boots[int(0.975 * len(boots))]
    return {
        "metric": metric,
        "n_benchmarks": len(block_names),
        "n_boot": n_boot,
        "observed": observed,
        "ci95": [lo, hi],
        "seed": seed,
        "note": "95% percentile CI from benchmark-block resampling",
    }


# --------------------------------------------------------------------------- #
# Whole report
# --------------------------------------------------------------------------- #
def build_report(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    overall = aggregate(candidates, lambda c: "__all__")
    by_benchmark = aggregate(candidates, lambda c: c["benchmark"],
                             include_missing_groups=[f"B{i}" for i in range(1, 14)])
    by_cohort = aggregate(candidates, lambda c: c["cohort"])
    macro_conservative = macro_average(by_benchmark, "conservative")
    macro_conditional = macro_average(by_benchmark, "conditional")
    boot_conservative = benchmark_bootstrap(candidates, "conservative")
    boot_conditional = benchmark_bootstrap(candidates, "conditional")
    return {
        "schema_version": "jss-holdout-stats-1",
        "overall": overall["__all__"],
        "by_benchmark": by_benchmark,
        "by_cohort": by_cohort,
        "macro_conservative": macro_conservative,
        "macro_conditional": macro_conditional,
        "bootstrap_conservative": boot_conservative,
        "bootstrap_conditional": boot_conditional,
    }


def render(report: dict[str, Any]) -> str:
    lines = ["# Held-out retention statistics"]
    o = report["overall"]
    lines.append(f"\n## Overall (n={o['n']})")
    lines.append(f"- adapter-validated (retention_tested): {o['n_ok']}")
    lines.append(f"- adapter_unverified: {o['n_unverified']}")
    lines.append(f"- not_parameterizable: {o['n_not_parameterizable']}")
    lines.append(f"- retained on all h1-h3: {o['n_retained']}")
    lines.append(f"- conservative retention: {o['conservative']}")
    lines.append(f"- conditional retention: {o['conditional']}")
    lines.append(f"- bounds [lower, upper]: {o['bounds']}")

    lines.append("\n## By benchmark")
    lines.append("| B | n | tested | unverified | notparam | retained | conservative | conditional | bounds |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    empty = []
    for bm in [f"B{i}" for i in range(1, 14)]:
        b = report["by_benchmark"].get(bm)
        if not b or b["n"] == 0:
            empty.append(bm)
            lines.append(f"| {bm} | 0 (empty cell) | | | | | | | |")
            continue
        lo, hi = (b["bounds"]["lower"], b["bounds"]["upper"])
        lines.append(
            f"| {bm} | {b['n']} | {b['n_ok']} | {b['n_unverified']} | "
            f"{b['n_not_parameterizable']} | {b['n_retained']} | "
            f"{b['conservative']:.3f} | "
            f"{('%.3f' % b['conditional']) if b['conditional'] is not None else 'n/a'} | "
            f"[{lo:.3f},{hi:.3f}] |"
        )
    lines.append(f"- zero-pass cells (among 13): {empty}")

    lines.append("\n## Macro-average (benchmark-clustered)")
    mc = report["macro_conservative"]
    lines.append(f"- conservative macro-average = {mc['macro_mean']:.3f} over {mc['n_benchmarks']} benchmarks")
    lines.append(f"- bootstrap 95% CI = {report['bootstrap_conservative']['ci95']}")
    lines.append(f"- conditional macro-average = {report['macro_conditional']['macro_mean']}")
    lines.append(f"- bootstrap 95% CI = {report['bootstrap_conditional']['ci95']}")
    return "\n".join(lines)


# --------------------------------------------------------------------------- #
# Self-test on synthetic data with known answers
# --------------------------------------------------------------------------- #
def selftest() -> None:
    def mk(benchmark, status, passes: list[bool] | None, cohort="C"):
        c = {"cohort": cohort, "method": "m", "benchmark": benchmark,
             "experiment_id": "x", "repeat": 1, "adapter_status": status}
        if passes is None:
            c["inputs"] = {}
        else:
            c["inputs"] = {k: {"pass": p} for k, p in zip(("h1", "h2", "h3"), passes)}
        return c

    # Benchmark A: 4 tested (3 full-retain), 0 unverified
    # Benchmark B: 1 tested full-retain, 2 unverified
    data = [
        mk("A", OK, [True, True, True]),
        mk("A", OK, [True, True, True]),
        mk("A", OK, [True, True, True]),
        mk("A", OK, [True, False, True]),
        mk("B", OK, [True, True, True]),
        mk("B", UNVERIFIED, None),
        mk("B", NOT_PARAM, None),
    ]
    expect = {
        "A_conservative": 3 / 4,   # 3 retained / 4 total
        "A_conditional": 3 / 4,    # 3 retained / 4 tested
        "A_nretained": 3,
        "B_conservative": 1 / 3,   # 1 / 3 (unverified+notparam count fail)
        "B_conditional": 1,        # 1 / 1 tested
        "B_bounds": (1 / 3, 1.0),  # lower all-fail unverifiable, upper all-pass
    }
    rep = build_report(data)
    A = rep["by_benchmark"]["A"]
    B = rep["by_benchmark"]["B"]
    checks = {
        "A_conservative": A["conservative"],
        "A_conditional": A["conditional"],
        "A_nretained": A["n_retained"],
        "B_conservative": B["conservative"],
        "B_conditional": B["conditional"],
        "B_bounds": (B["bounds"]["lower"], B["bounds"]["upper"]),
    }
    ok = True
    for k, want in expect.items():
        got = checks[k]
        match = got == want
        ok = ok and match
        print(f"[selftest] {k:16s} expected={want!s:20s} got={got!s:20s} {'ok' if match else 'FAIL'}")
    # macro-average A,B conservative = (0.75 + 0.3333)/2
    mac = rep["macro_conservative"]["macro_mean"]
    want_mac = (0.75 + 1 / 3) / 2
    print(f"[selftest] macro_mean expected={want_mac:.4f} got={mac:.4f} {'ok' if abs(mac - want_mac) < 1e-9 else 'FAIL'}")
    ok = ok and abs(mac - want_mac) < 1e-9
    boot = benchmark_bootstrap(data, "conservative", n_boot=500)
    assert boot["ci95"], "bootstrap returned no CI"
    print(f"[selftest] benchmark bootstrap CI = {boot['ci95']}")
    print("[selftest] ALL PASS" if ok else "[selftest] SELF-TEST FAILED")
    if not ok:
        sys.exit(1)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--out", default=ROOT / "analysis" / "jss_holdout_retention.json")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not IN_JSON.exists():
        print(f"held-out results not present yet: {IN_JSON}")
        print("run `--selftest` to verify the statistics machinery on synthetic data.")
        return
    data = json.loads(IN_JSON.read_text(encoding="utf-8"),
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    rep = build_report(data["candidates"])
    Path(args.out).write_text(
        json.dumps(rep, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    print(render(rep))
    print(f"\nwritten: {args.out}")


if __name__ == "__main__":
    main()
