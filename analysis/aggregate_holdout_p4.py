# -*- coding: utf-8 -*-
"""P4: aggregate P3b holdout results into conservative/conditional/bounds retention.

Reads:
  - holdout_results.json                    (P3b per candidate x variant outcomes)
  - jss_holdout_adapter_scan_v2.json        (per-benchmark candidate totals +三分)
Writes:
  - aggregate_holdout_p4.json               (per-benchmark + macro/micro metrics)

Identure per 规划 8.8:
  - Conservative: kept / total_canonical_candidates   (not_parameterizable/unverified = not kept)
  - Conditional : kept / retention_tested
  - bounds      : lower = conservative; upper = (kept + retented-unverified extra) / total,
                  where untested candidates are counted as kept.
A candidate is "kept" iff it is retention_tested AND passes on ALL of h1,h2,h3.
Benchmark without candidates are reported as "no_candidates" and excluded from macro means.
"""
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

PAPER = Path(__file__).resolve().parent
OUT_JSON = PAPER / "holdout_results.json"
SCAN = PAPER / "jss_holdout_adapter_scan_v2.json"
MERGED_OUT = PAPER / "aggregate_holdout_p4.json"

VARIANTS = ["h1", "h2", "h3"]


def load_results():
    return json.load(open(OUT_JSON, encoding="utf-8"))["results"]


def load_scan():
    """Return per-benchmark {total, retention_tested, unverified, not_parameterizable}."""
    scan = json.load(open(SCAN, encoding="utf-8"))
    agg = defaultdict(lambda: {"total": 0, "retention_tested": 0,
                               "adapter_unverified": 0, "not_parameterizable": 0})
    for grp in scan["records"]:
        for c in grp.get("details", []):
            b = c["benchmark"]
            agg[b]["total"] += 1
            st = c.get("final_status")
            if st in agg[b]:
                agg[b][st] += 1
            else:
                agg[b]["retention_tested"]  # noop guard
    return agg


def full_kept_by_benchmark(results):
    """candidate kept iff retention candidate passes on all 3 variants."""
    by_cand = defaultdict(dict)
    for r in results:
        by_cand[(r["benchmark"], r["candidate_path"])][r["variant"]] = r["outcome"]
    kept = defaultdict(int)
    ret_candidates = defaultdict(int)
    for (b, _cp), per_v in by_cand.items():
        ret_candidates[b] += 1
        if all(per_v.get(v) == "pass" for v in VARIANTS):
            kept[b] += 1
    return kept, ret_candidates


def main():
    results = load_results()
    scan = load_scan()
    kept, ret_candidates = full_kept_by_benchmark(results)

    bench_rows = []
    for b in sorted(scan):
        t = scan[b]
        tot, ret = t["total"], t["retention_tested"]
        if tot == 0:
            bench_rows.append({
                "benchmark": b, "no_candidates": True,
            })
            continue
        k = kept.get(b, 0)
        unv = t["adapter_unverified"]
        npara = t["not_parameterizable"]
        coverage = ret / (ret + unv + npara) if (ret + unv + npara) else None
        conserv = k / tot
        cond = k / ret if ret else None
        upper = (k + (tot - ret)) / tot
        bench_rows.append({
            "benchmark": b,
            "total": tot,
            "retention_tested": ret,
            "adapter_unverified": unv,
            "not_parameterizable": npara,
            "kept_full": k,
            "adapter_coverage": round(coverage, 4) if coverage is not None else None,
            "conservative": round(conserv, 4),
            "conditional": round(cond, 4) if cond is not None else None,
            "upper_bound": round(upper, 4),
            "lower_bound": round(conserv, 4),
        })

    # macro-average over benchmarks with candidates
    active = [r for r in bench_rows if not r.get("no_candidates")]
    def macro(key):
        vals = [r[key] for r in active if r.get(key) is not None]
        return round(sum(vals) / len(vals), 4) if vals else None
    macro_conserv = macro("conservative")
    macro_cond = macro("conditional")
    macro_upper = macro("upper_bound")

    # micro (pooled) aggregates
    T = sum(r["total"] for r in active)
    K = sum(r["kept_full"] for r in active)
    R = sum(r["retention_tested"] for r in active)
    U = sum(r["adapter_unverified"] for r in active)
    NP = sum(r["not_parameterizable"] for r in active)
    micro = {
        "total": T, "retention_tested": R, "kept_full": K,
        "adapter_coverage": round(R / (R + U + NP), 4),
        "conservative": round(K / T, 4),
        "conditional": round(K / R, 4),
        "upper_bound": round((K + (T - R)) / T, 4),
        "lower_bound": round(K / T, 4),
    }

    # benchmark-level bootstrap for macro conservative (sampling unit = benchmark)
    import random
    random.seed(42)
    B = 10000
    means = []
    for _ in range(B):
        draws = [active[random.randint(0, len(active) - 1)]["conservative"]
                 for __ in range(len(active))]
        means.append(sum(draws) / len(draws))
    means.sort()
    lo = round(means[int(0.025 * B)], 4)
    hi = round(means[int(0.975 * B)], 4)

    doc = {
        "schema": "jss-holdout-p4-1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "note": (
            "P4 retention aggregation. Conservative = kept_full / total canonical candidates "
            "(not_parameterizable/unverified not kept). Conditional = kept_full / retention_tested. "
            "bounds: lower = conservative; upper = (kept + untested) / total. "
            "candidate kept = retention_tested AND passes h1 AND h2 AND h3. "
            "Macro = mean over benchmarks with candidates; CI via benchmark-level bootstrap."
        ),
        "variants": VARIANTS,
        "per_benchmark": bench_rows,
        "macro_conservative": macro_conserv,
        "macro_conservative_bootstrap95ci": [lo, hi],
        "macro_conditional": macro_cond,
        "macro_upper_bound": macro_upper,
        "micro_pooled": micro,
    }
    MERGED_OUT.write_text(json.dumps(doc, ensure_ascii=False, indent=2), encoding="utf-8")

    print("=== per-benchmark ===")
    hdr = f"{'bench':5} {'tot':>4} {'ret':>4} {'kept':>4} {'consrv':>6} {'cond':>6} {'upper':>6}"
    print(hdr)
    for r in bench_rows:
        if r.get("no_candidates"):
            print(f"{r['benchmark']:5}  no_candidates")
            continue
        print(f"{r['benchmark']:5} {r['total']:>4} {r['retention_tested']:>4} {r['kept_full']:>4} "
              f"{r['conservative']:>6} {str(r['conditional']):>6} {r['upper_bound']:>6}")
    print("=== macro (benchmark-averaged) ===")
    print("conservative:", macro_conserv, "bootstrap95ci:", [lo, hi])
    print("conditional :", macro_cond)
    print("upper bound :", macro_upper)
    print("=== micro (pooled) ===")
    print(micro)


if __name__ == "__main__":
    main()