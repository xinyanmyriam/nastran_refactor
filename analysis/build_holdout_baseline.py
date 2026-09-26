"""Materialize the authoritative held-out retest population.

P0.5 (unified h0 baseline) of the held-out plan is met by the frozen canonical
rescore (`jss_v4_canonical_rescore.json`): every one of the 520 ledger
candidates was recompiled, executed and re-scored in the current frozen
toolchain, producing 173 equivalence-pass records with zero pass flips.

This script extracts exactly those 173 canonical-pass observations into a
frozen retest-population manifest (`jss_v4_holdout_baseline.json`) that the
adapter (P3a) and held-out execution (P3b) consume. Each record keeps the
canonical key, candidate/oracle hashes and the frozen h0 stdout SHA-256 used by
the h0 reconstruction gate.

Usage:
    E:\\...\\refactor\\.venv\\Scripts\\python.exe build_holdout_baseline.py
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "analysis"

RESCORE = PAPER / "jss_v4_canonical_rescore.json"
OUT_JSON = PAPER / "jss_v4_holdout_baseline.json"
EXPECTED_EQUIVALENT = 173


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def strict_value(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): strict_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [strict_value(v) for v in value]
    if hasattr(value, "item"):
        return strict_value(value.item())

    # keep numpy-safe JSON emission (no NaN/Infinity)
    try:
        import numbers
        if isinstance(value, numbers.Real) and value != value:
            return None
    except Exception:
        pass
    return value


def main() -> None:
    rescore = json.loads(RESCORE.read_text(encoding="utf-8"),
                         parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    obs = rescore.get("observations")
    if not isinstance(obs, list):
        raise SystemExit("rescore: no observations")
    passes = [r for r in obs if bool(r.get("equivalence"))]
    if len(passes) != EXPECTED_EQUIVALENT:
        raise SystemExit(f"expected {EXPECTED_EQUIVALENT} equivalence-pass records, got {len(passes)}")

    records: list[dict[str, Any]] = []
    for i, r in enumerate(sorted(passes, key=lambda x: (x["benchmark"], x["cohort"], x["repeat"])), 1):
        cand = ROOT / Path(r["candidate_path"].replace("/", os.sep))
        if not cand.is_file():
            raise SystemExit(f"candidate missing: {r['candidate_path']}")
        actual = sha256(cand)
        if actual != r["candidate_hash"]:
            raise SystemExit(f"candidate hash mismatch: {r['candidate_path']}")

        verification = r.get("verification") or {}
        records.append({
            "index": i,
            "canonical_key": r["canonical_key"],
            "cohort": r["cohort"],
            "model": r["model"],
            "method": r["method"],
            "benchmark": r["benchmark"],
            "experiment_id": r["experiment_id"],
            "repeat": r["repeat"],
            "candidate_path": r["candidate_path"],
            "candidate_hash": r["candidate_hash"],
            "oracle_path": r["oracle_path"],
            "oracle_hash": r["oracle_hash"],
            "stdout_sha256": r["stdout_sha256"],
            "h0_max_rel_diff": verification.get("max_rel_diff"),
        })

    def summarize(rowset: list[dict[str, Any]]) -> dict[str, Any]:
        by_cohort: dict[str, int] = {}
        by_benchmark: dict[str, int] = {}
        for rec in rowset:
            by_cohort[rec["cohort"]] = by_cohort.get(rec["cohort"], 0) + 1
            by_benchmark[rec["benchmark"]] = by_benchmark.get(rec["benchmark"], 0) + 1
        return {
            "total": len(rowset),
            "by_cohort": dict(sorted(by_cohort.items())),
            "by_benchmark": dict(sorted(by_benchmark.items())),
        }

    manifest = {
        "schema_version": "jss-v4-holdout-baseline-1",
        "built_from": "analysis/jss_v4_canonical_rescore.json",
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "note": (
            "Authoritative held-out retest population: candidates that passed the "
            "unified current-environment h0 rescore with zero pass flips. Adapter "
            "(P3a) and held-out execution (P3b) operate only on this population."
        ),
        "rescore_provenance": {
            "complete": bool(rescore.get("complete")),
            "ledger_path": rescore.get("ledger_path"),
            "ledger_hash": rescore.get("ledger_hash"),
            "mode": rescore.get("mode"),
        },
        "summary": summarize(records),
        "records": records,
    }
    text = json.dumps(strict_value(manifest), ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    tmp = OUT_JSON.with_suffix(".json.tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, OUT_JSON)
    manifest_hash = sha256(OUT_JSON)
    print(f"baseline written: {OUT_JSON}  ({len(records)} records, sha256={manifest_hash})")
    print(json.dumps(manifest["summary"], ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
