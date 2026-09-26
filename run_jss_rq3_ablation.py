"""Run the targeted JSS RQ3 tool ablation with resumable progress.

Formal design (default):
    5 configurations x B4/B7/B11 x 5 valid repetitions = 75 runs.

The runner randomizes configuration order within each repeat/benchmark block,
keeps infrastructure-invalid attempts out of the denominator, and writes an
incremental progress.json after every attempt.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import random
import sys
import time
from collections import defaultdict
from datetime import datetime
from pathlib import Path

from experiment_logger import ExperimentLogger, create_experiment_id
from config import MODELS
from run_experiments_v3 import (
    API_REQUEUE_LIMIT,
    PROJECT_ROOT,
    RQ3_METHOD_CONFIGS,
    _api_backoff,
    load_benchmark_config,
    run_single_experiment,
)

DEFAULT_BENCHMARKS = ("B4", "B7", "B11")
DEFAULT_METHODS = (
    "rq3_agent_full",
    "rq3_agent_no_oracle",
    "rq3_agent_feedback_only",
    "rq3_direct_r5",
    "rq3_direct_r15",
)
DEFAULT_MODEL = "deepseek-v4-flash"
DEFAULT_SEED = 20260818


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _redacted_model_config(model: str) -> dict:
    if model not in MODELS:
        raise SystemExit(f"Unknown model: {model}")
    cfg = MODELS[model]
    return {
        "provider": cfg.get("provider"),
        "model": cfg.get("model"),
        "base_url": cfg.get("base_url"),
        "max_tokens": cfg.get("max_tokens"),
    }


def _frozen_hashes(benchmarks: list[str], model: str) -> dict[str, str]:
    """Hash every executable protocol input except secret-bearing config.py."""
    paths = [
        PROJECT_ROOT / "run_experiments_v3.py",
        PROJECT_ROOT / "run_multilang_v2.py",
        PROJECT_ROOT / "tool_augmented_agent.py",
        PROJECT_ROOT / "llm_client.py",
        PROJECT_ROOT / "experiment_logger.py",
        PROJECT_ROOT / "benchmarks_config.json",
        Path(__file__).resolve(),
    ]
    for benchmark in benchmarks:
        cfg = load_benchmark_config(benchmark)
        for relative in cfg.get("source_files", []):
            paths.append(PROJECT_ROOT / relative)
        reference = cfg.get("reference_output_path")
        if reference:
            paths.append(PROJECT_ROOT / reference)
    unique = sorted({path.resolve() for path in paths}, key=str)
    missing = [str(path) for path in unique if not path.exists()]
    if missing:
        raise SystemExit(f"Missing frozen protocol inputs: {missing}")
    hashes = {
        path.relative_to(PROJECT_ROOT).as_posix(): _sha256(path)
        for path in unique
    }
    redacted = json.dumps(
        _redacted_model_config(model), sort_keys=True,
        separators=(",", ":")).encode("utf-8")
    hashes[f"redacted_model_config:{model}"] = hashlib.sha256(
        redacted).hexdigest()
    return hashes


def _acquire_lock(run_dir: Path):
    """Hold an OS-level exclusive campaign lock until this process exits."""
    path = run_dir / ".campaign.lock"
    handle = path.open("a+b")
    handle.seek(0, os.SEEK_END)
    if handle.tell() == 0:
        handle.write(b"0")
        handle.flush()
    handle.seek(0)
    try:
        if os.name == "nt":
            import msvcrt
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except (OSError, ImportError) as exc:
        handle.close()
        raise SystemExit(
            f"Campaign already active or lock unavailable: {path}") from exc
    return handle


def _release_lock(handle) -> None:
    handle.seek(0)
    if os.name == "nt":
        import msvcrt
        msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
    else:
        import fcntl
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    handle.close()


def _serializable_method_config(method: str) -> dict:
    cfg = RQ3_METHOD_CONFIGS[method]
    return {
        "base_method": cfg["base_method"],
        "enabled_tools": sorted(cfg["enabled_tools"]),
        "max_outer_repairs": cfg["max_outer_repairs"],
        "agent_max_steps": 15 if cfg["base_method"] == "agent" else None,
        "agent_max_estimated_tokens": (100_000
                                       if cfg["base_method"] == "agent"
                                       else None),
        "token_accounting": (
            "character_count_divided_by_3; pre-call gate; final-call "
            "overshoot possible"
            if cfg["base_method"] == "agent" else None),
        "knowledge_mode": "off" if cfg["base_method"] == "agent" else None,
        "checkpoints_enabled": False if cfg["base_method"] == "agent" else None,
    }


def _atomic_json(path: Path, data: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, default=str), encoding="utf-8")
    tmp.replace(path)


def _schedule(benchmarks: list[str], methods: list[str], repeats: int,
              seed: int) -> list[tuple[int, str, str]]:
    rng = random.Random(seed)
    schedule: list[tuple[int, str, str]] = []
    for repeat in range(1, repeats + 1):
        benchmark_order = list(benchmarks)
        rng.shuffle(benchmark_order)
        for benchmark in benchmark_order:
            method_order = list(methods)
            rng.shuffle(method_order)
            schedule.extend((repeat, benchmark, method)
                            for method in method_order)
    return schedule


def _summary(valid_results: list[dict], methods: list[str],
             benchmarks: list[str]) -> dict:
    cells: dict[str, dict] = {}
    for benchmark in benchmarks:
        for method in methods:
            rows = [r for r in valid_results
                    if r["benchmark"] == benchmark and r["method"] == method]
            cells[f"{benchmark}|{method}"] = {
                "valid": len(rows),
                "pass": sum(bool(r.get("passed")) for r in rows),
                "compile": sum(bool(r.get("compile_ok")) for r in rows),
                "run": sum(bool(r.get("run_ok")) for r in rows),
                "numerical": sum(bool(r.get("numerical_ok")) for r in rows),
            }
    return cells


def _validate_resumed_results(rows: list[dict], protocol: dict) -> None:
    """Fail closed on malformed, duplicated, or cross-protocol observations."""
    seen: set[tuple[str, str, int]] = set()
    for index, row in enumerate(rows):
        if row.get("valid") is not True:
            raise SystemExit(f"Resume valid_results[{index}] is not valid=true")
        benchmark = row.get("benchmark")
        method = row.get("method")
        model = row.get("model")
        try:
            repeat = int(row.get("repeat"))
        except (TypeError, ValueError) as exc:
            raise SystemExit(
                f"Resume valid_results[{index}] has invalid repeat") from exc
        if benchmark not in protocol["benchmarks"]:
            raise SystemExit(f"Resume row has out-of-protocol benchmark: {benchmark}")
        if method not in protocol["methods"]:
            raise SystemExit(f"Resume row has out-of-protocol method: {method}")
        if model != protocol["model"]:
            raise SystemExit(f"Resume row has out-of-protocol model: {model}")
        if not 1 <= repeat <= protocol["repeats"]:
            raise SystemExit(f"Resume row has out-of-range repeat: {repeat}")
        if row.get("run_tag") != protocol["run_tag"]:
            raise SystemExit("Resume row run_tag does not match campaign")
        if row.get("rq3_method_config") != protocol["method_configs"][method]:
            raise SystemExit("Resume row method config does not match campaign")
        key = (benchmark, method, repeat)
        if key in seen:
            raise SystemExit(f"Duplicate valid resume cell: {key}")
        seen.add(key)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Resumable targeted JSS RQ3 tool ablation")
    parser.add_argument("--benchmarks", default=",".join(DEFAULT_BENCHMARKS))
    parser.add_argument("--methods", default=",".join(DEFAULT_METHODS))
    parser.add_argument("--model", default=DEFAULT_MODEL)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--run-tag", default="jss_rq3_n5")
    parser.add_argument("--output-root", default="results/jss_rq3_tool_ablation")
    parser.add_argument(
        "--campaign-role", default="targeted_rq3",
        help="machine-readable study role stored in the frozen protocol")
    parser.add_argument(
        "--selection",
        default="purposive mechanism probes; not a population estimate",
        help="predeclared sampling/selection statement stored verbatim")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate and print the randomized schedule only")
    args = parser.parse_args()

    benchmarks = [x.strip() for x in args.benchmarks.split(",") if x.strip()]
    methods = [x.strip() for x in args.methods.split(",") if x.strip()]
    unknown = [m for m in methods if m not in RQ3_METHOD_CONFIGS]
    if unknown:
        raise SystemExit(f"Unknown RQ3 methods: {unknown}")
    for benchmark in benchmarks:
        load_benchmark_config(benchmark)
    if args.repeats < 1:
        raise SystemExit("--repeats must be positive")

    root = Path(args.output_root)
    if not root.is_absolute():
        root = PROJECT_ROOT / root
    run_dir = root / args.run_tag
    run_dir.mkdir(parents=True, exist_ok=True)
    campaign_lock = _acquire_lock(run_dir)
    progress_file = run_dir / "progress.json"
    frozen_hashes = _frozen_hashes(benchmarks, args.model)
    redacted_model = _redacted_model_config(args.model)

    protocol = {
        "run_tag": args.run_tag,
        "created_at": datetime.now().isoformat(),
        "campaign_role": args.campaign_role,
        "model": args.model,
        "benchmarks": benchmarks,
        "methods": methods,
        "repeats": args.repeats,
        "seed": args.seed,
        "target_valid_runs": len(benchmarks) * len(methods) * args.repeats,
        "method_configs": {m: _serializable_method_config(m) for m in methods},
        "selection": args.selection,
        "redacted_model_config": redacted_model,
        "runtime_environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "os_name": os.name,
        },
        "frozen_sha256": frozen_hashes,
    }

    attempts: list[dict] = []
    valid_results: list[dict] = []
    if progress_file.exists():
        if not args.resume:
            _release_lock(campaign_lock)
            raise SystemExit(
                f"{progress_file} exists; use --resume or a different --run-tag")
        old = json.loads(progress_file.read_text(encoding="utf-8"))
        attempts = old.get("attempts", [])
        valid_results = old.get("valid_results", [])
        old_protocol = old.get("protocol", {})
        for field in (
            "run_tag", "campaign_role", "model", "benchmarks", "methods",
            "repeats", "seed", "target_valid_runs", "method_configs",
            "selection", "redacted_model_config", "runtime_environment",
            "frozen_sha256",
        ):
            if old_protocol.get(field) != protocol.get(field):
                _release_lock(campaign_lock)
                raise SystemExit(f"Resume protocol mismatch for {field}")
        _validate_resumed_results(valid_results, old_protocol)
        protocol["created_at"] = old_protocol.get("created_at",
                                                   protocol["created_at"])

    completed = {
        (r["benchmark"], r["method"], int(r["repeat"]))
        for r in valid_results if r.get("valid", True)
    }
    schedule = _schedule(benchmarks, methods, args.repeats, args.seed)
    total = len(schedule)

    print("=" * 72)
    print("JSS RESUMABLE FORMAL CAMPAIGN")
    print("=" * 72)
    print(f"campaign_role: {args.campaign_role}")
    print(f"run_dir: {run_dir}")
    print(f"model: {args.model}")
    print(f"benchmarks: {benchmarks}")
    print(f"methods: {methods}")
    print(f"valid repeats/cell: {args.repeats}")
    print(f"target valid runs: {total}; already complete: {len(completed)}")
    if args.dry_run:
        for sequence, (repeat, benchmark, method) in enumerate(schedule,
                                                               start=1):
            print(f"{sequence:03d}: r{repeat} {benchmark} {method}")
        _release_lock(campaign_lock)
        return

    for sequence, (repeat, benchmark, method) in enumerate(schedule, start=1):
        key = (benchmark, method, repeat)
        if key in completed:
            continue

        cfg = _serializable_method_config(method)
        infra_attempt = 0
        result = None
        while infra_attempt <= API_REQUEUE_LIMIT:
            infra_attempt += 1
            eid = create_experiment_id(
                f"rq3_{args.run_tag}_{method}_r{repeat}_a{infra_attempt}",
                benchmark)
            logger = ExperimentLogger(eid, base_dir=str(run_dir / "runs"))
            logger.log_artifact("rq3_protocol", {
                "campaign_role": args.campaign_role,
                "run_tag": args.run_tag,
                "selection": args.selection,
                "frozen_sha256": frozen_hashes,
                "redacted_model_config": redacted_model,
                "benchmark": load_benchmark_config(benchmark),
                "method": method,
                "method_config": cfg,
                "repeat": repeat,
                "infra_attempt": infra_attempt,
                "sequence": sequence,
                "seed": args.seed,
                "formal": "smoke" not in args.run_tag.lower(),
            })

            print(f"[{sequence}/{total}] {benchmark} {method} "
                  f"r{repeat} attempt {infra_attempt}")
            result = run_single_experiment(
                benchmark, args.model, method, repeat, logger,
                max_repair_iterations=cfg["max_outer_repairs"])
            result.update({
                "sequence": sequence,
                "infra_attempt": infra_attempt,
                "run_tag": args.run_tag,
                "rq3_method_config": cfg,
                "experiment_id": eid,
            })
            attempts.append(result)
            logger.finalize(
                bool(result.get("valid", True) and result.get("passed")),
                f"{benchmark} {method} r{repeat}: "
                f"{'PASS' if result.get('passed') else 'FAIL'} "
                f"valid={result.get('valid', True)}")

            if result.get("valid", True):
                valid_results.append(result)
                completed.add(key)
                break
            if result.get("abort_reason") != "api_exhausted":
                break
            if infra_attempt <= API_REQUEUE_LIMIT:
                _api_backoff(infra_attempt + 2, logger)

            progress = {
                "protocol": protocol,
                "updated_at": datetime.now().isoformat(),
                "attempts": attempts,
                "valid_results": valid_results,
                "summary": _summary(valid_results, methods, benchmarks),
            }
            _atomic_json(progress_file, progress)

        progress = {
            "protocol": protocol,
            "updated_at": datetime.now().isoformat(),
            "attempts": attempts,
            "valid_results": valid_results,
            "summary": _summary(valid_results, methods, benchmarks),
        }
        _atomic_json(progress_file, progress)

        if result is not None:
            status = ("INVALID" if not result.get("valid", True)
                      else "PASS" if result.get("passed") else "FAIL")
            print(f"  -> {status}; completed valid {len(completed)}/{total}")

    missing = [
        {"repeat": repeat, "benchmark": benchmark, "method": method}
        for repeat, benchmark, method in schedule
        if (benchmark, method, repeat) not in completed
    ]
    final = {
        "protocol": protocol,
        "updated_at": datetime.now().isoformat(),
        "complete": not missing,
        "missing": missing,
        "attempts": attempts,
        "valid_results": valid_results,
        "summary": _summary(valid_results, methods, benchmarks),
    }
    _atomic_json(progress_file, final)
    _atomic_json(run_dir / "final_summary.json", final)

    print("=" * 72)
    print(f"complete: {not missing}; valid results: {len(valid_results)}/{total}")
    for cell, values in final["summary"].items():
        print(f"{cell}: {values['pass']}/{values['valid']}")
    if missing:
        print(f"missing cells: {len(missing)}; resume with the same command + --resume")
    _release_lock(campaign_lock)


if __name__ == "__main__":
    main()
