"""
Unified Experiment Pipeline v3 - Config-Driven 13-Benchmark Runner.
Supports: matrix symmetry, vector, eigenvalue, timeseries, scalar verification.
Usage: python run_experiments_v3.py --benchmark all --model deepseek-v4-flash --repeat 3

v3.1 additions:
  --method multiagent_v2  : Multi-agent with source annotation mining
  --method multiagent_v3  : v2 plus three measured fixes - drop the dead Python
                            code_generation/equivalence stages, widen the
                            pipeline-to-generator channel beyond 2500 chars,
                            and stop seeding the pipeline with the hand-written
                            description. v2 -> v3 isolates these for ablation.
  --method multiagent_kb  : Multi-agent with source mining + knowledge base
  --seed-kb               : Seed knowledge base with Level 3 insights before running
  --clear-kb              : Clear knowledge base before running
"""

import argparse, json, os, re, sys, tempfile, textwrap
from collections import Counter
from datetime import datetime
from itertools import groupby
from pathlib import Path
from typing import Optional
import numpy as np

PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))
from config import MODELS
from experiment_logger import ExperimentLogger, create_experiment_id
from llm_client import LLMClient
from agents.orchestrator import RefactoringOrchestrator
from knowledge_base import KnowledgeBase
# Reuse compile/run infrastructure from v2
from run_multilang_v2 import (
    compile_and_run_cpp, compute_reference as compute_analytical_reference,
    parse_matrix_from_output, find_windows_sdk,
)

CONFIG_FILE = PROJECT_ROOT / "benchmarks_config.json"
MAX_REPAIR_ITERATIONS = 5

# Confirmatory RQ3 tool-ablation profiles. These labels are deliberately
# separate from the historical multiagent_* capability axes: the JSS main
# comparison uses the autonomous tool-augmented Agent, not the fixed staged
# orchestrator. Agent cells receive no outer repair after their 15-step loop.
RQ3_METHOD_CONFIGS = {
    "rq3_agent_full": {
        "base_method": "agent",
        "enabled_tools": {"compile_and_run", "compare_with_reference",
                          "read_source", "inspect_conventions"},
        "max_outer_repairs": 0,
    },
    "rq3_agent_no_oracle": {
        "base_method": "agent",
        "enabled_tools": {"compile_and_run", "read_source",
                          "inspect_conventions"},
        "max_outer_repairs": 0,
    },
    "rq3_agent_feedback_only": {
        "base_method": "agent",
        "enabled_tools": {"compile_and_run", "compare_with_reference"},
        "max_outer_repairs": 0,
    },
    "rq3_direct_r5": {
        "base_method": "direct",
        "enabled_tools": set(),
        "max_outer_repairs": 5,
    },
    "rq3_direct_r15": {
        "base_method": "direct",
        "enabled_tools": set(),
        # One initial generation plus at most 14 repairs = 15 LLM calls.
        "max_outer_repairs": 14,
    },
}


def rq3_method_config(method: str) -> Optional[dict]:
    """Return a copy-safe RQ3 method definition, if this is an RQ3 label."""
    cfg = RQ3_METHOD_CONFIGS.get(method)
    if cfg is None:
        return None
    return {**cfg, "enabled_tools": set(cfg["enabled_tools"])}

# ============================================================
# Multi-agent capability axes
# ============================================================
#
# The multi-agent arm used to be four hand-picked bundles (multiagent,
# _v2, _kb, _v3) with no way to express "everything on" and no way to
# attribute an effect to one capability. The capabilities are in fact
# independent, so they are declared as axes here and the named methods
# become presets over them. That makes a factorial ablation expressible,
# which is what an attribution question needs.
#
# Order matters: it fixes the bit order of the `abl<bits>` labels.
MA_AXES = [
    # code,     description
    ("mining", "run the source_mining stage over the called subroutines and "
               "append the conventions it recovers to the generation prompt. "
               "Reads source the entry routine does not contain, e.g. that "
               "GMMATD stores every matrix row-major, which is documented in "
               "gmmatd.f and nowhere in ktrplt.f"),
    ("kb",     "append knowledge-base insights harvested from runs that "
               "passed on OTHER benchmarks. Cross-benchmark transfer, so "
               "beware of leakage if harvest and evaluation share data"),
    ("lean",   "skip the orchestrator's code_generation and equivalence "
               "stages. The former asks CodeGenAgent for Python while the "
               "benchmark compiles C++, and its output is discarded; the "
               "latter then checks that discarded Python against zero test "
               "cases. Measured at 55.8% of a multiagent repeat's tokens"),
    ("wide",   "widen the pipeline-to-generator channel: design budget "
               "2500 -> 6000 chars and additionally forward physics_spec and "
               "the recovered call graph, neither of which reached the "
               "generator at all before"),
    ("srcctx", "derive the pipeline's root context from the source instead of "
               "config['description']. That field is a hand-written label and "
               "can be wrong: B4's asserted a two-triangle decomposition while "
               "KQDMEM uses four overlapping corner triangles at half "
               "thickness, seeding the whole spec chain with the answer the "
               "oracle rejects"),
    ("selfmine", "extract the implementation conventions of the routine BEING "
                 "TRANSLATED - DATA tables, literal scale factors, loop guards "
                 "and writes into argument slots - and state them explicitly. "
                 "mining only looks outward at callees, so for B4 it returned "
                 "seven conventions about GMMATD/SMA1B/KTRMEM and nothing "
                 "about the four-triangle DATA M map, ECPT(8)/2.0 or JNOT, all "
                 "of which are in kqdmem.f itself and all of which decide the "
                 "answer. Extraction is deterministic, so it cannot invent a "
                 "convention"),
]
MA_CODES = [c for c, _ in MA_AXES]

# A method label names the modules it switches on, joined by '+' in MA_AXES
# order, e.g. multiagent_mining+lean+wide. Two shorthands:
#   multiagent        every axis off. Kept verbatim because the existing clean
#                     baseline (22/35) was recorded under this directory name.
#   multiagent_full   every axis on.
# Opaque version numbers are gone; the label now says what ran.
MA_BASELINE = "multiagent"
MA_FULL = "multiagent_full"

# Retired labels, still accepted so archived run directories and older
# commands keep resolving. v2/v3 said nothing about what was enabled.
MA_ALIASES = {
    "multiagent_v2": "multiagent_mining",
    "multiagent_kb": "multiagent_mining+kb",
    "multiagent_v3": "multiagent_mining+lean+wide+srcctx",
    # what multiagent_full meant before the selfmine axis existed
    "multiagent_full4": "multiagent_mining+kb+lean+wide+srcctx",
}


def _axes_from_label(method: str) -> set:
    if method == MA_BASELINE:
        return set()
    if method == MA_FULL:
        return set(MA_CODES)
    method = MA_ALIASES.get(method, method)
    if method == MA_FULL:
        return set(MA_CODES)
    if not method.startswith("multiagent_"):
        raise ValueError(f"not a multi-agent method: {method!r}")
    body = method[len("multiagent_"):]
    # legacy bitmask form, e.g. multiagent_abl10110
    m = re.fullmatch(r"abl([01]{%d})" % len(MA_CODES), body)
    if m:
        return {c for c, b in zip(MA_CODES, m.group(1)) if b == "1"}
    parts = [p for p in body.split("+") if p]
    unknown = [p for p in parts if p not in MA_CODES]
    if unknown:
        raise ValueError(
            f"unknown module(s) {unknown} in method {method!r}. "
            f"Valid modules: {MA_CODES}. Compose them with '+', "
            f"e.g. multiagent_mining+wide, or use {MA_BASELINE} / {MA_FULL}")
    return set(parts)


def resolve_ma_flags(method: str) -> dict:
    """Map a method label to the capability dict."""
    return {c: (c in _axes_from_label(method)) for c in MA_CODES}


def ma_label(flags: dict) -> str:
    """Canonical label for a flag set: the module names it switches on."""
    on = [c for c in MA_CODES if flags.get(c)]
    if not on:
        return MA_BASELINE
    if len(on) == len(MA_CODES):
        return MA_FULL
    return "multiagent_" + "+".join(on)


def format_ma_flags(flags: dict) -> str:
    return " ".join(f"{c}={'on' if flags.get(c) else 'off'}" for c in MA_CODES)


# Named entry points offered on the CLI and in --ma-axes. Any other
# combination is addressable by composing module names.
MA_PRESETS = {
    name: _axes_from_label(name) for name in (
        MA_BASELINE,
        "multiagent_mining",
        "multiagent_mining+kb",
        "multiagent_selfmine",
        "multiagent_mining+lean+wide+srcctx",
        # canonical order is MA_AXES order, so lean precedes srcctx precedes
        # selfmine; ma_label() emits this form and the round-trip must hold
        "multiagent_lean+srcctx+selfmine",
        MA_FULL,
    )
}

# Global knowledge base instance (initialized in main if needed)
_knowledge_base: Optional[KnowledgeBase] = None


def load_all_benchmarks() -> dict:
    with open(CONFIG_FILE, "r", encoding="utf-8") as f:
        return json.load(f)

def load_benchmark_config(benchmark_id: str) -> dict:
    for b in load_all_benchmarks()["benchmarks"]:
        if b["id"] == benchmark_id:
            return b
    raise ValueError(f"Unknown benchmark: {benchmark_id}")

def read_fortran_source(config: dict) -> str:
    code = ""
    for src in config["source_files"]:
        p = PROJECT_ROOT / src
        if p.exists():
            code += p.read_text(encoding="utf-8", errors="replace") + "\n\n"
    return code


# A self-contained Eigen program with main() that prints a matrix is never this
# small. The repair loop used to accept anything over 50 bytes, so a 66-byte
# fragment - prose, a stray fence, a snippet - became the new candidate and
# every later iteration repaired the fragment instead of the working code.
# Measured across the pre-guard cohort (run directories stamped before
# 2026-08-22): 129 of 1843 iterations (7.0%) took a candidate under 500 bytes,
# and 30 repeats ended on one. See audit_degenerate_artifacts.py.
MIN_PLAUSIBLE_CODE = 500


def _is_plausible_program(code: str) -> Optional[str]:
    """Return a reason string if this cannot be a complete program."""
    if not code or not code.strip():
        return "empty response"
    if code.startswith("ERROR:"):
        return "generator returned an error marker"
    n = len(code)
    if n < MIN_PLAUSIBLE_CODE:
        return f"only {n} bytes, below the {MIN_PLAUSIBLE_CODE}-byte floor"
    if "int main" not in code:
        return "no int main(), so it is a fragment rather than a program"
    return None


def _accept_candidate(new_code: str, current: str, best: str,
                      logger, tag: str) -> str:
    """
    Decide what the next iteration should work on.

    Rejecting a degenerate candidate is not enough: the loop must also go back
    to the best known-good code, otherwise a collapse early in the sequence
    poisons every later repair. Sequences like
        17609 -> 17609 -> 17675 -> 17679 -> 16668 -> 147
    were common before this guard.
    """
    reason = _is_plausible_program(new_code)
    if reason is None:
        return new_code
    fallback = best if _is_plausible_program(best) is None else current
    logger.log_info(f"    rejected {tag} candidate: {reason}; "
                    f"keeping the previous best program")
    return fallback


def _leading_comment_block(fortran_code: str, max_chars: int = 4000) -> str:
    """
    The comment header the original authors wrote, up to the first executable
    statement. For the NASTRAN element kernels this is where the ECPT layout is
    documented, so it is the most reliable non-executable description of the
    routine's interface that exists.
    """
    out = []
    for line in fortran_code.splitlines():
        s = line.rstrip()
        if not s:
            continue
        if s[0] in "Cc*" or s.lstrip().upper().startswith("SUBROUTINE"):
            out.append(s)
            if sum(len(x) + 1 for x in out) > max_chars:
                break
        elif out:
            break
    return "\n".join(out)[:max_chars]


def _orchestrator_context(config: dict, fortran_code: str) -> str:
    """
    Root context for the multi-agent pipeline (multiagent_v3).

    The earlier variants passed f"{config['description']}. Target: C++17/Eigen."
    as the pipeline's only context. That field is a hand-written label, and it
    can be wrong: B4's said "2-triangle decomposition" while KQDMEM uses four
    overlapping corner triangles at half thickness, so the comprehension and
    physics_spec stages were seeded with the very answer the oracle rejects,
    and 20 of 47 archived products reproduced it.

    So the algorithm is described only by the source. What we assert here is
    limited to the verification interface, which is ground truth from the
    harness, plus the original authors' own comment header. The description is
    still passed through but explicitly marked as unverified so no downstream
    agent treats it as authoritative.
    """
    header = _leading_comment_block(fortran_code)
    parts = [
        "Task: translate the NASTRAN-95 Fortran below to C++17/Eigen, "
        "preserving numerical behaviour.",
        "",
        "Verification interface. These are FIXED by the harness - by the "
        "reference matrix and the Fortran test driver that produced it - and "
        "are not open to reinterpretation:",
        f"  entry routine   : {Path(config['source_files'][0]).stem.upper()}",
        f"  source files    : {', '.join(Path(s).name for s in config['source_files'])}",
        f"  verification    : {config.get('verification_type')}",
    ]
    if config.get("num_nodes") and config.get("dof_per_node"):
        parts.append(f"  nodes           : {config['num_nodes']}")
        parts.append(f"  DOF per node    : {config['dof_per_node']}")
    if config.get("matrix_size"):
        parts.append(f"  compared matrix : {config['matrix_size']}"
                     f"x{config['matrix_size']}")
    if config.get("compared_dof_note"):
        parts.append(f"  scope           : {config['compared_dof_note']}")
    if config.get("heat_branch_active") is False:
        parts.append(
            "  HEAT branch     : NOT exercised. The driver sets HEAT = .FALSE., "
            "so every heat-conduction path in the source is dead code here. Do "
            "not let it influence the specification or the generated code.")
    parts += [
        "",
        "The following is an unverified human label. Do NOT treat it as a "
        "specification of the ALGORITHM - one of these labels was found to "
        "name the wrong decomposition entirely - but the interface numbers "
        f"above override it either way: \"{config.get('description', '')}\"",
        "",
        "The algorithm is defined solely by the Fortran source. Where the "
        "source and the label disagree about the algorithm, the source "
        "governs. Pay particular attention to literal constants, DATA "
        "statements, and branch conditions, because those encode "
        "implementation conventions that textbook formulations do not share.",
    ]
    if header:
        parts += ["", "Original comment header from the source:", header]
    return "\n".join(parts)


def _extract_entry_conventions(fortran_code: str) -> list[dict]:
    """
    Pull the implementation conventions out of the routine BEING TRANSLATED.

    The source_mining stage looks outward at callees, so for B4 it recovered
    seven conventions about GMMATD, SMA1B and KTRMEM and not one word about the
    three facts that actually decide correctness, all of which live in
    kqdmem.f itself:

        DATA M / 1, 2, 4, 2, 3, 1, 3, 4, 2, 4, 1, 3 /   four overlapping
                                                        corner triangles
        ECPT(7) = ECPT(8)/2.0                           thickness halved
        IF (J .EQ. JNOT) GO TO 150                      one triangle skipped
                                                        per pivot

    Measured consequence: the halving reached generated code in 0 of 23
    attempts, and the multiagent_full repeat that did implement the map and the
    skip was still off by a factor of two.

    Extraction is deterministic - regex over the source, no model call - so it
    is reproducible and cannot invent a convention that is not there.
    """
    out = []
    lines = fortran_code.splitlines()

    def add(kind, lineno, text, why, role):
        out.append({"kind": kind, "line": lineno, "role": role,
                    "code": text.strip(), "why": why})

    # Names that mean "this routine is being called once per pivot grid point".
    # Guards mentioning them describe the legacy CALL STRUCTURE, not the
    # element, and must not be reproduced by code that builds the whole matrix
    # in one pass. Measured on B4: the two repeats that implemented the JNOT
    # skip literally left two reference cells at zero and scored 0.50, while
    # the three that ignored it scored 2.7e-03, 2.7e-03 and 2.45e-07.
    PIVOTISH = re.compile(r"\b(NPVT|NPIVOT|JNOT|IPVT|NPT\d?)\b", re.I)
    # Degeneracy and error checks. Real code, but they say nothing about the
    # numerical convention, and they crowded out the facts that matter.
    DEGENERATE = re.compile(r"\.EQ\.\s*0\.0|\.LE\.\s*0\.0|\.EQ\.\s*0\b"
                            r"|VECL|DETERM|\bZERO\b", re.I)

    for i, raw in enumerate(lines, 1):
        # skip comments and continuation-only lines
        if not raw.strip() or raw[0] in "Cc*!":
            continue
        code = raw[6:] if len(raw) > 6 else raw
        s = code.strip()
        if not s:
            continue

        # DATA statements: hardcoded tables, often node maps or coefficients
        if re.match(r"DATA\s+", s, re.I):
            add("table", i, s,
                "a hardcoded table. If it maps node or DOF indices it defines "
                "the element decomposition, which will not match any textbook "
                "default",
                "implement")
            continue

        # arithmetic that scales a quantity by a literal other than 1
        m = re.search(r"=\s*([A-Z0-9_()]+)\s*([*/])\s*([0-9]*\.?[0-9]+)", s, re.I)
        if m and m.group(3) not in ("1", "1.", "1.0"):
            add("scale", i, s,
                f"scales by {m.group(2)}{m.group(3)}. A factor like this is an "
                "implementation convention, not a physical constant, and is "
                "the single most frequently dropped detail",
                "implement")
            continue

        # guards that skip work
        if re.search(r"IF\s*\(.*\)\s*(GO\s*TO|CYCLE|RETURN)", s, re.I):
            if not re.search(r"\.EQ\.|\.NE\.|\.LT\.|\.GT\.", s, re.I):
                continue
            if DEGENERATE.search(s):
                # a zero-length or singularity check; carries no convention
                continue
            if PIVOTISH.search(s):
                add("pivot_guard", i, s,
                    "gates work on WHICH PIVOT this call is for. The legacy "
                    "routine is invoked once per grid point and each call "
                    "handles one row block, so this condition belongs to the "
                    "call structure, not to the element",
                    "context")
            else:
                add("guard", i, s,
                    "a conditional that skips or exits. Check whether it "
                    "changes which sub-elements contribute to the result",
                    "implement")
            continue

        # Writes that REPACK one argument slot from another, e.g.
        #     ECPT(5) = ECPT(6)
        # These shift a routine's argument block into the layout its callee
        # expects, so the callee reads a different quantity than the name
        # suggests. Plain arithmetic into a slot, IVEC(1) = ECPT(15) - ECPT(11)
        # and the like, is ordinary computation and is deliberately excluded:
        # it made up 38 of 49 hits on B4 and pushed the three facts that matter
        # past the injection budget.
        m = re.match(r"([A-Z]+)\(\s*(\d+)\s*\)\s*=\s*"
                     r"([A-Z]+)\(\s*(\d+)\s*\)\s*$", s, re.I)
        if m and m.group(1).upper() == m.group(3).upper():
            add("repack", i, s,
                f"moves {m.group(1)}({m.group(4)}) into "
                f"{m.group(1)}({m.group(2)}). The argument block is being "
                "rearranged for a called routine, so that callee reads this "
                "slot as a different quantity. Use this to work out WHAT VALUE "
                "the callee receives; do not reproduce the shuffling itself",
                "context")

    return out


def _format_entry_conventions(items: list[dict], budget: int = 3500) -> str:
    if not items:
        return ""
    order = ["scale", "table", "guard", "repack", "pivot_guard"]

    def key(d):
        return (order.index(d["kind"]) if d["kind"] in order else 9, d["line"])

    seen = set()
    impl, ctx = [], []
    for d in sorted(items, key=key):
        k = (d["kind"], d["code"])
        if k in seen:
            continue
        seen.add(k)
        line = f"  [line {d['line']}] {d['code']}\n      -> {d['why']}"
        (impl if d.get("role", "implement") == "implement" else ctx).append(line)

    parts = ["\nIMPLEMENTATION CONVENTIONS FOUND IN THE ROUTINE YOU ARE "
             "TRANSLATING.\nExtracted mechanically from the source above and "
             "split by what you should do with them."]

    if impl:
        parts.append(
            "\nA. REPRODUCE THESE. Each changes the numerical result and no "
            "textbook formulation has them:\n" + "\n".join(impl))
    if ctx:
        parts.append(
            "\nB. DO NOT REPRODUCE THESE. They exist because the legacy "
            "routine is driven one pivot grid point at a time and receives its "
            "arguments through a shared array. Your program computes the whole "
            "matrix in a single pass, so it has no pivot loop and no argument "
            "block. Read them only to work out what the original computes, "
            "then implement that. Copying them literally drops real "
            "contributions and leaves reference entries at zero:\n"
            + "\n".join(ctx))

    return "\n".join(parts)[:budget]


def _spec_contradicts_harness(spec: dict, config: dict) -> Optional[str]:
    """
    Is the recovered specification inconsistent with the fixed interface?

    Added after the wide-channel cell scored 0/10 on B4 with rel=inf and rel>1.
    The pipeline had produced

        dof_per_node : "3 for structural (translational only), 1 for heat"
        total_dof    : "12 for structural, 4 for heat"

    for a benchmark compared as 8x8, and forwarding that text told the
    generator to build the wrong-sized matrix. A contradictory specification is
    worse than none, so the caller drops it instead of injecting it.
    """
    if not isinstance(spec, dict):
        return None

    def as_int(v):
        if isinstance(v, bool):
            return None
        if isinstance(v, int):
            return v
        if isinstance(v, float) and v.is_integer():
            return int(v)
        if isinstance(v, str):
            nums = re.findall(r"\d+", v)
            # a field that names more than one number is hedging
            if len(set(nums)) > 1:
                return "AMBIGUOUS"
            if nums:
                return int(nums[0])
        return None

    for field, expected in (("dof_per_node", config.get("dof_per_node")),
                            ("num_nodes", config.get("num_nodes")),
                            ("total_dof", config.get("compared_dof"))):
        if expected is None:
            continue
        got = as_int(spec.get(field))
        if got == "AMBIGUOUS":
            return (f"{field} is ambiguous ({spec.get(field)!r}); the harness "
                    f"fixes it at {expected}")
        if got is not None and got != expected:
            return (f"{field} is {got}, the harness fixes it at {expected}")
    return None

# ============================================================
# Reference Loading
# ============================================================

def load_reference(config: dict) -> Optional[np.ndarray]:
    """Load reference - analytical for B1-B4, file-based for B5+."""
    if config.get("reference_method") == "compute_analytical":
        key_map = {"B1": "krod", "B2": "kbar", "B3": "ktrmem", "B4": "ktriqd"}
        k = key_map.get(config["id"])
        return compute_analytical_reference(k) if k else None
    ref_path = config.get("reference_output_path")
    if not ref_path:
        return None
    full = PROJECT_ROOT / ref_path
    if not full.exists():
        print(f"  WARNING: Reference not found: {full}")
        return None
    text = full.read_text(encoding="utf-8", errors="replace")
    return _parse_reference(text, config)


# Semantic slots for scalar_compare, each with accepted aliases. Values are
# placed BY SLOT, never in key-encounter order: the previous extractor appended
# whichever of a flat key list it found, so a program emitting
# {"stress","force","torsional_stress"} yielded [torsional, stress, force] --
# the correct answer rotated by one, scoring rel_diff = 1.0. See
# diagnose_b13.py and 论文一_实验清单.md.
_SCALAR_SLOTS = [
    ("axial_stress", "stress", "sigma", "sigma_axial"),
    ("axial_force", "force", "N", "axial_load"),
    ("torsional_stress", "torsion", "tau", "shear_stress"),
]


def _extract_scalars(data: dict) -> Optional[np.ndarray]:
    """
    Pull scalar_compare values into fixed semantic positions.

    Returns None if no slot can be resolved. Slots that are absent become NaN
    so that a partial answer cannot masquerade as a shorter-but-aligned array;
    the verifier treats NaN as a mismatch rather than silently truncating.
    """
    vals, found = [], 0
    for aliases in _SCALAR_SLOTS:
        v = None
        for k in aliases:
            if k in data:
                try:
                    v = float(data[k])
                except (TypeError, ValueError):
                    v = None
                break
        if v is None:
            vals.append(float("nan"))
        else:
            vals.append(v)
            found += 1
    if found == 0:
        return None
    return np.array(vals, dtype=float)


def _parse_reference(text: str, config: dict) -> Optional[np.ndarray]:
    """Parse reference file into numpy array based on verification type."""
    vtype = config["verification_type"]
    # Try JSON first (preferred format)
    try:
        data = json.loads(text)
        if vtype == "vector_compare":
            for k in ["solution", "x", "vector", "values", "result"]:
                if k in data:
                    return np.array(data[k], dtype=float)
        elif vtype == "scalar_compare":
            got = _extract_scalars(data)
            if got is not None:
                return got
        elif vtype == "eigenvalue_compare":
            for k in ["eigenvalues", "evals", "values"]:
                if k in data:
                    return np.sort(np.array(data[k], dtype=float))
        elif vtype == "timeseries_compare":
            for k in ["timeseries", "solution", "values"]:
                if k in data:
                    return np.array(data[k], dtype=float)
        elif vtype == "matrix_symmetry":
            for k in ["stiffness_matrix", "matrix", "K", "inverse"]:
                if k in data:
                    return np.array(data[k], dtype=float)
        # Generic fallback for JSON
        for k in ["solution", "values", "result", "matrix", "eigenvalues"]:
            if k in data:
                return np.array(data[k], dtype=float)
    except (json.JSONDecodeError, ValueError, TypeError):
        pass

    # Row-based numeric parsing (for text reference files)
    rows = []
    for line in text.strip().split("\n"):
        line = line.strip()
        if not line or line[0] in "#!":
            continue
        nums = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][-+]?\d+)?", line)
        if nums:
            rows.append([float(n.replace("D", "E").replace("d", "e")) for n in nums])
    if not rows:
        return None
    if vtype == "matrix_symmetry" and config.get("matrix_size"):
        sz = config["matrix_size"]
        mat_rows = [r[:sz] for r in rows if len(r) >= sz]
        if len(mat_rows) >= sz:
            return np.array(mat_rows[:sz], dtype=float)
    if vtype == "eigenvalue_compare":
        return np.sort(np.array([v for row in rows for v in row], dtype=float))
    if vtype == "timeseries_compare":
        # Timeseries may have rows of different lengths; pad or truncate
        max_cols = max(len(r) for r in rows) if rows else 0
        padded = [r + [0.0]*(max_cols - len(r)) for r in rows]
        return np.array(padded, dtype=float)
    return np.array([v for row in rows for v in row], dtype=float)

# ============================================================
# Output Parsing & Verification
# ============================================================

def parse_output(stdout: str, config: dict) -> Optional[np.ndarray]:
    """Parse program stdout based on verification type."""
    vtype = config["verification_type"]
    if vtype == "matrix_symmetry":
        return parse_matrix_from_output(stdout, config["matrix_size"])

    # Try JSON extraction for all other types
    try:
        brace = stdout.find("{")
        if brace >= 0:
            depth = 0
            for i in range(brace, len(stdout)):
                if stdout[i] == "{": depth += 1
                elif stdout[i] == "}":
                    depth -= 1
                    if depth == 0:
                        data = json.loads(stdout[brace:i+1])
                        if vtype == "scalar_compare":
                            got = _extract_scalars(data)
                            if got is not None:
                                return got
                        if vtype == "vector_compare":
                            for k in ["solution", "x", "vector", "values", "result"]:
                                if k in data:
                                    return np.array(data[k], dtype=float)
                        if vtype == "eigenvalue_compare":
                            for k in ["eigenvalues", "evals", "values"]:
                                if k in data:
                                    return np.sort(np.array(data[k], dtype=float))
                        if vtype == "timeseries_compare":
                            for k in ["timeseries", "solution", "values"]:
                                if k in data:
                                    return np.array(data[k], dtype=float)
                        # Generic
                        for k in ["solution", "x", "vector", "eigenvalues",
                                   "evals", "stress", "force", "values", "result"]:
                            if k in data:
                                v = data[k]
                                if isinstance(v, (int, float)):
                                    return np.array([v], dtype=float)
                                return np.array(v, dtype=float)
                        break
    except (json.JSONDecodeError, ValueError):
        pass

    # Fallback: numeric row parsing
    rows = []
    for line in stdout.strip().split("\n"):
        line = line.strip()
        if not line or line[0] in "#!{[":
            continue
        nums = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?", line)
        if nums:
            rows.append([float(n) for n in nums])

    if not rows:
        return None

    if vtype == "timeseries_compare":
        # Rows may have different lengths; pad to max
        max_cols = max(len(r) for r in rows)
        padded = [r + [0.0]*(max_cols - len(r)) for r in rows]
        return np.array(padded, dtype=float)
    if vtype == "eigenvalue_compare":
        vals = sorted(v for row in rows for v in row)
        return np.array(vals, dtype=float)
    # vector/scalar: flatten
    return np.array([v for row in rows for v in row], dtype=float)

def verify_output(computed: Optional[np.ndarray], reference: Optional[np.ndarray],
                  config: dict) -> dict:
    """Type-specific verification. Returns {passed, max_rel_diff, ...}."""
    if computed is None:
        return {"passed": False, "error": "Could not parse output"}
    if reference is None:
        return {"passed": False, "error": "No reference available"}

    tol_map = load_all_benchmarks().get("tolerance", {})
    tol = tol_map.get(config["verification_type"], 1e-4)
    vtype = config["verification_type"]

    if vtype == "matrix_symmetry":
        # Ensure both are 2D
        if computed.ndim == 1:
            sz = config.get("matrix_size")
            if sz and len(computed) >= sz*sz:
                computed = computed[:sz*sz].reshape(sz, sz)
            else:
                return {"passed": False, "error": f"Cannot reshape computed (len={len(computed)}) to {sz}x{sz}"}
        if reference.ndim == 1:
            sz = config.get("matrix_size")
            if sz and len(reference) >= sz*sz:
                reference = reference[:sz*sz].reshape(sz, sz)
            else:
                return {"passed": False, "error": f"Cannot reshape reference (len={len(reference)}) to {sz}x{sz}"}

        # Handle shape mismatches between computed and reference
        if computed.shape != reference.shape:
            cs, rs = computed.shape[0], reference.shape[0]
            if cs == 12 and rs == 6:
                idx = [0, 1, 2, 6, 7, 8]
                computed = computed[np.ix_(idx, idx)]
            elif cs == 6 and rs == 12:
                idx = [0, 1, 2, 6, 7, 8]
                reference = reference[np.ix_(idx, idx)]
            elif cs < rs:
                reference = reference[:cs, :cs]
            elif cs > rs:
                computed = computed[:rs, :rs]
            if computed.shape != reference.shape:
                return {"passed": False,
                        "error": f"Shape mismatch: {computed.shape} vs {reference.shape}"}
        max_ref = max(np.max(np.abs(reference)), 1e-15)
        abs_d = float(np.max(np.abs(computed - reference)))
        rel_d = abs_d / max_ref
        return {"passed": rel_d < tol or abs_d < 1e-10,
                "max_abs_diff": abs_d, "max_rel_diff": rel_d,
                "is_symmetric": bool(np.allclose(computed, computed.T, atol=1e-10))}

    elif vtype == "eigenvalue_compare":
        c = np.sort(computed.flatten())
        r = np.sort(reference.flatten())
        # Eigenvalue outputs have an exact contract. Silently truncating to
        # min(len(computed), len(reference)) let partial or extra outputs pass
        # on a matching subset, which is invalid for an equivalence claim.
        if c.size != r.size:
            return {
                "passed": False,
                "contract_violation": True,
                "reason": (f"expected exactly {r.size} eigenvalues, got "
                           f"{c.size}"),
                "computed_shape": [int(c.size)],
                "reference_shape": [int(r.size)],
                "max_abs_diff": None,
                "max_rel_diff": None,
            }
        if not np.isfinite(c).all() or not np.isfinite(r).all():
            return {
                "passed": False,
                "contract_violation": True,
                "reason": "eigenvalue output contains NaN or Inf",
                "computed_shape": [int(c.size)],
                "reference_shape": [int(r.size)],
                "max_abs_diff": None,
                "max_rel_diff": None,
            }
    elif vtype == "timeseries_compare":
        # Do NOT silently truncate to the smaller shape. Truncation compared
        # misaligned rows and columns (step 0..10 of (time,u1,u2) against
        # snapshots 0..10 of (u1,u2,u3)) and reported the result as a
        # numerical error, hiding the fact that the output contract was
        # never satisfied. Report the mismatch instead.
        cshape = computed.shape if computed.ndim > 1 else (computed.shape[0], 1)
        rshape = reference.shape if reference.ndim > 1 else (reference.shape[0], 1)
        if cshape != rshape:
            return {
                "passed": False,
                "contract_violation": True,
                "reason": (f"output shape {cshape} does not match the "
                           f"reference shape {rshape}; the comparison is "
                           f"ill-defined and this run carries no numerical "
                           f"information"),
                "computed_shape": list(cshape),
                "reference_shape": list(rshape),
                "max_abs_diff": None,
                "max_rel_diff": None,
            }
        c = computed if computed.ndim > 1 else computed.reshape(-1, 1)
        r = reference if reference.ndim > 1 else reference.reshape(-1, 1)
    else:  # vector_compare, scalar_compare
        cf = computed.flatten()
        rf = reference.flatten()
        if vtype == "scalar_compare":
            # Slots are semantic, so a missing one is a contract violation
            # rather than something to truncate around.
            if cf.size != rf.size or np.isnan(cf).any():
                missing = [i for i in range(min(cf.size, rf.size))
                           if i < cf.size and np.isnan(cf[i])]
                return {
                    "passed": False,
                    "contract_violation": True,
                    "reason": (f"expected {rf.size} named scalars, resolved "
                               f"{int(np.count_nonzero(~np.isnan(cf)))}"
                               + (f"; unresolved slots {missing}" if missing else "")),
                    "computed_shape": [int(cf.size)],
                    "reference_shape": [int(rf.size)],
                    "max_abs_diff": None,
                    "max_rel_diff": None,
                }
            c, r = cf, rf
        else:
            n = min(cf.size, rf.size)
            c, r = cf[:n], rf[:n]

    max_ref = max(np.max(np.abs(r)), 1e-15)
    abs_d = float(np.max(np.abs(c - r)))
    rel_d = abs_d / max_ref
    return {"passed": rel_d < tol or abs_d < 1e-10,
            "max_abs_diff": abs_d, "max_rel_diff": rel_d}

# ============================================================
# Prompt Generation
# ============================================================

def _timeseries_spec(config: dict) -> str:
    """
    Output contract for time-history benchmarks.

    The contract must be stated exactly. Previously the prompt asked for
    "one row per timestep: time u1 u2 u3", so translations emitted 101 rows
    of 4 columns while the reference held 11 sampled rows of 3 displacements.
    The verifier truncated both to the smaller shape and compared misaligned
    data, turning a specification gap into an apparent numerical failure.
    """
    c = config.get("output_contract") or {}
    n = c.get("n_samples")
    dof = c.get("n_dof", 3)
    stride = c.get("sample_stride")
    total = c.get("total_steps")
    name = config["name"]

    if not (n and stride and total):
        # No contract declared: fall back to an explicit JSON request.
        return (f'Print the displacement history as JSON:\n'
                f'{{"test":"{name}","timeseries":[[u1,u2,u3], ...]}}\n'
                f'One inner list per reported step. Do not include a time '
                f'or step column.')

    steps = ", ".join(str(i * stride) for i in range(min(n, 4))) + ", ..., " + str(total)
    return (
        f'Print the displacement history as JSON:\n'
        f'{{"test":"{name}","timeseries":[[u1,u2,u3], ...]}}\n'
        f'Report EXACTLY {n} rows: the state at steps {steps} '
        f'(i.e. every {stride}th step of {total}, including step 0 and step {total}).\n'
        f'Each row holds exactly {dof} displacement values.\n'
        f'Do NOT include a time column, a step column, or the intermediate '
        f'steps. Use scientific notation. The JSON must be valid and parseable.'
    )


def _output_spec(config: dict) -> str:
    """Output format specification per verification type."""
    vtype = config["verification_type"]
    name = config["name"]
    sz = config.get("matrix_size")
    specs = {
        "matrix_symmetry": f'Print the {sz}x{sz} stiffness matrix as JSON to stdout:\n{{"stiffness_matrix":[[row1_values],[row2_values],...]}}\nUse scientific notation (e.g., 1.5e+09). The JSON must be valid and parseable.',
        "vector_compare": f'Print solution vector as JSON: {{"test":"{name}","solution":[x0,x1,...]}}',
        "eigenvalue_compare": f'Print sorted eigenvalues as JSON: {{"test":"{name}","eigenvalues":[l0,l1,...]}}',
        "timeseries_compare": _timeseries_spec(config),
        # Key names must match the reference and _TEST_PARAMS exactly. They
        # previously did not: this spec asked for "stress"/"force" while the
        # test parameters and the reference used "axial_stress"/"axial_force".
        # A model that satisfied both emitted a mix, and the extractor -- which
        # appended keys in ITS OWN fixed order -- returned a rotated array, so
        # a correct answer scored rel_diff = 1.0.
        "scalar_compare": (
            f'Print scalars as JSON:\n'
            f'{{"test":"{name}","axial_stress":val,"axial_force":val,'
            f'"torsional_stress":val}}\n'
            f'Use exactly these three keys. Use scientific notation. '
            f'The JSON must be valid and parseable.'),
    }
    return specs.get(vtype, "Print results to stdout.")


# Test parameters that match the reference solutions exactly
_TEST_PARAMS = {
    "B1": """Test case (USE THESE EXACT VALUES):
  Node A: (0, 0, 0), Node B: (2, 0, 0), L=2.0 m
  E = 200e9 Pa, A = 0.01 m^2, G = 76.923e9 Pa, J = 5e-6 m^4
  Expected K(1,1) = E*A/L = 1e9 N/m
  Output: 12x12 stiffness matrix (6 DOF per node: ux,uy,uz,rx,ry,rz)""",
    "B2": """Test case (USE THESE EXACT VALUES):
  Node A: (0, 0, 0), Node B: (2, 0, 0), L=2.0 m
  E = 200e9 Pa, G = 76.923e9 Pa, A = 0.01 m^2
  Iy = 8.333e-6 m^4, Iz = 8.333e-6 m^4, J = 1.667e-5 m^4
  Output: 12x12 stiffness matrix (6 DOF per node)""",
    "B3": """Test case (USE THESE EXACT VALUES):
  Nodes: (0,0) m, (2,0) m, (1,1.5) m
  E = 2.1e11 Pa, nu = 0.3, thickness t = 0.01 m
  Output: 6x6 stiffness matrix (2 DOF per node: ux, uy)""",
    "B4": """Test case (USE THESE EXACT VALUES):
  Nodes: (0,0), (2,0), (2,1.5), (0,1.5)
  E = 200e9 Pa, nu = 0.3, thickness t = 0.01 m
  Output: 8x8 stiffness matrix (2 DOF per node: ux, uy)""",
    "B5": """Test case (USE THESE EXACT VALUES):
  Triangle: A=(0,0,0), B=(1,0,0), C=(0,1,0)
  E = 200e9, nu = 0.3, plate thickness t = 0.01
  I = t^3/12 = 8.333e-8 (moment of inertia per unit width)
  Output: 9x9 bending stiffness matrix (3 DOF per node: w, theta_x, theta_y)""",
    "B6": """Test case (USE THESE EXACT VALUES):
  Quad: A=(0,0,0), B=(1,0,0), C=(1,1,0), D=(0,1,0)
  E = 200e9, nu = 0.3, plate thickness t = 0.01
  I = t^3/12 = 8.333e-8
  Output: 12x12 bending stiffness matrix (3 DOF per node: w, theta_x, theta_y)""",
    "B7": """Test case (USE THESE EXACT VALUES):
  Tetrahedron: N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
  E = 200e9, nu = 0.3
  Output: 12x12 stiffness matrix (3 DOF per node: ux, uy, uz)""",
    "B8": """Test case (USE THESE EXACT VALUES):
  Wedge (6 nodes): N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
  E = 200e9, nu = 0.3
  Output: 18x18 stiffness matrix (3 DOF per node)""",
    "B9": """Test case (USE THESE EXACT VALUES):
  90-degree elbow in XZ plane:
  Grid A: (1, 0, 0), Grid B: (0, 0, 1)
  Reference vector V = (0, 1, 0)
  Radius of curvature R = 1.0, Angle = 90 degrees
  E = 200e9, nu = 0.3, G = E/(2*(1+nu))
  Circular pipe: outer radius = 0.05, wall thickness = 0.005
  A = pi*(ro^2-ri^2), I1 = I2 = pi/4*(ro^4-ri^4), J = 2*I1
  Output: 12x12 stiffness matrix (6 DOF per node)""",
    "B10": """Test case (USE THESE EXACT VALUES):
  4x4 SPD matrix K = [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
  Solve Kx = b where b = K*[1,2,3,4]^T
  Also compute K_inverse and determinant.
  Output JSON: {"solution":[x1,x2,x3,x4],"determinant":value,"inverse":[[...],...]}""",
    "B11": """Test case (USE THESE EXACT VALUES):
  5x5 tridiagonal matrix: diagonal=[2,2,2,2,2], off-diagonal=[-1,-1,-1,-1]
  Compute ALL 5 eigenvalues using QR iteration.
  Analytical: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
  Output JSON: {"eigenvalues":[sorted list of 5 eigenvalues]}""",
    "B12": """Test case (USE THESE EXACT VALUES):
  3-DOF spring-mass system (both ends fixed):
  K = [2,-1,0; -1,2,-1; 0,-1,2], M = I (identity), C = 0
  Initial: u0 = [1, 0, 0], v0 = [0, 0, 0]
  Integrator: CENTRAL DIFFERENCE (explicit), dt=0.01, 100 steps.
    This is what NASTRAN TRD1C implements, and the reference was produced by
    running it. Do NOT use Newmark-beta: it is a different scheme and will
    not reproduce the reference.
      K_eff   = M/dt^2
      u_{n+1} = K_eff^-1 * [ (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} ]
    Derive the first step from the given initial conditions.""",
    "B13": """Test case (USE THESE EXACT VALUES):
  Rod along x-axis: L=2.0, A=0.01, E=2.1e11, nu=0.3
  J = 5e-6 (polar moment), C = J/max_radius = 0.005
  Displacements: node_a = (0,0,0,0,0,0), node_b = (0.001,0,0,0,0,0) (pure axial)
  Compute: axial_force, axial_stress, torsional_stress
  Output JSON: {"axial_stress":value,"axial_force":value,"torsional_stress":value}""",
}


def build_generation_prompt(config: dict, fortran_code: str) -> str:
    """Build initial C++ generation prompt."""
    test_params = _TEST_PARAMS.get(config["id"], "Use physically reasonable test parameters.")
    return f"""Translate the following NASTRAN-95 Fortran 77 subroutine(s) to modern C++17.

Element/Routine: {config['name']} - {config['description']}

Requirements:
1. Use Eigen library (#include <Eigen/Dense>) for matrices/vectors.
2. Single self-contained file with main() that runs the test case below.
3. Compilable with MSVC: cl.exe /std:c++17 /EHsc /O2
4. Double precision throughout. No external deps beyond Eigen + stdlib.

{test_params}

Output specification:
{_output_spec(config)}

Fortran source:
```fortran
{fortran_code}
```

Generate COMPLETE compilable C++ code.
"""


def _run_root_cause_analysis(llm, config: dict, code: str, 
                             computed: np.ndarray, reference: np.ndarray,
                             logger) -> str:
    """
    Automated root cause analysis when repair loop is stuck.
    
    Analyzes the structural pattern of numerical errors between computed
    and reference matrices, hypothesizes the code-level bug, and generates
    a targeted fix. This replicates the human expert analysis process.
    """
    benchmark_id = config["id"]
    sz = config.get("matrix_size", 12)
    
    # Analyze error structure
    analysis_lines = []
    max_ref = max(np.max(np.abs(reference)), 1e-15)
    
    # Reshape to matrix if possible
    if sz is None or sz == 0:
        # Non-matrix benchmark (eigenvalue, vector, etc.) — skip structural analysis
        analysis_lines.append("Non-matrix output — structural pattern analysis not applicable.")
        analysis_lines.append(f"Computed norm: {np.linalg.norm(computed):.4e}")
        analysis_lines.append(f"Reference norm: {np.linalg.norm(reference):.4e}")
        if np.linalg.norm(reference) > 0:
            ratio = np.linalg.norm(computed) / np.linalg.norm(reference)
            analysis_lines.append(f"Norm ratio (computed/reference): {ratio:.4f}")
        analysis_text = "\n".join(analysis_lines)
        # Build simple prompt
        prompt = f"""You are debugging a numerical computation.
The code compiles and runs but produces wrong results.

ELEMENT: {config['name']} ({config['description']})

ERROR ANALYSIS:
{analysis_text}

CURRENT CODE:
```cpp
{code}
```

Fix the numerical error. Return the COMPLETE fixed C++ file.
"""
        resp = llm.generate(
            "You are an expert debugger. Fix the numerical error.",
            prompt, f"root_cause_{benchmark_id}")
        return resp if resp else ""
    
    c = computed.reshape(sz, sz) if computed.ndim == 1 else computed
    r = reference.reshape(sz, sz) if reference.ndim == 1 else reference
    
    # Handle size mismatch: if computed is larger, truncate to reference size
    if c.shape != r.shape:
        min_sz = min(c.shape[0], r.shape[0])
        c = c[:min_sz, :min_sz]
        r = r[:min_sz, :min_sz]
        sz = min_sz
    
    # Update max_ref from reshaped matrix
    max_ref = max(np.max(np.abs(r)), 1e-15)
    
    # Block-level analysis (for element stiffness matrices)
    # 1. Block-level analysis (for element stiffness matrices)
    if sz % 3 == 0:
        n_nodes = sz // 3
        analysis_lines.append(f"Block analysis ({n_nodes} nodes, 3 DOF each):")
        for bi in range(n_nodes):
            for bj in range(n_nodes):
                block_c = c[3*bi:3*bi+3, 3*bj:3*bj+3]
                block_r = r[3*bi:3*bi+3, 3*bj:3*bj+3]
                block_err = np.max(np.abs(block_c - block_r)) / max_ref
                if block_err > 0.01:
                    # Check if this block is the transpose of what it should be
                    trans_err = np.max(np.abs(block_c - block_r.T)) / max_ref
                    is_transposed = trans_err < block_err * 0.1
                    
                    # Check if it matches some other block
                    swap_match = ""
                    for bk in range(n_nodes):
                        for bl in range(n_nodes):
                            if (bk, bl) != (bi, bj):
                                other_r = r[3*bk:3*bk+3, 3*bl:3*bl+3]
                                if np.max(np.abs(block_c - other_r)) / max_ref < 0.001:
                                    swap_match = f"matches reference block [{bk},{bl}]"
                                    break
                    
                    status = "TRANSPOSED" if is_transposed else f"err={block_err:.2e}"
                    analysis_lines.append(
                        f"  Block[{bi},{bj}]: {status} {swap_match}")
    
    # 2. Check if entire matrix is its own transpose
    if c.shape == r.shape:
        full_trans_err = np.max(np.abs(c - r.T)) / max_ref
        if full_trans_err < 0.01:
            analysis_lines.append(f"\nPATTERN: Entire computed matrix = reference TRANSPOSED (err={full_trans_err:.2e})")
            analysis_lines.append("This means the assembly loop has i,j indices swapped.")
    
    # 3. Check if off-diagonal blocks are individually transposed
    if sz % 3 == 0:
        n_nodes = sz // 3
        trans_blocks = 0
        total_off_diag = 0
        for bi in range(n_nodes):
            for bj in range(n_nodes):
                if bi != bj:
                    total_off_diag += 1
                    block_c = c[3*bi:3*bi+3, 3*bj:3*bj+3]
                    block_r = r[3*bi:3*bi+3, 3*bj:3*bj+3]
                    block_err = np.max(np.abs(block_c - block_r)) / max_ref
                    trans_err = np.max(np.abs(block_c - block_r.T)) / max_ref
                    if block_err > 0.01 and trans_err < block_err * 0.1:
                        trans_blocks += 1
        
        if trans_blocks > total_off_diag * 0.5:
            analysis_lines.append(
                f"\nPATTERN: {trans_blocks}/{total_off_diag} off-diagonal blocks are TRANSPOSED")
            analysis_lines.append(
                "ROOT CAUSE: The assembly convention is wrong.")
            analysis_lines.append(
                "In NASTRAN, K[i,j] = V * B[j]^T * D * B[i] (row=target, col=source).")
            analysis_lines.append(
                "FIX: Swap the indices: K.block<3,3>(3*i, 3*j) = V * B[j].transpose() * D * B[i]")
            analysis_lines.append(
                "Or equivalently: after computing K normally, set K = K.transpose()")
    
    # 4. Check for uniform scale factor
    nonzero_mask = np.abs(r) > max_ref * 0.001
    if np.sum(nonzero_mask) > 5:
        ratios = c[nonzero_mask] / r[nonzero_mask]
        mean_ratio = np.mean(ratios)
        std_ratio = np.std(ratios)
        cv = std_ratio / (abs(mean_ratio) + 1e-15)  # coefficient of variation
        
        if cv < 0.3 and abs(mean_ratio - 1.0) > 0.05:
            # Uniform scale factor detected
            analysis_lines.append(f"\nPATTERN: UNIFORM SCALE FACTOR detected!")
            analysis_lines.append(f"  Computed values are approximately {mean_ratio:.4f}x the expected values (CV={cv:.3f})")
            analysis_lines.append(f"  This means your result is consistently off by a factor of ~{mean_ratio:.3f}")
            analysis_lines.append(f"  The CORRECT result should be ~{1.0/mean_ratio:.3f}x your current output.")
            analysis_lines.append(f"\n  COMMON CAUSES for plate bending elements:")
            analysis_lines.append(f"  - Missing or double-counted thickness factor: I = t^3/12")
            analysis_lines.append(f"  - Wrong area factor in K^x: should multiply by 4*Area (Area of sub-triangle, NOT total plate)")
            analysis_lines.append(f"  - Sub-triangle area = total_plate_area / 3")
            analysis_lines.append(f"  - Missing T-matrix rotation when assembling sub-triangle contributions")
            analysis_lines.append(f"  - D matrix: D = I * G where I = t^3/12 and G = material stiffness (3x3 for plate)")
            
            if 0.45 < mean_ratio < 0.55 or 1.8 < mean_ratio < 2.2:
                analysis_lines.append(f"\n  SPECIFIC: ratio ~= {mean_ratio:.2f} suggests a factor-of-2 error.")
                analysis_lines.append(f"  CHECK: Are you using 4*Area or 2*Area? The correct factor is 4*Area.")
                analysis_lines.append(f"  CHECK: Sub-triangle area should be (x_b * y_c / 2), then K^x *= 4*area.")
            elif 0.3 < mean_ratio < 0.35 or 2.8 < mean_ratio < 3.5:
                analysis_lines.append(f"\n  SPECIFIC: ratio ~= {mean_ratio:.2f} suggests a factor-of-3 error.")
                analysis_lines.append(f"  CHECK: Are you dividing by 3 (for 3 sub-triangles) when you shouldn't be?")
                analysis_lines.append(f"  The sub-triangles each contribute their own stiffness; no division by 3 needed.")
        
        elif cv >= 0.3 and abs(mean_ratio - 1.0) > 0.3:
            # Non-uniform error but systematically wrong — likely assembly mapping issue
            analysis_lines.append(f"\nPATTERN: NON-UNIFORM SCALE ERROR (mean_ratio={mean_ratio:.2f}, CV={cv:.2f})")
            analysis_lines.append(f"  Some blocks are more wrong than others — suggests assembly mapping error.")
            analysis_lines.append(f"  COMMON CAUSES:")
            analysis_lines.append(f"  - Sub-triangle contributions assembled to wrong global DOF positions")
            analysis_lines.append(f"  - Missing coordinate rotation (T-matrix) when transforming sub-triangle K to global frame")
            analysis_lines.append(f"  - Wrong node-to-DOF mapping in static condensation")
            analysis_lines.append(f"  - The 3x3 T-matrix rotation T^T * K_block * T must be applied to EACH 3x3 block")
            analysis_lines.append(f"    before adding to the global matrix (this rotates from sub-triangle local to plate local coords)")
    
    # 5. Show specific element comparisons
    analysis_lines.append("\nTop wrong elements:")
    diffs = np.abs(c - r)
    flat_idx = np.argsort(diffs.flatten())[-5:][::-1]
    for idx in flat_idx:
        ri, ci = idx // sz, idx % sz
        if diffs[ri, ci] > max_ref * 0.01:
            analysis_lines.append(
                f"  K[{ri},{ci}]: computed={c[ri,ci]:.6e}, expected={r[ri,ci]:.6e}")
    
    analysis_text = "\n".join(analysis_lines)
    
    # Build the root cause analysis prompt
    prompt = f"""You are debugging a finite element stiffness matrix computation.
The code compiles and runs but produces numerically wrong results.

ELEMENT: {config['name']} ({config['description']})
MATRIX SIZE: {sz}x{sz}

ERROR ANALYSIS (automated):
{analysis_text}

CURRENT CODE:
```cpp
{code}
```

Based on the error pattern analysis above, identify the ROOT CAUSE of the numerical error.
Then provide the COMPLETE FIXED C++ code.

Focus on:
1. Assembly index convention (K[i,j] = B[i]^T*D*B[j] vs B[j]^T*D*B[i])
2. Matrix transpose issues in coordinate transforms
3. Row-major vs column-major storage when translating from Fortran

Return the COMPLETE fixed C++ file.
"""
    
    logger.log_info(f"    Root cause analysis triggered (stuck pattern detected)")
    logger.log_info(f"    Analysis: {analysis_text[:200]}")
    
    resp = llm.generate(
        "You are an expert debugger specializing in FEM stiffness matrix computations. "
        "Analyze the error pattern and fix the root cause.",
        prompt,
        f"root_cause_{benchmark_id}")
    
    return resp if resp else ""


def build_repair_prompt(config: dict, code: str, error_info: dict) -> str:
    """Build repair prompt based on failure stage."""
    stage = error_info.get("stage", "unknown")
    err = error_info.get("error", "")[:1500]

    if stage == "compile":
        return f"""C++ code for {config['name']} failed to compile.
Error:\n```\n{err}\n```\nCode:\n```cpp\n{code}\n```\nFix ALL errors. Return COMPLETE file."""

    elif stage == "runtime":
        return f"""C++ code for {config['name']} crashed at runtime.
Error:\n```\n{err}\n```\nCode:\n```cpp\n{code}\n```\nFix runtime error. Return COMPLETE file."""

    elif stage == "contract":
        # The numbers may well be right; only the emitted shape is wrong.
        # Telling the model to hunt for an arithmetic bug wastes the iteration.
        num = error_info.get("numerical", {})
        return f"""C++ code for {config['name']} compiles and runs, but the OUTPUT SHAPE
does not match what the verifier expects. The arithmetic may be correct --
do not change the numerical algorithm unless you find an actual defect.

Problem: {err}
Emitted shape:  {num.get('computed_shape')}
Expected shape: {num.get('reference_shape')}

Required output format:
{_output_spec(config)}

Code:\n```cpp\n{code}\n```

Adjust ONLY the output stage so the emitted shape matches exactly.
Return the COMPLETE fixed C++ file."""

    else:  # numerical
        num = error_info.get("numerical", {})
        ref = error_info.get("reference_preview", "")
        
        # Enhanced numerical diff feedback
        diff_details = _build_numerical_diff_details(error_info, config)
        
        return f"""C++ code for {config['name']} compiles and runs but gives WRONG numerical results.
Max relative error: {num.get('max_rel_diff','N/A')}
Reference (first 6 values): {ref}

{diff_details}

INSTRUCTIONS: Based on the error analysis above, identify and fix the bug.
Common causes: wrong sign, swapped indices, missing transpose, incorrect DOF mapping,
row-major vs column-major confusion, missing scale factor, wrong formula coefficient.
Code:\n```cpp\n{code}\n```\nReturn COMPLETE fixed file."""


def _build_numerical_diff_details(error_info: dict, config: dict) -> str:
    """Build detailed numerical comparison for repair feedback."""
    num = error_info.get("numerical", {})
    computed_arr = error_info.get("computed_array")
    reference_arr = error_info.get("reference_array")
    
    if computed_arr is None or reference_arr is None:
        return "No detailed comparison available."
    
    lines = ["DETAILED ERROR ANALYSIS:"]
    
    c = np.asarray(computed_arr, dtype=float)
    r = np.asarray(reference_arr, dtype=float)
    
    sz = config.get("matrix_size")
    
    # Reshape to matrix if applicable
    if sz and c.ndim == 1 and len(c) >= sz * sz:
        c = c[:sz*sz].reshape(sz, sz)
    if sz and r.ndim == 1 and len(r) >= sz * sz:
        r = r[:sz*sz].reshape(sz, sz)
    
    if c.shape != r.shape:
        lines.append(f"  Shape mismatch: computed {c.shape} vs reference {r.shape}")
        return "\n".join(lines)
    
    # Check for NaN/Inf
    if np.any(np.isnan(c)):
        nan_count = np.sum(np.isnan(c))
        lines.append(f"  WARNING: Output contains {nan_count} NaN values!")
        lines.append("  This usually means division by zero or singular matrix inversion.")
        lines.append("  Check all matrix inversions and divisions for degenerate cases.")
        return "\n".join(lines)
    
    if np.any(np.isinf(c)):
        lines.append("  WARNING: Output contains Inf values! Check for overflow or division by zero.")
        return "\n".join(lines)
    
    # Check all-zeros
    if np.allclose(c, 0, atol=1e-20):
        lines.append("  OUTPUT IS ALL ZEROS! The computation produces no meaningful result.")
        lines.append("  Check: is the main computation actually being called? Are results written to output?")
        return "\n".join(lines)
    
    # Element-wise comparison for matrices
    if c.ndim == 2:
        abs_diff = np.abs(c - r)
        max_ref = max(np.max(np.abs(r)), 1e-15)
        
        # Find top 5 worst elements
        flat_diff = abs_diff.flatten()
        worst_indices = np.argsort(flat_diff)[-5:][::-1]

        # Baseline error, needed by several pattern checks below. Previously
        # this was assigned only inside the square-matrix transpose check but
        # read again in the DOF-swap check, which is guarded on row count
        # alone. Any non-square output with >=6 rows -- e.g. B12's 11x3 time
        # history -- therefore raised UnboundLocalError, so those benchmarks
        # could never receive a numerical repair at all.
        orig_err = float(np.max(flat_diff) / max_ref)
        
        lines.append(f"  Matrix size: {c.shape[0]}x{c.shape[1]}")
        lines.append(f"  Overall max relative error: {np.max(flat_diff)/max_ref:.4e}")
        lines.append(f"  Worst elements (row, col): computed vs expected:")
        
        for idx in worst_indices:
            row, col = idx // c.shape[1], idx % c.shape[1]
            cv = c[row, col]
            rv = r[row, col]
            if abs(rv) > 1e-10:
                rel = abs(cv - rv) / abs(rv)
                lines.append(f"    K[{row},{col}]: {cv:.6e} vs {rv:.6e} (rel_err={rel:.2e})")
            else:
                lines.append(f"    K[{row},{col}]: {cv:.6e} vs {rv:.6e} (ref~0)")
        
        # Pattern detection: check for transpose
        if c.shape[0] == c.shape[1]:
            trans_err = np.max(np.abs(c - r.T)) / max_ref
            if trans_err < orig_err * 0.3:
                lines.append(f"\n  PATTERN: Matrix may be TRANSPOSED (error with .T = {trans_err:.2e} vs {orig_err:.2e})")
                lines.append("  FIX: Check if you need to transpose the result, or swap row/column assembly.")
        
        # Pattern detection: check for sign error
        sign_match = np.sign(c) == np.sign(r)
        sign_mismatch_rate = 1.0 - np.mean(sign_match[np.abs(r) > 1e-10])
        if sign_mismatch_rate > 0.3:
            lines.append(f"\n  PATTERN: {sign_mismatch_rate*100:.0f}% of elements have WRONG SIGN")
            lines.append("  FIX: Check for missing negative sign in coordinate transform or force direction.")
        
        # Pattern detection: check for uniform scale factor
        nonzero_mask = np.abs(r) > 1e-10
        if np.sum(nonzero_mask) > 5:
            ratios = c[nonzero_mask] / r[nonzero_mask]
            mean_ratio = np.mean(ratios)
            std_ratio = np.std(ratios)
            if std_ratio / (abs(mean_ratio) + 1e-15) < 0.1 and abs(mean_ratio - 1.0) > 0.01:
                lines.append(f"\n  PATTERN: Uniform SCALE FACTOR detected: computed/reference = {mean_ratio:.4f}")
                lines.append(f"  FIX: You may be missing a multiplication by {1.0/mean_ratio:.4f}")
        
        # Near-miss with directional structure.
        #
        # None of the detectors above fire on a small structured error: the
        # scale detector needs |mean ratio - 1| > 0.01 and this sits at 5e-4,
        # while the worst-element list shows cells that all look almost right.
        # Two B4 repeats stalled at rel 2.74e-03 for five iterations with no
        # actionable feedback, even though the deviation was perfectly ordered:
        # every x DOF 0.22% high, every y DOF 0.21% low. Reporting per-direction
        # means turns that into something a model can chase, and says plainly
        # that the structure is already right so it should not rewrite it.
        dof = config.get("dof_per_node")
        if (nonzero_mask.sum() > 5 and 1e-9 < orig_err < 1e-2
                and isinstance(dof, int) and 2 <= dof <= 6
                and c.shape[0] % dof == 0 and c.shape[1] % dof == 0):
            ratio_full = np.where(nonzero_mask, c / np.where(nonzero_mask, r, 1),
                                  np.nan)
            names = {2: ["x", "y"], 3: ["x", "y", "z"],
                     6: ["x", "y", "z", "rx", "ry", "rz"]}.get(
                         dof, [f"d{k}" for k in range(dof)])
            # Group by the direction pair (row DOF, column DOF). Averaging a
            # whole row instead mixes the like-direction and cross-direction
            # terms and cancels the very structure worth reporting: on the
            # stalled B4 artifact the xx cells were +0.22% and the yy cells
            # -0.21%, which a row-wise mean flattened to +0.05%.
            grid = {}
            for a in range(dof):
                for b in range(dof):
                    sub = ratio_full[a::dof, b::dof]
                    fin = sub[np.isfinite(sub)]
                    if fin.size:
                        grid[(a, b)] = (float(fin.mean()), float(fin.std()))
            if grid:
                lines.append(
                    f"\n  PATTERN: NEAR MISS. The overall error is "
                    f"{orig_err:.2e}, so the structure, the assembly and the "
                    f"sign conventions are already right. Do NOT rewrite the "
                    f"algorithm; find the single coefficient that is slightly "
                    f"off.")
                lines.append("  Mean computed/reference ratio by DOF-direction "
                             "pair:")
                for (a, b), (mean, std) in sorted(grid.items()):
                    tag = "same" if a == b else "cross"
                    lines.append(
                        f"    {names[a]}-{names[b]:<3s} ({tag:5s}): "
                        f"{mean:.6f}  ({(mean-1)*100:+.3f}%)")
                same = [grid[(k, k)][0] for k in range(dof) if (k, k) in grid]
                spread = max(same) - min(same) if len(same) > 1 else 0.0
                if spread > 1e-3:
                    lines.append(
                        f"  The like-direction terms disagree with each other "
                        f"by {spread*100:.2f} percentage points, so the error "
                        f"is DIRECTION DEPENDENT. Look at whatever treats the "
                        f"axes asymmetrically: the b and c coefficients of the "
                        f"strain-displacement matrix, the per-sub-element "
                        f"area, or the row order of the shear term.")
                else:
                    lines.append(
                        "  The like-direction terms agree, so the error is "
                        "direction independent. Look for one scalar: a "
                        "material constant, the thickness, or an area factor.")

        # Pattern detection: check for DOF swap (columns or rows swapped).
        # Only meaningful for square matrices: the check permutes rows AND the
        # corresponding columns, which requires both axes to index the same DOF.
        if c.shape[0] == c.shape[1] and c.shape[0] >= 6:
            best_swap = None
            best_swap_err = orig_err
            for i in range(min(c.shape[0], 12)):
                for j in range(i+1, min(c.shape[0], 12)):
                    # Check if swapping rows/cols i and j reduces error
                    c_swap = c.copy()
                    c_swap[[i,j]] = c_swap[[j,i]]
                    c_swap[:,[i,j]] = c_swap[:,[j,i]]
                    swap_err = np.max(np.abs(c_swap - r)) / max_ref
                    if swap_err < best_swap_err * 0.3:
                        best_swap = (i, j, swap_err)
                        best_swap_err = swap_err
            if best_swap:
                i, j, swap_err = best_swap
                lines.append(f"\n  PATTERN: Swapping DOF indices {i} and {j} significantly reduces error ({swap_err:.2e} vs {orig_err:.2e})")
                lines.append(f"  FIX: Check if DOF {i} and DOF {j} are assigned to wrong positions in assembly.")
    
    else:
        # Vector comparison
        lines.append(f"  Vector length: {len(c)}")
        abs_diff = np.abs(c - r)
        worst_idx = np.argmax(abs_diff)
        lines.append(f"  Worst element [{worst_idx}]: computed={c[worst_idx]:.6e} vs ref={r[worst_idx]:.6e}")
    
    return "\n".join(lines)

# ============================================================
# Code Generation & Extraction
# ============================================================

def extract_code(response: str) -> str:
    """Extract C++ code from LLM response."""
    for marker in ["```cpp", "```c++", "```c"]:
        if marker in response:
            return response.split(marker, 1)[1].split("```")[0].strip()
    if "```" in response:
        parts = response.split("```", 2)
        if len(parts) >= 3:
            code = parts[1]
            nl = code.find("\n")
            if 0 < nl < 20:
                code = code[nl+1:]
            return code.strip()
    return response.strip()


def generate_code(model: str, method: str, fortran_code: str,
                  config: dict, logger) -> str:
    """Generate C++ via direct prompt, tool agent, or staged pipeline."""
    llm = LLMClient(model=model, logger=logger)
    prompt = build_generation_prompt(config, fortran_code)
    rq3_cfg = rq3_method_config(method)
    base_method = rq3_cfg["base_method"] if rq3_cfg else method

    if base_method == "direct":
        resp = llm.generate(
            system_prompt="You are an expert in Fortran 77, C++17/Eigen, and FEA.",
            user_prompt=prompt, agent_name=f"direct_{config['id']}")
        if not resp:
            return ""
        return extract_code(resp)

    if base_method == "physics_first":
        return generate_code_physics_first(model, fortran_code, config, logger)

    if base_method == "agent":
        from tool_augmented_agent import run_tool_augmented
        agent_options = {}
        if rq3_cfg:
            agent_options = {
                "enabled_tools": rq3_cfg["enabled_tools"],
                "knowledge_mode": "off",
                "enable_checkpoints": False,
                "max_steps": 15,
                "max_tokens": 100_000,
            }
        return run_tool_augmented(model, config, fortran_code, logger, llm,
                                  **agent_options)

    if base_method == "phased_agent":
        from phased_agent import run_phased_agent
        return run_phased_agent(model, config, fortran_code, logger, llm)

    # Multi-agent capabilities are five INDEPENDENT axes (see MA_AXES). The
    # named methods are presets over those axes, not separate pipelines.
    flags = resolve_ma_flags(method)
    enable_mining = flags["mining"]
    use_kb = flags["kb"] and _knowledge_base is not None

    # 'wide' widens the channel from the pipeline to the generator.
    #
    # Measured on the archive: the orchestrator's agents emit a median 199,825
    # characters, of which only json.dumps(modern_design)[:2500] could ever
    # reach the C++ generator - at most 1.3%, and as paraphrase rather than
    # source. physics_spec and architecture never reached it at all. So
    # `multiagent` was structurally close to `direct` plus a 2.5 KB appendix,
    # which is the mechanism behind the two arms being indistinguishable.
    DESIGN_BUDGET = 6000 if flags["wide"] else 2500
    SPEC_BUDGET = 3000 if flags["wide"] else 0
    ARCH_BUDGET = 1500 if flags["wide"] else 0

    logger.log_info(f"  MA capabilities: {format_ma_flags(flags)}")

    try:
        nastran_src = str(PROJECT_ROOT / "nastran" / "NASTRAN-95" / "mis")
        orch = RefactoringOrchestrator(
            llm_client=llm, logger=logger,
            output_dir=str(PROJECT_ROOT / "results" / logger.experiment_id),
            enable_source_mining=enable_mining,
            knowledge_base=_knowledge_base if use_kb else None,
            nastran_source_dir=nastran_src if enable_mining else None,
            skip_unused_stages=flags["lean"])
        state = orch.run(
            source_files={f"{config['name'].lower()}.f": fortran_code},
            context=(_orchestrator_context(config, fortran_code)
                     if flags["srcctx"]
                     else f"{config['description']}. Target: C++17/Eigen."),
            reference_results=[],
            benchmark_config=(config if (use_kb or enable_mining
                                         or flags["srcctx"]) else None))

        # RefactoringOrchestrator records stage exceptions in state.errors and
        # may return a partial state instead of raising. Such a run did not
        # execute the declared Multi-Agent method and must not enter its cell.
        fatal_stages = {"architecture", "physics_spec", "modern_design"}
        if not enable_mining:
            fatal_stages.add("comprehension")
        fatal_errors = [
            error for error in (getattr(state, "errors", None) or [])
            if error.get("stage") in fatal_stages
        ]
        if fatal_errors or not getattr(state, "modern_design", None):
            failed = ", ".join(
                sorted({error.get("stage", "unknown")
                        for error in fatal_errors}))
            logger.log_info(
                "  Multi-agent required stages incomplete"
                + (f": {failed}" if failed else ": modern_design missing"))
            return "ERROR: MULTI_AGENT_PIPELINE"

        if state and hasattr(state, 'modern_design') and state.modern_design:
            design_str = json.dumps(state.modern_design, indent=2,
                                    default=str)[:DESIGN_BUDGET]
            prompt += f"\n\nDesign context:\n{design_str}"

        if flags["selfmine"]:
            conv = _extract_entry_conventions(fortran_code)
            block = _format_entry_conventions(conv)
            if block:
                prompt += "\n" + block
                kinds = {}
                for c in conv:
                    kinds[c["kind"]] = kinds.get(c["kind"], 0) + 1
                logger.log_info(f"  selfmine: {len(conv)} conventions "
                                f"from the entry routine {kinds}")

        if flags["wide"]:
            spec = getattr(state, "physics_spec", None)
            if spec:
                why = _spec_contradicts_harness(spec, config)
                if why:
                    # A specification that disagrees with the fixed interface
                    # is worse than no specification: it instructs the
                    # generator to build the wrong thing.
                    logger.log_info(
                        f"  Dropping recovered spec, contradicts harness: {why}")
                else:
                    prompt += ("\n\nRecovered specification (from the Fortran, "
                               "not from a textbook):\n"
                               + json.dumps(spec, indent=2,
                                            default=str)[:SPEC_BUDGET])
            arch = getattr(state, "architecture", None)
            if arch and arch.get("call_graph"):
                prompt += ("\n\nCall graph recovered from the source:\n"
                           + json.dumps(arch.get("call_graph"), indent=2,
                                        default=str)[:ARCH_BUDGET])
        
        # v2: Inject source mining conventions into the generation prompt
        if enable_mining and state and state.source_mining:
            mining = state.source_mining
            conventions = mining.get("critical_conventions", [])
            storage_warnings = mining.get("storage_warnings", [])
            if conventions or storage_warnings:
                prompt += "\n\n⚠️ CRITICAL CONVENTIONS (from called subroutine analysis):\n"
                for w in storage_warnings:
                    prompt += f"- STORAGE: {w}\n"
                for c in conventions:
                    prompt += f"- {c}\n"
                prompt += "You MUST respect these conventions for numerical correctness.\n"
        
        # v2: Inject KB context into generation prompt
        if use_kb and state and state.kb_context:
            prompt += f"\n\n{state.kb_context}"
    
    except Exception as e:
        # A Multi-Agent infrastructure failure must never be relabelled as a
        # Direct result. Return an explicit marker so the repetition is
        # recorded as invalid and re-queued by the caller.
        logger.log_info(f"  Multi-agent pipeline error: {e}. Marking repetition invalid.")
        return "ERROR: MULTI_AGENT_PIPELINE"

    resp = llm.generate(
        system_prompt="You are an expert C++17/Eigen developer for FEA.",
        user_prompt=prompt, agent_name=f"codegen_{config['id']}")
    if not resp:
        return ""
    return extract_code(resp)


def generate_code_physics_first(model: str, fortran_code: str,
                                config: dict, logger) -> str:
    """
    Physics-First strategy (EAAI-inspired):
    1. Extract physics specification from Fortran source
    2. Extract numerical method specification
    3. Generate C++ from specs ONLY (without seeing Fortran source)
    
    This separates "understanding legacy code" from "writing new code",
    preventing knowledge conflation.
    """
    llm = LLMClient(model=model, logger=logger)
    test_params = _TEST_PARAMS.get(config["id"], "")
    
    # Stage 1: Extract physics specification from Fortran
    physics_prompt = f"""Analyze this NASTRAN-95 Fortran 77 subroutine and extract its PHYSICS SPECIFICATION.

Do NOT translate the code. Instead, document WHAT it computes in mathematical terms:
- Governing equations (e.g., K = integral of B^T * D * B)
- Constitutive relations (stress-strain law)
- Degrees of freedom per node
- Coordinate transformations used
- Any special conditions (pin flags, offsets, modes)

Fortran source ({config['name']}, {config['total_lines']} lines):
```fortran
{fortran_code}
```

Output a structured physics specification in this format:
ELEMENT: [name]
TYPE: [rod/beam/membrane/plate/solid/etc]
NODES: [number] nodes, [DOF] DOF per node
MATRIX SIZE: [N]x[N]
GOVERNING EQUATION: [mathematical formula for stiffness matrix]
MATERIAL MODEL: [isotropic/anisotropic, parameters needed]
COORDINATE TRANSFORM: [description if any]
SPECIAL FEATURES: [any non-standard logic]
KEY FORMULAS: [list the actual mathematical expressions used]
"""
    
    logger.log_info("  [Physics-First] Stage 1: Extracting physics spec...")
    physics_spec = llm.generate(
        "You are an expert in structural mechanics and FEM theory.",
        physics_prompt, f"physics_extract_{config['id']}")
    
    if not physics_spec:
        logger.log_info("  [Physics-First] Physics extraction failed (API issue)")
        return ""
    
    # Stage 2: Extract numerical method details
    numerics_prompt = f"""Based on this Fortran subroutine, extract the NUMERICAL METHOD details:
- Integration scheme (analytical, Gauss quadrature, etc.)
- Matrix assembly strategy (direct formula, B-matrix approach, etc.)
- Any numerical tricks (reduced integration, selective, etc.)
- Order of operations

Fortran source ({config['name']}):
```fortran
{fortran_code[:3000]}
```

Output the numerical algorithm as pseudocode (NOT Fortran, NOT C++).
"""
    
    logger.log_info("  [Physics-First] Stage 2: Extracting numerical method...")
    numerics_spec = llm.generate(
        "You are an expert in computational mechanics numerical methods.",
        numerics_prompt, f"numerics_extract_{config['id']}")
    
    if not numerics_spec:
        numerics_spec = "(numerical method extraction failed)"
    
    # Stage 3: Generate C++ from specs ONLY (no Fortran source!)
    codegen_prompt = f"""Implement the following finite element in C++17 using Eigen.

PHYSICS SPECIFICATION:
{physics_spec[:3000]}

NUMERICAL METHOD:
{numerics_spec[:2000]}

{test_params}

Requirements:
1. Use Eigen library (#include <Eigen/Dense>)
2. Single self-contained file with main() that runs the test case
3. Compilable with MSVC: cl.exe /std:c++17 /EHsc /O2
4. Double precision throughout

Output specification:
{_output_spec(config)}

IMPORTANT: Implement from the physics/math specification above.
Do NOT try to translate Fortran line-by-line.
Use clean modern C++ with descriptive variable names.

Generate COMPLETE compilable C++ code.
"""
    
    logger.log_info("  [Physics-First] Stage 3: Generating C++ from specs...")
    resp = llm.generate(
        "You are an expert C++17/Eigen developer specializing in FEM implementation.",
        codegen_prompt, f"codegen_physics_{config['id']}")
    
    if not resp:
        return ""
    return extract_code(resp)

# ============================================================
# Single Experiment Runner (compile-run-verify with repair loop)
# ============================================================

# Number of consecutive empty API responses tolerated within one repetition
# before it is abandoned as invalid. Raised from 5 after the contamination
# audit showed abandoned repetitions being scored as model failures.
API_RETRY_BUDGET = 8

# How many times a repetition abandoned for API reasons may be re-queued
# before it is recorded as invalid and left out of the pass rate.
API_REQUEUE_LIMIT = 2


def _api_backoff(attempt: int, logger, base: int = 15, cap: int = 240) -> None:
    """Exponential backoff between API retries, capped."""
    import time
    delay = min(base * (2 ** max(attempt - 1, 0)), cap)
    logger.log_info(f"    backing off {delay}s before retry")
    time.sleep(delay)


def run_single_experiment(benchmark_id: str, model: str, method: str,
                          repeat_num: int, logger,
                          max_repair_iterations: Optional[int] = None) -> dict:
    """Full experiment iteration: generate, compile, run, verify, repair."""
    method_cfg = rq3_method_config(method)
    if max_repair_iterations is None:
        max_repair_iterations = (method_cfg["max_outer_repairs"]
                                 if method_cfg else MAX_REPAIR_ITERATIONS)
    if max_repair_iterations < 0:
        raise ValueError("max_repair_iterations must be non-negative")
    config = load_benchmark_config(benchmark_id)
    fortran_code = read_fortran_source(config)
    if not fortran_code.strip():
        return {"passed": False, "benchmark": benchmark_id, "model": model,
                "method": method, "repeat": repeat_num,
                "error": "Source not found", "iterations": 0}

    reference = load_reference(config)
    logger.log_info(f"--- {benchmark_id} {config['name']} | {model} | {method} | run {repeat_num} ---")
    logger.log_stage_start(f"{benchmark_id}_{method}_r{repeat_num}")

    code = generate_code(model, method, fortran_code, config, logger)
    if not code or len(code) < 50 or code.startswith("ERROR:"):
        # Retry once if initial generation or the selected method pipeline fails.
        import time; time.sleep(10)
        logger.log_info("  Initial generation failed, retrying after 10s...")
        code = generate_code(model, method, fortran_code, config, logger)
    if not code or len(code) < 50 or code.startswith("ERROR:"):
        method_failed = bool(code and code.startswith("ERROR: MULTI_AGENT_PIPELINE"))
        abort_reason = "method_pipeline_error" if method_failed else "api_exhausted"
        logger.log_info(
            f"  ABORT[{abort_reason}]: initial generation never returned "
            "an output from the requested method")
        logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", False)
        return _aborted(benchmark_id, model, method, repeat_num, abort_reason,
                        error=("Multi-agent pipeline failed"
                               if method_failed
                               else "Code generation failed (API issue)"),
                        iterations=0)

    llm = LLMClient(model=model, logger=logger)
    iteration = 0
    api_failures = 0
    max_api_failures = API_RETRY_BUDGET
    
    # Regression prevention: track best code and its error
    best_code = code
    best_rel_diff = float('inf')
    best_compile_ok = False
    best_run_ok = False

    while iteration <= max_repair_iterations:
        logger.log_info(f"  iter {iteration}: compile & run")
        logger.log_artifact(f"code_{benchmark_id}_r{repeat_num}_i{iteration}.cpp", code, is_json=False)

        with tempfile.TemporaryDirectory(prefix="v3_") as td:
            result = compile_and_run_cpp(code, Path(td))

        if not result["success"]:
            if iteration >= max_repair_iterations:
                logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", False)
                return _result(benchmark_id, model, method, repeat_num, False,
                               compile_ok=(result["stage"] != "compile"),
                               iterations=iteration,
                               error=f"{result['stage']}: {result['error'][:150]}")
            # Repair from the best code that is known to compile and run, not
            # from whatever the last iteration left behind. Regression
            # prevention used to fire only when a candidate compiled, ran, and
            # scored worse; a candidate that failed to compile never reached
            # that branch, so one bad iteration poisoned all the rest.
            repair_base = code
            if best_compile_ok and best_run_ok and best_code is not code \
                    and _is_plausible_program(best_code) is None:
                logger.log_info("    repairing from the best code that ran, "
                                f"not the failing candidate "
                                f"(stage={result['stage']})")
                repair_base = best_code
                code = best_code
            rp = build_repair_prompt(config, repair_base, result)
            resp = llm.generate("Fix the C++ code. Return COMPLETE file.", rp,
                                f"repair_{benchmark_id}_i{iteration+1}")
            if not resp or len(resp.strip()) < 30 or resp.startswith("ERROR:"):
                # API failure - don't consume iteration
                api_failures += 1
                logger.log_info(f"  API failure {api_failures}/{max_api_failures} (not counting as iteration)")
                if api_failures >= max_api_failures:
                    logger.log_info(f"  ABORT[api_exhausted]: no candidate evaluated "
                                    f"after {api_failures} empty responses")
                    logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", False)
                    return _aborted(benchmark_id, model, method, repeat_num,
                                    "api_exhausted",
                                    compile_ok=(result["stage"] != "compile"),
                                    error="API failure (repeated timeouts)",
                                    iterations=iteration)
                _api_backoff(api_failures, logger)
                continue
            nc = extract_code(resp)
            code = _accept_candidate(nc, code, best_code, logger,
                                     f"repair_i{iteration+1}")
            if _is_plausible_program(nc) is None:
                api_failures = 0  # Reset on success
            iteration += 1
            continue

        # Verify
        computed = parse_output(result.get("stdout", ""), config)
        vr = verify_output(computed, reference, config)
        if vr.get("passed"):
            logger.log_info(f"  PASS rel_diff={vr.get('max_rel_diff',0):.2e}")
            logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", True)
            
            # v2: Accumulate knowledge on success
            if _knowledge_base and method in ("multiagent_v2", "multiagent_v3",
                                              "multiagent_kb"):
                try:
                    _knowledge_base.add_from_successful_run(
                        benchmark_id=benchmark_id,
                        subroutine_name=config.get("name", ""),
                        model=model,
                        method=method,
                        category=config.get("category", ""),
                        called_subroutines=_extract_calls_from_fortran(fortran_code),
                        final_code=code,
                    )
                    logger.log_info(f"  KB: Accumulated insights from successful run")
                except Exception as e:
                    logger.log_info(f"  KB: Failed to accumulate ({e})")
            
            return _result(benchmark_id, model, method, repeat_num, True,
                           compile_ok=True, run_ok=True, numerical_ok=True,
                           iterations=iteration, max_rel_diff=vr.get("max_rel_diff", 0))

        if iteration >= max_repair_iterations:
            # Report using best code metrics (regression prevention). A
            # contract violation carries no relative difference, so treat it
            # as +inf instead of comparing None with a float.
            latest_diff = vr.get("max_rel_diff")
            if latest_diff is None:
                latest_diff = float('inf')
            final_diff = min(best_rel_diff, latest_diff)
            logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", False)
            return _result(benchmark_id, model, method, repeat_num, False,
                           compile_ok=True, run_ok=True, numerical_ok=False,
                           iterations=iteration, max_rel_diff=final_diff)

        # Regression prevention: track best numerical result.
        # A contract violation carries no rel_diff, so normalise to +inf --
        # comparing None against a float raises TypeError in Python 3.
        current_diff = vr.get("max_rel_diff")
        if current_diff is None:
            current_diff = float('inf')
        if current_diff < best_rel_diff:
            best_rel_diff = current_diff
            best_code = code
            best_compile_ok = True
            best_run_ok = True
            logger.log_info(f"    new best: rel_diff={current_diff:.2e}")
        elif iteration > 0 and current_diff > best_rel_diff * 2:
            # Regression detected: revert to best code for next repair attempt
            logger.log_info(f"    regression detected (rel={current_diff:.2e} > best={best_rel_diff:.2e}), reverting to best code")
            code = best_code

        # Decide repair strategy: regular repair vs root cause analysis
        # Trigger root cause analysis if stuck (best hasn't improved after 2+ numerical attempts)
        numerical_attempts = iteration  # approximate: how many iterations we've done
        use_root_cause = (
            numerical_attempts >= 2 and 
            best_rel_diff < float('inf') and
            best_rel_diff > 0.001 and  # still failing (not just floating point noise)
            computed is not None and reference is not None
        )
        
        if use_root_cause:
            # Root cause analysis: deeper LLM-driven debugging
            resp = _run_root_cause_analysis(llm, config, best_code, computed, reference, logger)
            if resp and len(resp.strip()) > 50:
                nc = extract_code(resp)
                code = _accept_candidate(nc, code, best_code, logger,
                                         f"root_cause_i{iteration+1}")
                if _is_plausible_program(nc) is None:
                    api_failures = 0
            iteration += 1
            continue

        # Regular numerical repair.
        # A contract violation is reported as such rather than as a numerical
        # deviation: the model must be told the SHAPE is wrong, otherwise it
        # searches for an arithmetic bug that does not exist.
        ref_preview = str(reference.flatten()[:6].tolist()) if reference is not None else "N/A"
        if vr.get("contract_violation"):
            logger.log_info(f"    contract violation: {vr.get('reason','')}")
            ei = {"stage": "contract",
                  "error": vr.get("reason", "output shape mismatch"),
                  "numerical": vr, "reference_preview": ref_preview,
                  "computed_array": None, "reference_array": None}
            tag = f"repair_contract_{benchmark_id}_i{iteration+1}"
            sysmsg = ("The output SHAPE is wrong, not the arithmetic. "
                      "Emit exactly the required rows and columns.")
        else:
            ei = {"stage": "numerical", "error": str(vr.get("max_rel_diff")),
                  "numerical": vr, "reference_preview": ref_preview,
                  "computed_array": computed.tolist() if computed is not None else None,
                  "reference_array": reference.tolist() if reference is not None else None}
            tag = f"repair_num_{benchmark_id}_i{iteration+1}"
            sysmsg = "Fix numerical error."
        rp = build_repair_prompt(config, code, ei)
        resp = llm.generate(sysmsg, rp, tag)
        if not resp or len(resp.strip()) < 30 or resp.startswith("ERROR:"):
            # API failure - don't consume iteration
            api_failures += 1
            logger.log_info(f"  API failure {api_failures}/{max_api_failures} (not counting as iteration)")
            if api_failures >= max_api_failures:
                # NOTE: unlike the compile-stage abort, a candidate WAS
                # evaluated here (it compiled, ran, and was compared), so the
                # best observed rel_diff is retained for diagnostics even
                # though the repetition itself is marked invalid.
                logger.log_info(f"  ABORT[api_exhausted]: repair unavailable "
                                f"after {api_failures} empty responses")
                logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", False)
                return _aborted(benchmark_id, model, method, repeat_num,
                                "api_exhausted",
                                compile_ok=True, run_ok=True, numerical_ok=False,
                                max_rel_diff=(None if best_rel_diff == float('inf')
                                              else best_rel_diff),
                                error="API failure (repeated timeouts)",
                                iterations=iteration)
            _api_backoff(api_failures, logger)
            continue
        nc = extract_code(resp)
        code = _accept_candidate(nc, code, best_code, logger,
                                 f"numeric_repair_i{iteration+1}")
        if _is_plausible_program(nc) is None:
            api_failures = 0
        iteration += 1

    logger.log_stage_end(f"{benchmark_id}_{method}_r{repeat_num}", False)
    return _result(benchmark_id, model, method, repeat_num, False,
                   iterations=max_repair_iterations, error="Unexpected exit")


def _result(bid, model, method, repeat, passed, **kw) -> dict:
    """
    Helper to build result dict.

    `valid` distinguishes an evaluated outcome from an aborted one. A run is
    INVALID when the harness never managed to evaluate a candidate
    translation -- currently only when the model API returns empty responses
    until the retry budget is exhausted. Invalid runs must be excluded from
    pass rates rather than counted as failures: doing otherwise attributes an
    infrastructure problem to the model. See analyze_api_contamination.py,
    which found 21 of 289 historical failures to be of this kind, unevenly
    distributed across backends (17.0% of one backend's failures against 0%
    of another's), which biases backend comparison specifically.
    """
    r = {"passed": passed, "valid": True, "abort_reason": None,
         "benchmark": bid, "model": model, "method": method, "repeat": repeat}
    r.update(kw)
    return r


def _aborted(bid, model, method, repeat, reason, **kw) -> dict:
    """Build a result marking the repetition as not evaluated (invalid)."""
    r = _result(bid, model, method, repeat, False, **kw)
    r["valid"] = False
    r["abort_reason"] = reason
    return r


def _extract_calls_from_fortran(fortran_code: str) -> list[str]:
    """Extract called subroutine names from Fortran source."""
    calls = []
    for match in re.finditer(r'^\s*CALL\s+(\w+)', fortran_code, re.MULTILINE | re.IGNORECASE):
        name = match.group(1).upper()
        if name not in calls:
            calls.append(name)
    return calls

# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="NASTRAN-95 Unified Benchmark Runner v3")
    parser.add_argument("--benchmark", default="all", help="all | B1 | B1,B5,B10")
    parser.add_argument("--model", default="deepseek-v4-flash", help="all | model name")
    parser.add_argument("--method", default="direct",
                        help="all | direct | agent | physics_first | "
                             + " | ".join(MA_PRESETS)
                             + " | any composition of modules, e.g. "
                             + "multiagent_mining+wide  (modules: "
                             + ", ".join(MA_CODES) + ")")
    parser.add_argument("--ma-axes", action="store_true",
                        help="list the multi-agent capability axes and the "
                             "presets over them, then exit")
    parser.add_argument("--ma-factorial", metavar="CODES",
                        help="run the full factorial over the given comma-"
                             "separated axes, e.g. --ma-factorial "
                             "mining,wide,srcctx runs 2^3 = 8 cells. Other "
                             "axes stay off unless named in --ma-base")
    parser.add_argument("--ma-base", metavar="CODES", default="",
                        help="axes held ON across every factorial cell")
    parser.add_argument("--lang", default="cpp", help="Target language (cpp)")
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--seed-kb", action="store_true",
                        help="Seed knowledge base with Level 3 insights before running")
    parser.add_argument("--clear-kb", action="store_true",
                        help="Clear knowledge base before running")
    parser.add_argument("--kb-path", default=None,
                        help="Path to knowledge base JSON file")
    args = parser.parse_args()

    # Initialize knowledge base for v2 methods
    global _knowledge_base
    if args.ma_axes:
        print("Multi-agent capability axes\n" + "=" * 70)
        for code, desc in MA_AXES:
            print(f"\n  {code}")
            for line in textwrap.wrap(desc, 66):
                print(f"      {line}")
        print("\n\nPresets\n" + "=" * 70)
        print(f"  {'preset':18s} " + " ".join(f"{c:>7s}" for c in MA_CODES))
        for name in sorted(MA_PRESETS, key=lambda n: len(MA_PRESETS[n])):
            f = resolve_ma_flags(name)
            print(f"  {name:18s} "
                  + " ".join(f"{'on' if f[c] else '-':>7s}" for c in MA_CODES))
        print(f"\n  Any other combination is addressable by composing module")
        print(f"  names with '+', in any order:")
        for ex in ("multiagent_mining+wide", "multiagent_lean",
                   "multiagent_srcctx+kb"):
            print(f"    {ex:34s} {format_ma_flags(resolve_ma_flags(ex))}")
        print(f"\n  Retired labels still accepted for archived runs:")
        for old, new in MA_ALIASES.items():
            print(f"    {old:18s} -> {new}")
        return

    needs_kb = (args.method == "all"
                or (args.method.startswith("multiagent")
                    and resolve_ma_flags(args.method)["kb"])
                or (args.ma_factorial and "kb" in
                    (args.ma_factorial + "," + args.ma_base)))
    if needs_kb or args.seed_kb or args.clear_kb:
        kb_path = args.kb_path or str(PROJECT_ROOT / "knowledge_base.json")
        _knowledge_base = KnowledgeBase(storage_path=kb_path)
        
        if args.clear_kb:
            _knowledge_base.clear()
            print("Knowledge base cleared.")
        
        if args.seed_kb:
            _knowledge_base.seed_from_level3()
            stats = _knowledge_base.get_stats()
            print(f"Knowledge base seeded: {stats['total']} insights")
        
        stats = _knowledge_base.get_stats()
        print(f"Knowledge Base: {stats.get('total', 0)} insights loaded")

    full_cfg = load_all_benchmarks()
    all_ids = [b["id"] for b in full_cfg["benchmarks"]]

    # Resolve selections
    bids = all_ids if args.benchmark == "all" else [x.strip() for x in args.benchmark.split(",")]
    models = full_cfg["models"] if args.model == "all" else [x.strip() for x in args.model.split(",")]
    if args.ma_factorial:
        axes = [c.strip() for c in args.ma_factorial.split(",") if c.strip()]
        base = {c.strip() for c in args.ma_base.split(",") if c.strip()}
        bad = [c for c in axes + sorted(base) if c not in MA_CODES]
        if bad:
            raise SystemExit(f"unknown axis {bad}; valid axes: {MA_CODES}")
        methods = []
        for mask in range(2 ** len(axes)):
            on = set(base) | {a for i, a in enumerate(axes)
                              if mask >> i & 1}
            methods.append(ma_label({c: (c in on) for c in MA_CODES}))
        seen = set()
        methods = [m for m in methods if not (m in seen or seen.add(m))]
        print(f"factorial over {axes} with base {sorted(base) or '[]'}"
              f" -> {len(methods)} cells")
    elif args.method == "all":
        methods = ["direct", "agent", "physics_first"] + list(MA_PRESETS)
    else:
        methods = [args.method]

    total_runs = len(bids) * len(models) * len(methods) * args.repeat
    print(f"{'='*70}\nNASTRAN-95 Unified Pipeline v3\n{'='*70}")
    print(f"Benchmarks: {bids}\nModels: {models}\nMethods: {methods}")
    print(f"Repeats: {args.repeat} | Total runs: {total_runs}")
    print(f"Time: {datetime.now().isoformat()}\n{'='*70}\n")

    all_results = []
    run_count = 0

    for bid in bids:
        bcfg = load_benchmark_config(bid)
        print(f"\n{'='*60}\n{bid}: {bcfg['name']} - {bcfg['description']}")
        print(f"  {bcfg['category']} | {bcfg['total_lines']} lines | "
              f"{bcfg['difficulty']} | {bcfg['verification_type']}\n{'='*60}")

        for mdl in models:
            for mth in methods:
                eid = create_experiment_id(f"v3_{mdl}_{mth}", bid)
                logger = ExperimentLogger(eid, base_dir=str(PROJECT_ROOT / "results"))

                # Snapshot the exact benchmark definition this run used.
                # The directory name encodes benchmark, model, method and
                # timestamp but nothing about source_files, reference_method or
                # tolerances, so two runs made either side of a config edit were
                # distinguishable only by their timestamps. That ambiguity is
                # what made earlier batches hard to pool; recording the config
                # removes the guesswork.
                logger.log_artifact("benchmark_config", {
                    "config": bcfg,
                    "effective_protocol": {
                        "requested_method": mth,
                        "model": mdl,
                        "max_repair_iterations": (
                            rq3_method_config(mth)["max_outer_repairs"]
                            if rq3_method_config(mth)
                            else MAX_REPAIR_ITERATIONS),
                    },
                    "ma_flags": (resolve_ma_flags(mth)
                                 if mth.startswith("multiagent") else None),
                    "source_file_line_counts": {
                        s: (len((PROJECT_ROOT / s).read_text(
                                encoding="utf-8", errors="replace").splitlines())
                            if (PROJECT_ROOT / s).exists() else None)
                        for s in bcfg.get("source_files", [])},
                })

                for ri in range(1, args.repeat + 1):
                    run_count += 1
                    # Re-queue a repetition abandoned for infrastructure
                    # reasons rather than recording it as a model failure.
                    requeueable_reasons = {"api_exhausted", "method_pipeline_error"}
                    r = run_single_experiment(bid, mdl, mth, ri, logger)
                    requeues = 0
                    while (not r.get("valid", True)
                           and r.get("abort_reason") in requeueable_reasons
                           and requeues < API_REQUEUE_LIMIT):
                        requeues += 1
                        print(f"  [{run_count}/{total_runs}] INVALID "
                              f"({r.get('abort_reason')}) {mdl} | {mth:10s} | "
                              f"r{ri}: re-queueing "
                              f"({requeues}/{API_REQUEUE_LIMIT})")
                        logger.log_info(f"  RE-QUEUE {requeues}/{API_REQUEUE_LIMIT} "
                                        f"for {bid} {mdl} {mth} r{ri}: "
                                        f"{r.get('abort_reason')}")
                        _api_backoff(requeues + 2, logger)
                        r = run_single_experiment(bid, mdl, mth, ri, logger)
                    r["api_requeues"] = requeues
                    # Persist the final result immediately. The gap audit reads
                    # these records preferentially, so an interrupted batch
                    # neither loses completed repeats nor counts invalid ones.
                    logger.log_artifact(
                        f"result_{bid}_{mth}_r{ri}", r)
                    all_results.append(r)

                    if not r.get("valid", True):
                        s = "INVALID"
                    else:
                        s = "PASS" if r["passed"] else "FAIL"
                    d = f"iter={r.get('iterations',0)}"
                    if r.get("max_rel_diff") is not None:
                        d += f" rel={r['max_rel_diff']:.2e}"
                    if not r.get("valid", True):
                        d += f" [{r.get('abort_reason')}]"
                    print(f"  [{run_count}/{total_runs}] {s} {mdl} | {mth:10s} | r{ri}: {d}")

                batch = all_results[-args.repeat:]
                n_valid = sum(1 for x in batch if x.get("valid", True))
                n_pass = sum(1 for x in batch if x["passed"])
                # Report the pass rate over evaluated repetitions only.
                logger.finalize(n_valid > 0 and n_pass == n_valid,
                                f"{bid} {mdl} {mth}: {n_pass}/{n_valid}"
                                + (f" ({args.repeat - n_valid} invalid)"
                                   if n_valid < args.repeat else ""))

    # Summary
    print(f"\n{'='*70}\nFINAL SUMMARY\n{'='*70}")
    print(f"{'Bench':<6}{'Model':<20}{'Method':<12}{'Compile':<9}{'Run':<7}"
          f"{'Verify':<8}{'Pass':<8}{'Inval'}")
    print("-" * 78)
    sr = sorted(all_results, key=lambda x: (x["benchmark"], x["model"], x["method"]))
    for (b, m, mt), grp in groupby(sr, key=lambda x: (x["benchmark"], x["model"], x["method"])):
        items = list(grp)
        # Rates are computed over EVALUATED repetitions. Invalid ones are
        # counted separately so a low rate can never be an artefact of the
        # harness silently converting aborts into failures.
        ev = [i for i in items if i.get("valid", True)]
        n = len(ev)
        inval = len(items) - n
        if n == 0:
            print(f"{b:<6}{m:<20}{mt:<12}{'--':<9}{'--':<7}{'--':<8}{'--':<8}{inval}")
            continue
        c = sum(1 for i in ev if i.get("compile_ok"))
        ru = sum(1 for i in ev if i.get("run_ok"))
        nu = sum(1 for i in ev if i.get("numerical_ok"))
        p = sum(1 for i in ev if i["passed"])
        print(f"{b:<6}{m:<20}{mt:<12}{f'{c}/{n}':<9}{f'{ru}/{n}':<7}"
              f"{f'{nu}/{n}':<8}{f'{p}/{n}':<8}{inval}")

    total = len(all_results)
    evaluated = [r for r in all_results if r.get("valid", True)]
    invalid = [r for r in all_results if not r.get("valid", True)]
    tp = sum(1 for r in evaluated if r["passed"])
    print("-" * 78)
    if evaluated:
        print(f"TOTAL: {len(evaluated)} evaluated, {tp} passed "
              f"({tp/len(evaluated)*100:.1f}%)")
    if invalid:
        by_reason = Counter(r.get("abort_reason") for r in invalid)
        print(f"INVALID (excluded from the rate): {len(invalid)} of {total} "
              f"-- {dict(by_reason)}")
        rq = sum(r.get("api_requeues", 0) for r in all_results)
        print(f"  repetitions re-queued before abandoning: {rq}")
        print(f"  NOTE: report these separately in the paper; do not fold "
              f"them into the pass rate.")

    # Save
    results_dir = PROJECT_ROOT / "results"
    results_dir.mkdir(exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_file = results_dir / f"v3_summary_{ts}.json"
    out_file.write_text(json.dumps({
        "timestamp": datetime.now().isoformat(),
        "config": {"benchmarks": bids, "models": models, "methods": methods,
                   "repeat": args.repeat, "max_repair": MAX_REPAIR_ITERATIONS},
        "results": all_results,
        "overall": {
            "total": total,
            "evaluated": len(evaluated),
            "invalid": len(invalid),
            "invalid_by_reason": dict(Counter(r.get("abort_reason") for r in invalid)),
            "pass": tp,
            "pass_rate_over_evaluated": (tp / len(evaluated)) if evaluated else None,
            "api_requeues": sum(r.get("api_requeues", 0) for r in all_results),
        },
    }, indent=2, default=str), encoding="utf-8")
    print(f"Saved: {out_file}")


if __name__ == "__main__":
    main()
