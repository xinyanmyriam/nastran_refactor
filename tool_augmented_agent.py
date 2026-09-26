"""
Tool-Augmented Agent for code translation.

Unlike the fixed multi-agent pipeline where humans pre-define the stage order,
this agent decides its own strategy: it can write code, compile it, compare
against the oracle, read specific source lines, or inspect conventions — in
whatever order it chooses.

The only variable vs Direct is tool access. Model, source, oracle, tolerance
are all the same.

Usage from run_experiments_v3.py:
    method = "agent"
"""
import json
import re
import sys
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))
from run_multilang_v2 import (
    compile_and_run_cpp, parse_matrix_from_output,
)


# ============================================================
# Tool definitions
# ============================================================

TOOL_ORDER = (
    "compile_and_run",
    "compare_with_reference",
    "read_source",
    "inspect_conventions",
    "verify_checkpoint",
)

TOOL_DETAILS = {
    "compile_and_run": (
        "[TOOL: compile_and_run]\n"
        "   Compile and run your C++ code. Put the COMPLETE code between "
        "```cpp and ``` fences before this line.\n"
        "   Returns: compilation result, execution result, stdout."
    ),
    "compare_with_reference": (
        "[TOOL: compare_with_reference]\n"
        "   Compare the most recent successful run output against the "
        "compiled Fortran oracle.\n"
        "   Returns: max_rel_diff, passed, and structured diagnostics.\n"
        "   You cannot see the reference values directly."
    ),
    "read_source": (
        "[TOOL: read_source FILENAME START_LINE END_LINE]\n"
        "   Re-read selected lines from the Fortran source files.\n"
        "   Example: [TOOL: read_source kqdmem.f 110 120]"
    ),
    "inspect_conventions": (
        "[TOOL: inspect_conventions]\n"
        "   Return mechanically extracted implementation conventions of "
        "the routine."
    ),
    "verify_checkpoint": (
        "[TOOL: verify_checkpoint MATRIX_NAME v1 v2 v3 ...]\n"
        "   Verify an intermediate matrix against a compiled-Fortran "
        "checkpoint. Call with no arguments to list checkpoints."
    ),
}


def _tool_descriptions(enabled_tools: set[str]) -> str:
    lines = [
        "You have access only to the tools listed below. Call one by "
        "outputting a line starting with [TOOL: tool_name].",
        "",
    ]
    number = 1
    for name in TOOL_ORDER:
        if name in enabled_tools:
            lines.append(f"{number}. {TOOL_DETAILS[name]}")
            lines.append("")
            number += 1
    lines.extend([
        "When you believe your code is correct, output it with ```cpp "
        "fences and write [DONE] on a separate line.",
        "If you submit code without an explicit tool call, only the enabled "
        "testing tools will be used automatically.",
    ])
    return "\n".join(lines)


class ToolAugmentedAgent:
    """
    An agent that translates Fortran to C++ using tools iteratively.

    Unlike the fixed pipeline, the model decides what to do at each step.
    """

    def __init__(self, llm_client, config: dict, fortran_code: str,
                 logger=None, max_steps: int = 15, max_tokens: int = 100_000,
                 enabled_tools: Optional[set[str]] = None,
                 knowledge_mode: str = "full",
                 enable_checkpoints: bool = True):
        self.llm = llm_client
        self.config = config
        self.fortran_code = fortran_code
        self.logger = logger
        self.max_steps = max_steps
        self.max_tokens = max_tokens
        self.enabled_tools = set(TOOL_ORDER if enabled_tools is None
                                 else enabled_tools)
        unknown_tools = self.enabled_tools.difference(TOOL_ORDER)
        if unknown_tools:
            raise ValueError(f"Unknown agent tools: {sorted(unknown_tools)}")
        if "compile_and_run" not in self.enabled_tools:
            raise ValueError("compile_and_run is required for tool-agent runs")
        if knowledge_mode not in {"full", "off"}:
            raise ValueError("knowledge_mode must be 'full' or 'off'")
        self.knowledge_mode = knowledge_mode
        self.enable_checkpoints = enable_checkpoints

        # State
        self.conversation = []
        self.steps = 0
        self.tokens_used = 0
        self.last_stdout = ""
        self.last_code = ""
        self.last_runnable_code = ""
        self.best_code = ""
        self.best_rel_diff = float("inf")
        self.best_passed = False
        self.tool_calls = []

        # Source files available for read_source
        self._source_files = {}
        for src in config.get("source_files", []):
            p = PROJECT_ROOT / src
            if p.exists():
                self._source_files[p.name] = p.read_text(
                    encoding="utf-8", errors="replace")

        # Conventions (same as selfmine)
        from run_experiments_v3 import (
            _extract_entry_conventions, _format_entry_conventions,
        )
        self._conventions = _format_entry_conventions(
            _extract_entry_conventions(fortran_code))

        # Knowledge is explicitly controllable so confirmatory ablations do
        # not ingest insights harvested from the benchmark under evaluation.
        self._kb_context = (self._load_kb_context(config)
                            if knowledge_mode == "full" else "")

        # Reference and benchmark-specific tolerance for comparison.
        from run_experiments_v3 import (
            load_all_benchmarks, load_reference, verify_output,
        )
        self._reference = load_reference(config)
        self._verify = lambda computed: verify_output(computed, self._reference,
                                                     config)
        tolerance_map = load_all_benchmarks().get("tolerance", {})
        self._tolerance = tolerance_map.get(
            config.get("verification_type"), 1e-4)

        # Checkpoints are disabled in the core RQ3 tool ablation.
        checkpoints_allowed = (enable_checkpoints and
                               "verify_checkpoint" in self.enabled_tools)
        self._checkpoints = (self._load_checkpoints(config)
                             if checkpoints_allowed else {})

        # Build the initial prompt (same content as Direct)
        from run_experiments_v3 import build_generation_prompt
        self._initial_prompt = build_generation_prompt(config, fortran_code)

    def _log(self, msg):
        if self.logger:
            self.logger.log_info(f"  [agent] {msg}")

    @staticmethod
    def _load_kb_context(config: dict) -> str:
        """
        Load relevant knowledge base entries for this benchmark.

        Returns formatted context string with cross-file insights (GMMATD
        row-major, pivot assembly patterns, SMA1B transpose, etc.) that
        cannot be discovered from the source code alone.
        """
        try:
            from knowledge_base import KnowledgeBase
            kb = KnowledgeBase(str(PROJECT_ROOT / "knowledge_base.json"))
            insights = kb.query_for_benchmark(config)
            # Also explicitly query for GMMATD row-major (critical for plate bending)
            called = config.get("known_callees", [])
            if "GMMATD" in called:
                extra = kb.query(called_subroutines=["GMMATD"],
                                category=config.get("category", ""),
                                keywords=["gmmatd", "row", "storage", "transpose"])
                # Merge without duplicates
                seen_ids = {i.id for i in insights}
                for e in extra:
                    if e.id not in seen_ids:
                        insights.append(e)
                        seen_ids.add(e.id)
            if insights:
                return kb.format_context(insights, max_length=5000)
        except Exception:
            pass
        return ""

    @staticmethod
    def _load_checkpoints(config: dict) -> dict:
        """
        Load intermediate checkpoint matrices for incremental verification.
        Only loads if a checkpoint file specific to this benchmark exists.
        """
        ref_path = config.get("reference_output_path", "")
        if not ref_path:
            return {}
        # Only look for checkpoint files in the SAME directory as the reference,
        # and only if named after the benchmark's source routine
        ref_dir = (PROJECT_ROOT / ref_path).parent
        source_stems = [Path(s).stem.lower()
                       for s in config.get("source_files", [])]
        for f in ref_dir.iterdir() if ref_dir.exists() else []:
            if "checkpoint" not in f.name.lower():
                continue
            # Must match one of the source file stems
            if any(stem in f.name.lower() for stem in source_stems):
                try:
                    data = json.loads(f.read_text(encoding="utf-8"))
                    if isinstance(data, dict):
                        return data
                except (json.JSONDecodeError, ValueError):
                    pass
        return {}

    def _build_staged_system_prompt(self) -> str:
        """When checkpoints exist, force step-by-step implementation."""
        cp_names = list(self._checkpoints.keys())
        cp_list = "\n".join(f"  {i+1}. {name}" for i, name in enumerate(cp_names))

        base = (
            "You are a code migration agent. Your goal: translate the given "
            "NASTRAN-95 Fortran routine to C++17/Eigen such that the output "
            f"matches the compiled Fortran oracle within "
            f"{self._tolerance:.3e} relative error.\n\n"
            "CRITICAL STRATEGY: This routine is complex. You MUST work step "
            "by step. Do NOT attempt to write the complete program at once - "
            "that approach has been tried repeatedly and fails every time on "
            "this routine.\n\n"
            "Instead, follow this incremental approach:\n\n"
            "STEP 1: Write a COMPLETE C++ program that computes ONLY the "
            "first intermediate matrix below, and prints its values so I can "
            "verify it against the compiled Fortran's intermediate output.\n\n"
            "STEP 2: Once step 1 passes, EXTEND your program to compute the "
            "next intermediate. I will verify again.\n\n"
            "Continue until all intermediates are verified, then assemble "
            "the final output.\n\n"
            "Available intermediate checkpoints (in computation order):\n"
            f"{cp_list}\n\n"
            "IMPORTANT RULES:\n"
            "- Each submission must be a COMPLETE, compilable C++ file with "
            "main().\n"
            "- To verify an intermediate, call:\n"
            "  [TOOL: verify_checkpoint <name> val1 val2 ...]\n"
            "  Pass values as a flat row-major list of the matrix entries.\n"
            "  Call with no values to see expected size and first few values.\n"
            "- Do NOT skip ahead. Verify each step before adding the next.\n"
            "- After ALL intermediates pass, compute the full 9x9 stiffness "
            "matrix and print as JSON:\n"
            '  {"stiffness_matrix": [[row1], [row2], ...]}\n'
            "  Then call [TOOL: compile_and_run] and [TOOL: compare_with_reference].\n\n"
            "OTHER TOOLS:\n"
            "- [TOOL: read_source FILENAME START END] - read Fortran lines\n"
            "- [TOOL: inspect_conventions] - see extracted implementation conventions\n"
            "- [TOOL: compile_and_run] - compile and run your C++ code\n"
            "- [TOOL: compare_with_reference] - compare final output vs oracle\n\n"
            "START by calling [TOOL: verify_checkpoint] with no arguments to "
            "see the available checkpoints and their sizes. Then implement "
            "the first one.\n"
        )

        # Inject knowledge base context if available
        if self._kb_context:
            base += "\n" + self._kb_context + "\n"

        return base

    def _build_system_prompt(self) -> str:
        # If checkpoints exist, use the historical staged strategy. Core RQ3
        # ablations disable checkpoints, so all compared cells use the same
        # base prompt and differ only in the advertised tool set.
        if self._checkpoints:
            return self._build_staged_system_prompt()

        parts = [
            "You are a code migration agent. Your goal: translate the given "
            "NASTRAN-95 Fortran routine to C++17/Eigen such that the output "
            f"matches the compiled Fortran oracle within the benchmark's "
            f"relative-error tolerance ({self._tolerance:.3e}).",
            "",
            "Strategy tips:",
            "- Start by writing your best complete, self-contained C++ "
            "program with main().",
        ]
        if "compile_and_run" in self.enabled_tools:
            parts.append("- Compile submitted code and use compiler/runtime "
                         "feedback for targeted fixes.")
        if "compare_with_reference" in self.enabled_tools:
            parts.append("- Compare runnable output against the compiled "
                         "oracle before declaring success.")
        if "read_source" in self.enabled_tools:
            parts.append("- Re-read specific source lines when diagnostics "
                         "suggest a local implementation convention.")
        if "inspect_conventions" in self.enabled_tools:
            parts.append("- Inspect mechanically extracted conventions when "
                         "scaling, indexing, or decomposition is uncertain.")
        if "verify_checkpoint" in self.enabled_tools:
            parts.append("- Use available checkpoints to verify intermediate "
                         "matrices step by step.")
        parts.extend([
            "- Do not repeatedly guess; use only the available feedback.",
            "",
            _tool_descriptions(self.enabled_tools),
        ])

        if self._kb_context:
            parts.extend(["", self._kb_context])

        return "\n".join(parts)

    def _execute_tool(self, tool_line: str, code_before: str) -> str:
        """Parse and execute an enabled tool call, returning text feedback."""
        tool_line = tool_line.strip()
        requested = next((name for name in TOOL_ORDER if name in tool_line),
                         None)
        if requested is None:
            return f"Unknown tool: {tool_line}"
        if requested not in self.enabled_tools:
            return (f"UNAVAILABLE: {requested} is disabled for this run. "
                    "Use only the tools listed in the system prompt.")

        if requested == "compile_and_run":
            return self._tool_compile_and_run(code_before)
        if requested == "compare_with_reference":
            return self._tool_compare()
        if requested == "read_source":
            return self._tool_read_source(tool_line)
        if requested == "inspect_conventions":
            return self._tool_conventions()
        return self._tool_verify_checkpoint(tool_line)

    def _tool_compile_and_run(self, code: str) -> str:
        if not code or len(code) < 100:
            return "ERROR: No valid C++ code found before the tool call."

        self.last_code = code
        with tempfile.TemporaryDirectory(prefix="agent_") as td:
            result = compile_and_run_cpp(code, Path(td))

        if not result["success"]:
            stage = result.get("stage", "unknown")
            error = result.get("error", "")[:800]
            return f"FAILED at {stage}:\n{error}"

        self.last_stdout = result.get("stdout", "")
        self.last_runnable_code = code
        lines = self.last_stdout.count("\n")
        return (f"SUCCESS: compiled and ran. Output: {lines} lines, "
                f"{len(self.last_stdout)} chars.\n"
                f"First 200 chars of stdout: {self.last_stdout[:200]}")

    def _tool_compare(self) -> str:
        if not self.last_stdout:
            return "ERROR: No successful run output to compare. Run compile_and_run first."

        # Use the same parser as the harness, which handles all verification
        # types (matrix_symmetry, vector_compare, scalar_compare, etc.)
        from run_experiments_v3 import parse_output
        computed = parse_output(self.last_stdout, self.config)
        if computed is None:
            vtype = self.config.get("verification_type", "unknown")
            if vtype == "vector_compare":
                hint = ("Make sure you print the solution vector as JSON: "
                        '{"solution": [v1, v2, ...]} or as space-separated '
                        "numbers on one line.")
            elif vtype == "scalar_compare":
                hint = ("Make sure you print the scalar results as JSON: "
                        '{"axial_stress": ..., "axial_force": ..., '
                        '"torsional_stress": ...}')
            elif vtype == "eigenvalue_compare":
                hint = ("Make sure you print eigenvalues as JSON: "
                        '{"eigenvalues": [e1, e2, ...]}')
            else:
                hint = ("Make sure you print the output in a parseable format "
                        "(JSON with appropriate keys, or rows of numbers).")
            return f"ERROR: Could not parse the output. {hint}"

        vr = self._verify(computed)
        raw_rel = vr.get("max_rel_diff")
        rel = (float(raw_rel) if raw_rel is not None else float("inf"))

        # Track the best oracle-scored candidate only when oracle feedback is
        # enabled. Contract violations have no meaningful residual.
        if rel < self.best_rel_diff:
            self.best_rel_diff = rel
            self.best_code = self.last_code

        if vr.get("passed"):
            self.best_passed = True
            return (f"PASSED! max_rel_diff = {rel:.3e} "
                    f"(threshold {self._tolerance:.3e}).\n"
                    "Your code is numerically equivalent to the legacy kernel.")

        # Build diagnostic
        rel_text = f"{rel:.3e}" if np.isfinite(rel) else "not defined"
        parts = [f"FAILED. max_rel_diff = {rel_text} "
                 f"(threshold {self._tolerance:.3e})."]
        if vr.get("contract_violation"):
            parts.append(f"Contract violation: {vr.get('reason', 'shape mismatch')}")

        if vr.get("is_symmetric") is not None:
            parts.append(f"Symmetric: {vr['is_symmetric']}")

        # For vector/scalar types, show simple element-wise comparison
        if (self._reference is not None and computed is not None
                and (computed.ndim == 1 or
                     (computed.ndim == 2 and
                      computed.shape != self._reference.shape))):
            ref_flat = self._reference.flatten()
            got_flat = computed.flatten()
            n = min(len(ref_flat), len(got_flat))
            if n > 0:
                parts.append(f"\nYour output ({len(got_flat)} values) vs "
                             f"expected ({len(ref_flat)} values):")
                for idx in range(min(n, 8)):
                    g, r = got_flat[idx], ref_flat[idx]
                    ratio = g / r if abs(r) > 1e-30 else float("inf")
                    parts.append(f"  [{idx}]: yours={g:.6e} expected="
                                 f"{r:.6e} ratio={ratio:.4f}")
                if len(got_flat) != len(ref_flat):
                    parts.append(f"\n  LENGTH MISMATCH: you output "
                                 f"{len(got_flat)} values but {len(ref_flat)} "
                                 f"are expected.")

        # Per-direction analysis if we have a matrix (not vector/scalar)
        if (self._reference is not None and computed is not None
                and computed.ndim == 2
                and computed.shape == self._reference.shape):
            ref = self._reference
            dof = self.config.get("dof_per_node")
            mask = np.abs(ref) > 1e-9 * np.abs(ref).max()
            ratio = np.where(mask, computed / np.where(mask, ref, 1), np.nan)

            # Worst cells
            diff = np.abs(computed - ref)
            flat_idx = np.argsort(diff.flatten())[::-1][:5]
            parts.append("\nWorst 5 cells (row, col): computed vs expected:")
            for idx in flat_idx:
                r, c = divmod(int(idx), computed.shape[1])
                if diff[r, c] > 0:
                    parts.append(f"  K[{r},{c}]: {computed[r,c]:.6e} vs "
                                 f"{ref[r,c]:.6e}")

            # Direction-grouped ratios
            if isinstance(dof, int) and dof in (2, 3, 6):
                names = {2: ["x", "y"], 3: ["x", "y", "z"],
                         6: ["x", "y", "z", "rx", "ry", "rz"]}[dof]
                parts.append("\nMean ratio (computed/reference) by direction:")
                for a in range(dof):
                    sub = ratio[a::dof, a::dof]
                    fin = sub[np.isfinite(sub)]
                    if fin.size:
                        parts.append(f"  {names[a]}-{names[a]}: "
                                     f"{float(fin.mean()):.6f}")

            # Overall scale
            fin_all = ratio[np.isfinite(ratio)]
            if fin_all.size:
                mean_r = float(fin_all.mean())
                # A near-zero mean ratio means the candidate is essentially the
                # zero matrix (or otherwise degenerate), so the reciprocal
                # "missing factor" hint is undefined. Reporting the raw ratio
                # without dividing avoids crashing the whole campaign on a
                # single degenerate candidate.
                if abs(mean_r - 1.0) > 0.01 and abs(mean_r) > 1e-9:
                    parts.append(f"\nOverall mean ratio: {mean_r:.4f} "
                                 f"(you may be missing a factor of "
                                 f"{1.0/mean_r:.4f})")
                elif abs(mean_r) <= 1e-9:
                    parts.append(
                        f"\nOverall mean ratio: {mean_r:.4f} "
                        f"(the output is essentially zero; check that the "
                        f"matrix is actually being populated)")
                elif abs(mean_r - 1.0) < 0.005 and rel < 0.01:
                    parts.append(
                        "\nNEAR MISS: structure is correct, a small "
                        "coefficient is off. Do NOT rewrite the algorithm.")

            # Sub-block transpose detection.
            # Many FEA stiffness matrices are composed of dof×dof sub-blocks,
            # one per node pair. A common assembly error transposes some of
            # these blocks, which flips the off-diagonal signs within the
            # block while the full matrix can still appear symmetric (because
            # K[i,j] and K[j,i] are both wrong). The whole-matrix transpose
            # check misses this. Automatically derives the block size from
            # dof_per_node so it works on any benchmark.
            #
            # Strategy: try transposing ALL off-diagonal blocks at once. On
            # B7 the error is not in one block but in all of them (assembly
            # index ordering), so the single-block test never triggered.
            if (isinstance(dof, int) and dof >= 2
                    and computed.shape[0] % dof == 0
                    and rel > 1e-4):
                n_nodes = computed.shape[0] // dof
                max_ref = float(np.abs(ref).max())

                # Test 1: transpose ALL off-diagonal blocks
                K_fixed = computed.copy()
                for bi in range(n_nodes):
                    for bj in range(n_nodes):
                        if bi != bj:
                            rs = slice(bi * dof, (bi + 1) * dof)
                            cs = slice(bj * dof, (bj + 1) * dof)
                            K_fixed[rs, cs] = computed[rs, cs].T
                rel_fixed_all = float(np.abs(K_fixed - ref).max() / max_ref)

                if rel_fixed_all < rel * 0.01:
                    # Massive improvement: the assembly transposes off-diag
                    parts.append(
                        f"\nPATTERN: BLOCK TRANSPOSE detected. Transposing "
                        f"ALL off-diagonal {dof}x{dof} sub-blocks reduces "
                        f"the error from {rel:.2e} to {rel_fixed_all:.2e}.")
                    parts.append(
                        f"FIX: In your assembly loop, when you write the "
                        f"sub-block for nodes i and j (i != j), the row and "
                        f"column DOF indices WITHIN the block are swapped. "
                        f"The most common cause: you compute Ke as a local "
                        f"element matrix and scatter it into the global K, "
                        f"but for the off-diagonal blocks you should use "
                        f"Ke_sub.transpose() instead of Ke_sub. "
                        f"Alternatively, check whether your B-matrix columns "
                        f"list the DOF in a different order than your "
                        f"assembly indexing assumes.")
                else:
                    # Test 2: try transposing individual blocks
                    best_fix = None
                    for bi in range(n_nodes):
                        for bj in range(n_nodes):
                            rs = slice(bi * dof, (bi + 1) * dof)
                            cs = slice(bj * dof, (bj + 1) * dof)
                            blk_c = computed[rs, cs]
                            blk_r = ref[rs, cs]
                            err_now = float(np.abs(blk_c - blk_r).max())
                            err_t = float(np.abs(blk_c.T - blk_r).max())
                            if (err_now > max_ref * 1e-3
                                    and err_t < err_now * 0.1):
                                imp = err_now / max(err_t, 1e-30)
                                if best_fix is None or imp > best_fix[3]:
                                    best_fix = (bi, bj, err_t / max_ref, imp)
                    if best_fix:
                        bi, bj, new_rel, imp = best_fix
                        parts.append(
                            f"\nPATTERN: BLOCK TRANSPOSE detected. The "
                            f"{dof}x{dof} sub-block at node pair ({bi},{bj}) "
                            f"appears transposed: transposing it reduces the "
                            f"block error by {imp:.0f}x.")
                        parts.append(
                            f"FIX: Check the assembly index order for node "
                            f"pair ({bi},{bj}). The row/column DOF within "
                            f"that block are swapped.")

        return "\n".join(parts)

    def _tool_read_source(self, tool_line: str) -> str:
        # Parse: [TOOL: read_source filename start end]
        # Models may format this in various ways, e.g.:
        #   [TOOL: read_source kqdmem.f 110 120]
        #   [TOOL: read_source file=kqdmem.f start=110 end=120]
        #   [TOOL: read_source kqdmem.f kqdmem.f 110 120]  (repeated name)
        cleaned = (tool_line.replace("[TOOL:", "").replace("]", "")
                   .replace("read_source", "").strip())
        # Remove keyword= prefixes
        cleaned = re.sub(r"\b(file|filename|start|end|start_line|end_line)\s*=\s*",
                         "", cleaned)
        parts = cleaned.split()

        # Separate filename tokens (non-numeric) from line numbers (numeric)
        names = []
        nums = []
        for p in parts:
            try:
                nums.append(int(p))
            except ValueError:
                names.append(p)

        fname = names[0] if names else ""
        start = nums[0] if len(nums) > 0 else 1
        end = nums[1] if len(nums) > 1 else start + 30

        # Find the file
        matched = None
        for name, content in self._source_files.items():
            if fname.lower() in name.lower():
                matched = (name, content)
                break

        if not matched:
            return (f"File '{fname}' not found. "
                    f"Available: {list(self._source_files.keys())}")

        name, content = matched
        lines = content.splitlines()
        start = max(1, min(start, len(lines)))
        end = min(end, len(lines))

        result_lines = []
        for i in range(start - 1, end):
            result_lines.append(f"{i+1:5d}| {lines[i]}")

        return f"--- {name} lines {start}-{end} ---\n" + "\n".join(result_lines)

    def _tool_conventions(self) -> str:
        if not self._conventions:
            return "No conventions extracted for this routine."
        return self._conventions

    def _tool_verify_checkpoint(self, tool_line: str) -> str:
        """
        Verify an intermediate matrix against the compiled Fortran checkpoint.

        This is the automation of the strategy that solved B5 by hand:
        translate step by step, verify each intermediate against the compiled
        Fortran's output, only assemble once all parts are confirmed.
        """
        if not self._checkpoints:
            return ("No checkpoints available for this benchmark. "
                    "Use compile_and_run + compare_with_reference instead.")

        # Parse arguments
        cleaned = (tool_line.replace("[TOOL:", "").replace("]", "")
                   .replace("verify_checkpoint", "").strip())
        parts = cleaned.split()

        # If no arguments, list available checkpoints
        if not parts:
            lines = ["Available checkpoints (intermediate matrices from the "
                     "compiled Fortran):"]
            for name, vals in self._checkpoints.items():
                arr = np.array(vals) if isinstance(vals, list) else None
                if arr is not None:
                    n = len(arr)
                    side = int(np.sqrt(n))
                    lines.append(f"  {name}: {side}x{side} ({n} values), "
                                 f"|max|={np.abs(arr).max():.3e}")
                else:
                    lines.append(f"  {name}: {vals}")
            lines.append("\nUsage: [TOOL: verify_checkpoint MATRIX_NAME "
                         "val1 val2 val3 ...]")
            lines.append("Pass values as a flat row-major list.")
            lines.append("\nStrategy: compute each intermediate matrix in "
                         "your code and verify it here before proceeding to "
                         "the next. This catches errors early rather than "
                         "discovering at the end that the final matrix is "
                         "wrong by a factor of 3000.")
            return "\n".join(lines)

        # Find the checkpoint name (first non-numeric token)
        name_parts = []
        val_parts = []
        for p in parts:
            try:
                float(p)
                val_parts.append(p)
            except ValueError:
                if not val_parts:
                    name_parts.append(p)
                else:
                    val_parts.append(p)  # might be e.g. "1.2e-3" split weird

        cp_name = " ".join(name_parts) if name_parts else ""

        # Match checkpoint name (fuzzy)
        matched_key = None
        for key in self._checkpoints:
            if cp_name.lower() in key.lower() or key.lower() in cp_name.lower():
                matched_key = key
                break
        if matched_key is None:
            # Try partial match
            for key in self._checkpoints:
                if any(w.lower() in key.lower() for w in name_parts):
                    matched_key = key
                    break

        if matched_key is None:
            return (f"Checkpoint '{cp_name}' not found. "
                    f"Available: {list(self._checkpoints.keys())}")

        ref_vals = np.array(self._checkpoints[matched_key], dtype=float)

        if not val_parts:
            # Just show what this checkpoint expects
            n = len(ref_vals)
            side = int(np.sqrt(n))
            return (f"Checkpoint '{matched_key}': expects {n} values "
                    f"({side}x{side} matrix, row-major).\n"
                    f"|max| = {np.abs(ref_vals).max():.6e}\n"
                    f"First 6 values: {ref_vals[:6].tolist()}\n"
                    f"Call with your computed values to check.")

        # Parse the submitted values
        try:
            submitted = np.array([float(v) for v in val_parts], dtype=float)
        except ValueError:
            return "ERROR: could not parse the values as numbers."

        if len(submitted) != len(ref_vals):
            return (f"ERROR: expected {len(ref_vals)} values for "
                    f"'{matched_key}', got {len(submitted)}.")

        # Compare
        max_ref = max(np.abs(ref_vals).max(), 1e-30)
        abs_diff = np.abs(submitted - ref_vals)
        rel_diff = float(abs_diff.max() / max_ref)

        if rel_diff < 1e-4:
            return (f"CHECKPOINT PASSED: '{matched_key}' matches to "
                    f"rel={rel_diff:.2e}. This intermediate is correct, "
                    f"proceed to the next step.")

        # Detailed diagnosis
        n = len(ref_vals)
        side = int(np.sqrt(n))
        worst_idx = int(np.argmax(abs_diff))
        worst_row = worst_idx // side if side > 0 else 0
        worst_col = worst_idx % side if side > 0 else 0

        parts_out = [
            f"CHECKPOINT FAILED: '{matched_key}' max_rel_diff = {rel_diff:.3e}",
            f"  Worst entry: [{worst_row},{worst_col}] "
            f"yours={submitted[worst_idx]:.6e} "
            f"expected={ref_vals[worst_idx]:.6e}",
        ]

        # Check for common patterns
        nonzero = np.abs(ref_vals) > 1e-9 * max_ref
        if nonzero.sum() > 3:
            ratio = submitted[nonzero] / ref_vals[nonzero]
            mean_ratio = float(np.mean(ratio))
            if abs(mean_ratio) > 0.01 and np.std(ratio) / abs(mean_ratio) < 0.1:
                parts_out.append(
                    f"  PATTERN: uniform scale factor {mean_ratio:.4f} "
                    f"(expected 1.0). You may be missing a multiply by "
                    f"{1.0/mean_ratio:.4f}")

            # Check transpose
            if side > 0:
                sub_t = submitted.copy()
                mat = sub_t.reshape(side, side)
                mat_t = mat.T.flatten()
                rel_t = float(np.abs(mat_t - ref_vals).max() / max_ref)
                if rel_t < rel_diff * 0.1:
                    parts_out.append(
                        f"  PATTERN: TRANSPOSED. Your matrix appears to be "
                        f"the transpose of what is expected (rel after "
                        f"transpose: {rel_t:.3e}). This usually means "
                        f"row-major vs column-major confusion in GMMATD calls.")

        parts_out.append(f"\n  First 6 expected: {ref_vals[:6].tolist()}")
        parts_out.append(f"  First 6 yours:    {submitted[:6].tolist()}")

        # Row-by-row comparison for the first few rows
        if side > 0 and side <= 6:
            parts_out.append(f"\n  Row-by-row (yours | expected | ratio):")
            got_mat = submitted.reshape(side, side)
            ref_mat = ref_vals.reshape(side, side)
            for row in range(min(3, side)):
                for col in range(side):
                    g = got_mat[row, col]
                    r = ref_mat[row, col]
                    rat = g / r if abs(r) > 1e-20 else float("inf")
                    if abs(g - r) > 1e-9 * max_ref:
                        parts_out.append(
                            f"    [{row},{col}]: {g:12.4e} | "
                            f"{r:12.4e} | ratio {rat:.4f}")

        return "\n".join(parts_out)

    def _extract_code(self, response: str) -> str:
        """Extract C++ code from the response."""
        if "```cpp" in response:
            parts = response.split("```cpp")
            if len(parts) > 1:
                code = parts[-1].split("```")[0].strip()
                if code and len(code) > 100:
                    return code
        if "```c++" in response:
            parts = response.split("```c++")
            if len(parts) > 1:
                code = parts[-1].split("```")[0].strip()
                if code and len(code) > 100:
                    return code
        return ""

    def _find_tool_call(self, response: str) -> Optional[str]:
        """Find a tool call in the response."""
        for line in response.splitlines():
            if line.strip().startswith("[TOOL:"):
                return line.strip()
            if line.strip() == "[DONE]":
                return "[DONE]"
        return None

    def run(self) -> dict:
        """
        Execute the agent loop.

        Returns dict with: passed, code, rel_diff, steps, tokens, tool_calls
        """
        system = self._build_system_prompt()
        self.conversation = [{"role": "user", "content": self._initial_prompt}]

        self._log(f"starting, max {self.max_steps} steps, "
                  f"source files: {list(self._source_files.keys())}, "
                  f"tools: {sorted(self.enabled_tools)}, "
                  f"knowledge={self.knowledge_mode}, "
                  f"checkpoints={'on' if self._checkpoints else 'off'}")

        while self.steps < self.max_steps and self.tokens_used < self.max_tokens:
            # Build the conversation for the LLM
            conv_text = "\n\n".join(
                f"{'USER' if m['role']=='user' else 'ASSISTANT' if m['role']=='assistant' else 'TOOL RESULT'}: "
                f"{m['content']}"
                for m in self.conversation[-12:]  # keep context bounded
            )

            response = self.llm.generate(
                system_prompt=system,
                user_prompt=conv_text,
                agent_name=f"agent_step{self.steps}")

            if not response:
                self._log("empty response, retrying")
                self.steps += 1
                continue

            # Estimate tokens (rough)
            self.tokens_used += len(conv_text) // 3 + len(response) // 3

            self.conversation.append({"role": "assistant", "content": response})

            # Extract code if present
            code = self._extract_code(response)

            # Check for tool call
            tool_call = self._find_tool_call(response)

            if tool_call == "[DONE]":
                self._log(f"agent declared done at step {self.steps}")
                if code:
                    self.last_code = code
                    if self.best_rel_diff == float("inf"):
                        self.best_code = code
                break

            elif tool_call:
                self._log(f"step {self.steps}: {tool_call[:60]}")
                self.tool_calls.append({"step": self.steps, "tool": tool_call})

                # Execute the tool
                tool_result = self._execute_tool(tool_call, code)
                self.conversation.append(
                    {"role": "user",
                     "content": f"[Tool Result]:\n{tool_result}"})

                # Check if we passed
                if "PASSED!" in tool_result:
                    self._log(f"PASSED at step {self.steps}")
                    break

            elif code:
                # Code without an explicit tool call: auto-test only with the
                # tools enabled for this ablation cell. In particular, a
                # no-oracle run must never receive numerical diagnostics here.
                self._log(f"step {self.steps}: code submitted, auto-testing")
                auto_label = ("auto_compile_and_compare"
                              if "compare_with_reference" in self.enabled_tools
                              else "auto_compile_only")
                self.tool_calls.append(
                    {"step": self.steps, "tool": auto_label})

                compile_result = self._tool_compile_and_run(code)
                if ("SUCCESS" in compile_result and
                        "compare_with_reference" in self.enabled_tools):
                    compare_result = self._tool_compare()
                    self.conversation.append(
                        {"role": "user",
                         "content": f"[Auto-test Result]:\n{compile_result}\n\n"
                                    f"[Comparison with oracle]:\n{compare_result}"})
                    if "PASSED!" in compare_result:
                        self._log(f"PASSED at step {self.steps}")
                        break
                else:
                    self.conversation.append(
                        {"role": "user",
                         "content": f"[Auto-test Result]:\n{compile_result}"})
            else:
                # No code, no tool — the model is just thinking/explaining
                self._log(f"step {self.steps}: reasoning (no action)")
                # Nudge it to do something
                self.conversation.append(
                    {"role": "user",
                     "content": "Please either write code (in ```cpp fences) "
                                "or call a tool. Do not just explain."})

            self.steps += 1

        # Final result. Oracle-disabled runs return the last runnable program;
        # the external harness scores it after the conversation without
        # leaking that score back to the model.
        passed = self.best_passed
        if "compare_with_reference" in self.enabled_tools:
            final_code = (self.best_code or self.last_runnable_code or
                          self.last_code)
        else:
            final_code = self.last_runnable_code or self.last_code
        return {
            "passed": passed,
            "code": final_code,
            "rel_diff": (self.best_rel_diff
                         if self.best_rel_diff < float("inf") else None),
            "steps": self.steps,
            "tokens_used": self.tokens_used,
            "tool_calls": self.tool_calls,
            "enabled_tools": sorted(self.enabled_tools),
            "knowledge_mode": self.knowledge_mode,
            "checkpoints_enabled": bool(self._checkpoints),
            "conversation_length": len(self.conversation),
        }


# ============================================================
# Integration point for run_experiments_v3.py
# ============================================================

def run_tool_augmented(model: str, config: dict, fortran_code: str,
                      logger, llm_client,
                      max_steps: int = 15, max_tokens: int = 100_000,
                      enabled_tools: Optional[set[str]] = None,
                      knowledge_mode: str = "full",
                      enable_checkpoints: bool = True) -> str:
    """Run the tool-augmented agent and return its selected final code."""
    agent = ToolAugmentedAgent(
        llm_client=llm_client,
        config=config,
        fortran_code=fortran_code,
        logger=logger,
        max_steps=max_steps,
        max_tokens=max_tokens,
        enabled_tools=enabled_tools,
        knowledge_mode=knowledge_mode,
        enable_checkpoints=enable_checkpoints,
    )
    result = agent.run()

    if logger:
        logger.log_info(f"  [agent] finished: passed={result['passed']} "
                        f"steps={result['steps']} "
                        f"rel_diff={result.get('rel_diff')} "
                        f"tools={len(result['tool_calls'])} "
                        f"enabled={result['enabled_tools']}")
        logger.log_artifact("agent_result", {
            "passed": result["passed"],
            "steps": result["steps"],
            "tokens_used": result["tokens_used"],
            "rel_diff": result.get("rel_diff"),
            "tool_calls": result["tool_calls"],
            "enabled_tools": result["enabled_tools"],
            "knowledge_mode": result["knowledge_mode"],
            "checkpoints_enabled": result["checkpoints_enabled"],
            "conversation_length": result["conversation_length"],
        })

    return result.get("code", "")
