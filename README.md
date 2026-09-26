# LLM-Assisted Legacy Fortran (NASTRAN) Migration — Code and Artifacts

This repository contains the code, frozen inputs, oracles, and raw run
artifacts for an empirical study that evaluates LLM-based migration of
NASTRAN-95 Fortran kernels to modern C++/Eigen. Three workflow packages are
compared — **Direct**, **Staged**, and **Full Agent** — and the single-input
accepted candidates are re-examined on a frozen set of held-out multi-input
variants (h1–h3).

## 1. Benchmarks

13 NASTRAN-95 kernels (plus the B5a KTRBSC configuration) are defined in
[`benchmarks_config.json`](benchmarks_config.json):

| ID | Kernel | Output contract |
|---|---|---|
| B1 | CROD | symmetric stiffness matrix |
| B2 | CBAR | symmetric stiffness matrix |
| B3 | CTRIA3 | symmetric stiffness matrix |
| B4 | CQUAD4 | symmetric stiffness matrix |
| B5/B5a | KTRPLT / KTRBSC | symmetric stiffness matrix |
| B6 | KQDPLT | symmetric stiffness matrix |
| B7 | KTETRA | symmetric stiffness matrix (12×12) |
| B8 | KSOLID | symmetric stiffness matrix |
| B9 | KELBOW | symmetric stiffness matrix |
| B10 | INVERD | solution vector (n = 4) |
| B11 | QRITER | eigenvalues (n = 5) |
| B12 | TRD1C | time series |
| B13 | SROD1 | scalar |

Acceptance tolerances (per verification type) are in
[`benchmarks_config.json`](benchmarks_config.json): matrix/vector/scalar
`1e-4`, eigenvalues `1e-6`, time series `1e-3`. A candidate passes when the
relative \(L_\infty\) error is below tolerance or the absolute error is below
`1e-10`. Matrix shape mismatches (12×12 vs 6×6) are resolved by DOF
submatrix extraction before comparison.

## 2. Repository layout

```
nastran_refactor/
├── run_experiments_v3.py         # Oracle verifier + main experiment runner
├── tool_augmented_agent.py       # Full Agent controller
├── run_jss_rq3_ablation.py       # RQ3 / Qwen Direct orchestration
├── run_multilang_v2.py           # Multi-language execution
├── config.py                     # Model endpoints (API keys redacted)
├── llm_client.py                 # OpenAI-compatible LLM client
├── experiment_logger.py          # Run logging utilities
├── knowledge_base.py / .json     # Migration knowledge base
├── benchmarks_config.json        # Benchmarks, tolerances, contracts
├── agents/                       # Full Agent multi-agent modules
├── nastran/NASTRAN-95/           # Legacy NASTRAN-95 Fortran sources (mis/, um/, ...)
├── source_fortran/test_drivers/  # Frozen Fortran drivers + reference outputs
├── holdout_manifest.json         # 13 × 4 held-out inputs with SHA-256
├── holdout_oracle_manifest.json  # SHA-256 of 52 held-out oracle records
├── holdout_inputs/               # 52 held-out input records (B{id}_h{k}.json)
├── holdout_oracles/              # 52 oracle records (kernel hash, driver/stub)
├── analysis/                     # Aggregated statistics, provenance, hold-out pipeline
│   ├── jss_v4_analysis.json          # Main-experiment statistics (520 observations)
│   ├── jss_v4_provenance_ledger.json # Frozen code/candidate/oracle hashes
│   ├── jss_v4_canonical_rescore.json # P0.5 zero-flip re-score
│   ├── jss_v4_holdout_baseline.json  # 173 canonical candidates
│   ├── jss_holdout_adapter_scan_v2.json
│   ├── holdout_results.json          # P3b per-variant outcomes
│   ├── aggregate_holdout_p4.json     # P4 retention metrics + bootstrap CIs
│   ├── jss_holdout_p3b_p4_provenance.json
│   ├── build_holdout_*.py            # Input/oracle/baseline builders
│   ├── holdout_adapter_v2.py         # Adapter scan
│   ├── run_holdout_p3b.py            # Offline revalidation
│   └── aggregate_holdout_p4.py       # Retention aggregation
└── results/                      # Raw run artifacts (frozen batches only)
    ├── paper1_gap_to_n10_all_models.json
    ├── jss_full_agent_main/jss_full_agent_main_n10_frozen_v2/
    └── jss_qwen37plus_main/jss_qwen37plus_direct_n10_frozen_v1/
```

The two `*_frozen_*` result directories contain per-run prompts, generated
candidates, compilation/verification outputs, and `final_summary.json` /
`rq3_protocol.json` summaries.

## 3. Environment

- **Python 3.9+** with: `numpy`, `scipy`, `openai`
  (the OpenAI SDK is used against OpenAI-compatible endpoints).
- **gfortran** with flags `-ffixed-form -std=legacy -w
  -fallow-argument-mismatch` for building legacy reference drivers.
- **MSVC (`cl.exe`)** and **Eigen 3.4.0** for compiling migrated C++
  candidates.
- API credentials: edit `config.py` and replace the four placeholders with
  your own keys (OpenRouter, API2D/OpenAI-compatible, DeepSeek, Qwen
  DashScope). **Never commit real keys.**

## 4. Running the experiments

`run_experiments_v3.py` generates candidates through the configured workflow,
compiles them against the frozen drivers in `source_fortran/`, and applies
the oracle defined in `benchmarks_config.json`. The Full Agent workflow is
driven by `tool_augmented_agent.py` together with the modules under
`agents/`; `run_jss_rq3_ablation.py` orchestrates the additional Qwen Direct
cohort.

The frozen held-out inputs and oracles used for multi-input revalidation are
provided in `holdout_inputs/`, `holdout_oracles/`, `holdout_manifest.json`,
and `holdout_oracle_manifest.json` (every record carries SHA-256 hashes and
the reference driver/stub sources). The analysis pipeline scripts and their
aggregated outputs live under `analysis/`:

1. **Build inputs and oracles** — `analysis/build_holdout_inputs.py`,
   `analysis/build_holdout_oracles.py`.
2. **Build the canonical baseline** — `analysis/build_holdout_baseline.py`
   (173 candidates; output `jss_v4_holdout_baseline.json`).
3. **Adapter scan** — `analysis/holdout_adapter_v2.py` classifies candidates
   as `clean`, `geometry`, or `not_parameterizable`
   (`jss_holdout_adapter_scan_v2.json`).
4. **P3b offline revalidation** — `analysis/run_holdout_p3b.py` re-tests
   parameterizable candidates on h1/h2/h3 (`holdout_results.json`).
5. **P4 aggregation** — `analysis/aggregate_holdout_p4.py` computes
   conservative, conditional, and upper-bound retention with benchmark-level
   bootstrap CIs (`aggregate_holdout_p4.json`).

A candidate is **kept** only if it passes all three variants h1–h3. The
scripts resolve the repository root from their own location, so they must be
run from inside `analysis/` (e.g. `python run_holdout_p3b.py`).

## 5. Provenance and integrity

- Every held-out input and oracle record is pinned with SHA-256 in
  `holdout_manifest.json` and `holdout_oracle_manifest.json`; oracle
  records additionally include per-kernel source byte hashes.
- The main-experiment provenance ledger and the P3b/P4 provenance record are
  `analysis/jss_v4_provenance_ledger.json` and
  `analysis/jss_holdout_p3b_p4_provenance.json`. The P0.5 re-score
  (`analysis/jss_v4_canonical_rescore.json`) records zero outcome flips
  (`true->true = 173`, `false->false = 347`).
- Legacy NASTRAN-95 kernel sources are included unmodified under
  `nastran/NASTRAN-95/`.
- Compiled binaries, NASTRAN PDF manuals, run outputs, and intermediate
  batches are excluded; binaries are rebuilt from the included sources.

## 6. Notes

- Model generations use temperature `0.0`; provider-side nondeterminism may
  still occur, so token counts and generations can vary across reruns.
- `config.py` in this repository contains placeholder keys only.
