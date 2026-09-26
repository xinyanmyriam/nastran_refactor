"""
Refactoring Workflow Orchestrator.

Coordinates the full refactoring pipeline:
Source Code → Comprehension → Source Mining → Architecture → Specification → Design → Generation → Verification

Implements the Refactoring Workflow Graph (RWG) with quality gates at each stage.
All operations are logged with timestamps via ExperimentLogger.

v2 additions:
- Source Annotation Mining stage (between Comprehension and Architecture)
- Knowledge Base integration (optional, for cross-benchmark learning)
"""

import json
import time
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .comprehension_agent import ComprehensionAgent, CodeIntent
from .architecture_agent import ArchitectureAgent, ModuleArchitecture
from .physics_spec_agent import PhysicsSpecAgent, PhysicsSpecification
from .modern_design_agent import ModernDesignAgent, ModernDesign
from .code_gen_agent import CodeGenAgent, GenerationResult
from .equivalence_agent import EquivalenceAgent, EquivalenceReport
from .source_mining_agent import SourceMiningAgent, MiningResult


@dataclass
class RefactoringState:
    """Complete state of the refactoring workflow."""
    
    # Input
    source_files: dict[str, str] = field(default_factory=dict)
    
    # Stage outputs
    code_intents: list[dict] = field(default_factory=list)
    source_mining: dict = field(default_factory=dict)  # v2: mining results
    architecture: dict = field(default_factory=dict)
    physics_spec: dict = field(default_factory=dict)
    modern_design: dict = field(default_factory=dict)
    generated_code: str = ""
    equivalence_report: dict = field(default_factory=dict)
    
    # v2: Additional context from knowledge base
    kb_context: str = ""
    
    # Metadata
    stage: str = "init"
    quality_gates: dict[str, bool] = field(default_factory=dict)
    timestamps: dict[str, float] = field(default_factory=dict)
    iterations: dict[str, int] = field(default_factory=dict)
    errors: list[dict] = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "source_files": list(self.source_files.keys()),
            "stage": self.stage,
            "quality_gates": self.quality_gates,
            "timestamps": self.timestamps,
            "iterations": self.iterations,
            "code_intents": self.code_intents,
            "source_mining": self.source_mining,
            "architecture": self.architecture,
            "physics_spec": self.physics_spec,
            "modern_design": self.modern_design,
            "generated_code_length": len(self.generated_code),
            "equivalence_report": self.equivalence_report,
            "kb_context": self.kb_context[:200] if self.kb_context else "",
            "errors": self.errors,
        }


class RefactoringOrchestrator:
    """Orchestrates the full NASTRAN refactoring pipeline with comprehensive logging."""
    
    def __init__(self, llm_client, logger, output_dir: str = "results",
                 enable_source_mining: bool = False,
                 knowledge_base=None,
                 nastran_source_dir: str = None,
                 skip_unused_stages: bool = False):
        """
        Args:
            llm_client: LLM client with generate(system_prompt, user_prompt) method.
            logger: ExperimentLogger instance for recording all operations.
            output_dir: Directory to save intermediate results.
            enable_source_mining: If True, run source annotation mining stage (v2).
            knowledge_base: Optional KnowledgeBase instance for cross-benchmark learning.
            nastran_source_dir: Path to NASTRAN-95 mis/ directory for source mining.
            skip_unused_stages: If True, omit the code_generation and equivalence
                stages. When this orchestrator is driven by
                run_experiments_v3.generate_code, neither stage contributes to
                the artifact that is finally compiled:

                  * code_generation asks CodeGenAgent for PYTHON
                    (agents/code_gen_agent.py: "produces Python
                    implementations", parses ```python fences) while the
                    benchmark compiles C++17 with MSVC and Eigen. The caller
                    discards state.generated_code and issues its own C++ call.
                    Measured on the archive, this stage is 55.8% of the token
                    spend of a `multiagent` repeat.
                  * equivalence is invoked with reference_results=[], so it
                    checks that discarded Python against zero test cases and
                    always reports not-equivalent.

                Defaults to False so the historical `multiagent` and
                `multiagent_v2` behaviour stays bit-for-bit reproducible; the
                existing 22/35 baseline belongs to that configuration.
        """
        self.llm_client = llm_client
        self.logger = logger
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # v2 configuration
        self.enable_source_mining = enable_source_mining
        self.knowledge_base = knowledge_base
        self.skip_unused_stages = skip_unused_stages
        
        # Initialize agents
        self.comprehension_agent = ComprehensionAgent(llm_client)
        self.architecture_agent = ArchitectureAgent(llm_client)
        self.physics_spec_agent = PhysicsSpecAgent(llm_client)
        self.modern_design_agent = ModernDesignAgent(llm_client)
        self.code_gen_agent = CodeGenAgent(llm_client)
        self.equivalence_agent = EquivalenceAgent(llm_client)
        
        # v2: Source mining agent (only initialized if enabled)
        self.source_mining_agent = None
        if enable_source_mining:
            self.source_mining_agent = SourceMiningAgent(
                llm_client, nastran_source_dir=nastran_source_dir
            )
        
        # State
        self.state = RefactoringState()
        
        self.logger.log_info(f"Orchestrator initialized with output_dir={output_dir}")
        if enable_source_mining:
            self.logger.log_info("  Source Mining: ENABLED")
        if knowledge_base:
            self.logger.log_info(f"  Knowledge Base: ENABLED ({knowledge_base.get_stats().get('total', 0)} insights)")
    
    def run(
        self,
        source_files: dict[str, str],
        context: str = "",
        reference_results: list[dict] = None,
        benchmark_config: dict = None,
    ) -> RefactoringState:
        """
        Execute the full refactoring pipeline.
        
        Args:
            source_files: Dict mapping filename -> Fortran source code.
            context: Additional context (documentation, etc.).
            reference_results: Known correct results for equivalence testing.
            benchmark_config: Optional benchmark config dict (for KB queries).
            
        Returns:
            Final RefactoringState with all artifacts.
        """
        self.state = RefactoringState(source_files=source_files)
        
        # Store benchmark config for use by source mining
        self._benchmark_config = benchmark_config
        
        self.logger.log_info(f"Pipeline started with {len(source_files)} source file(s)")
        for fname in source_files:
            self.logger.log_info(f"  Source: {fname} ({len(source_files[fname])} chars)")
        
        # v2: Query knowledge base for relevant context
        if self.knowledge_base and benchmark_config:
            self._inject_kb_context(benchmark_config)
        
        # Build stage list (v2: conditionally include source mining)
        stages = [
            ("comprehension", lambda: self._run_comprehension(context)),
        ]
        
        if self.enable_source_mining:
            stages.append(("source_mining", lambda: self._run_source_mining(context)))
        
        stages.extend([
            ("architecture", lambda: self._run_architecture_recovery(context)),
            ("physics_spec", lambda: self._run_physics_specification(context)),
            ("modern_design", lambda: self._run_modern_design()),
        ])

        # code_generation emits Python and equivalence checks it against an
        # empty reference set, so neither affects what gets compiled when the
        # caller supplies its own C++ generation step. See __init__ docstring.
        if not self.skip_unused_stages:
            stages.extend([
                ("code_generation", lambda: self._run_code_generation()),
                ("equivalence",
                 lambda: self._run_equivalence_verification(reference_results)),
            ])
        else:
            self.logger.log_info(
                "  Skipping code_generation and equivalence: the caller "
                "performs its own C++ generation, so both stages are dead work")
        
        for stage_name, stage_fn in stages:
            self.logger.log_stage_start(stage_name)
            
            try:
                stage_fn()
                success = True
            except Exception as e:
                success = False
                self.logger.log_error(stage_name, e, context=f"Stage {stage_name} raised exception")
                self.state.errors.append({
                    "stage": stage_name,
                    "error": str(e),
                    "traceback": traceback.format_exc(),
                })
            
            self.logger.log_stage_end(stage_name, success)
            
            # Quality gate check
            gate_passed = self._check_quality_gate(stage_name)
            
            if not gate_passed and stage_name not in ("equivalence", "source_mining"):
                # v2: If source mining is enabled, comprehension failure is non-fatal
                # (source mining uses regex + known_callees, not comprehension output)
                if stage_name == "comprehension" and self.enable_source_mining:
                    self.logger.log_info(f"  Comprehension failed but source mining enabled — continuing")
                    continue
                self.logger.log_info(f"Pipeline halted at stage: {stage_name}")
                break
        
        # Save final state
        self._save_state()
        
        return self.state
    
    def _inject_kb_context(self, benchmark_config: dict):
        """Query knowledge base and inject relevant context into state."""
        # Get called subroutines from all source files
        import re
        all_calls = []
        for code in self.state.source_files.values():
            for match in re.finditer(r'^\s*CALL\s+(\w+)', code, re.MULTILINE | re.IGNORECASE):
                name = match.group(1).upper()
                if name not in all_calls:
                    all_calls.append(name)
        
        insights = self.knowledge_base.query_for_benchmark(
            benchmark_config, called_subroutines=all_calls
        )
        
        if insights:
            self.state.kb_context = self.knowledge_base.format_context(insights, max_length=4000)
            self.logger.log_info(f"  KB: Found {len(insights)} relevant insights")
            self.logger.log_artifact("kb_context", {
                "num_insights": len(insights),
                "sources": [i.source_benchmark for i in insights],
                "context_preview": self.state.kb_context[:300],
            })
        else:
            self.logger.log_info("  KB: No relevant insights found")
    
    def _run_comprehension(self, context: str):
        """Stage 1: Understand the Fortran code."""
        self.state.stage = "comprehension"
        self.state.timestamps["comprehension_start"] = time.time()
        
        intents = []
        for filename, code in self.state.source_files.items():
            self.logger.log_info(f"Analyzing: {filename}")
            intent = self.comprehension_agent.analyze_subroutine(code, filename)
            intents.append(intent.to_dict())
            self.logger.log_info(f"  → Domain: {intent.physical_domain}")
            self.logger.log_info(f"  → Operation: {intent.physical_operation}")
            self.logger.log_debug(f"  → Variables: {len(intent.input_variables)} in, {len(intent.output_variables)} out")
        
        self.state.code_intents = intents
        self.state.timestamps["comprehension_end"] = time.time()
        
        self.logger.log_artifact("code_intents", intents)
    
    def _run_source_mining(self, context: str):
        """Stage 1.5 (v2): Mine annotations from called subroutines' source code.
        
        This stage addresses the key failure mode: LLMs miss implicit conventions
        (storage order, non-standard algorithms) that are documented only in callees.
        """
        self.state.stage = "source_mining"
        self.state.timestamps["source_mining_start"] = time.time()
        
        if not self.source_mining_agent:
            self.logger.log_info("Source mining skipped (agent not initialized)")
            self.state.timestamps["source_mining_end"] = time.time()
            return
        
        # Gather called subroutines from comprehension results
        called_subs = []
        for intent in self.state.code_intents:
            for sub in intent.get("called_subroutines", []):
                if sub.upper() not in called_subs:
                    called_subs.append(sub.upper())
        
        # Also merge in known_callees from benchmark config (if provided)
        if self._benchmark_config:
            for sub in self._benchmark_config.get("known_callees", []):
                if sub.upper() not in called_subs:
                    called_subs.append(sub.upper())
        
        # Fallback: if comprehension found nothing, extract directly via regex
        if not called_subs:
            import re as _re
            for code in self.state.source_files.values():
                for match in _re.finditer(r'^\s*CALL\s+(\w+)', code, _re.MULTILINE | _re.IGNORECASE):
                    name = match.group(1).upper()
                    if name not in called_subs:
                        called_subs.append(name)
        
        self.logger.log_info(f"Mining annotations for {len(called_subs)} called subroutines: {called_subs}")
        
        # Get the primary source code (concatenated)
        primary_code = '\n\n'.join(self.state.source_files.values())
        primary_filename = list(self.state.source_files.keys())[0] if self.state.source_files else "unknown"
        
        # Run mining
        mining_result = self.source_mining_agent.mine_annotations(
            primary_code=primary_code,
            primary_filename=primary_filename,
            called_subroutines=called_subs,
        )
        
        # Store results
        self.state.source_mining = mining_result.to_dict()
        
        # Log summary
        self.logger.log_info(f"  Resolved: {len(mining_result.resolved_sources)} sources")
        self.logger.log_info(f"  Unresolved: {mining_result.unresolved}")
        self.logger.log_info(f"  Storage warnings: {len(mining_result.storage_warnings)}")
        self.logger.log_info(f"  Critical conventions: {len(mining_result.critical_conventions)}")
        
        for w in mining_result.storage_warnings:
            self.logger.log_info(f"  ⚠️ {w}")
        for c in mining_result.critical_conventions:
            self.logger.log_info(f"  📋 {c}")
        
        self.state.timestamps["source_mining_end"] = time.time()
        self.logger.log_artifact("source_mining", self.state.source_mining)
    
    def _run_architecture_recovery(self, context: str):
        """Stage 2: Recover architecture."""
        self.state.stage = "architecture"
        self.state.timestamps["architecture_start"] = time.time()
        
        arch = self.architecture_agent.recover_architecture(
            self.state.source_files, context
        )
        self.state.architecture = arch.to_dict()
        
        self.logger.log_info(f"Call graph: {len(arch.call_graph)} subroutines")
        self.logger.log_info(f"COMMON blocks: {len(arch.common_blocks)}")
        self.logger.log_info(f"Refactoring order: {arch.refactoring_order}")
        
        self.state.timestamps["architecture_end"] = time.time()
        self.logger.log_artifact("architecture", self.state.architecture)
    
    def _run_physics_specification(self, context: str):
        """Stage 3: Generate physics specification."""
        self.state.stage = "physics_spec"
        self.state.timestamps["physics_spec_start"] = time.time()
        
        primary_intent = self.state.code_intents[0] if self.state.code_intents else {}
        
        # v2: Enrich context with mining results and KB insights
        enriched_context = context
        
        # Append source mining conventions
        mining = self.state.source_mining
        if mining:
            conventions = mining.get("critical_conventions", [])
            storage_warnings = mining.get("storage_warnings", [])
            if conventions or storage_warnings:
                enriched_context += "\n\n--- CRITICAL CONVENTIONS FROM CALLED SUBROUTINES ---\n"
                for w in storage_warnings:
                    enriched_context += f"⚠️ STORAGE: {w}\n"
                for c in conventions:
                    enriched_context += f"• {c}\n"
                enriched_context += "--- END CONVENTIONS ---\n"
        
        # Append KB context
        if self.state.kb_context:
            enriched_context += "\n\n" + self.state.kb_context
        
        spec = self.physics_spec_agent.generate_specification(primary_intent, enriched_context)
        self.state.physics_spec = spec.to_dict()
        
        self.logger.log_info(f"Domain: {spec.domain}")
        self.logger.log_info(f"Element: {spec.element_type}")
        self.logger.log_info(f"DOF: {spec.num_nodes} nodes × {spec.dof_per_node} DOF = {spec.total_dof}")
        self.logger.log_info(f"Constitutive: {spec.constitutive_law}")
        
        self.state.timestamps["physics_spec_end"] = time.time()
        self.logger.log_artifact("physics_spec", self.state.physics_spec)
    
    def _run_modern_design(self):
        """Stage 4: Design modern implementation."""
        self.state.stage = "modern_design"
        self.state.timestamps["modern_design_start"] = time.time()
        
        original_code = list(self.state.source_files.values())[0] if self.state.source_files else ""
        
        design = self.modern_design_agent.design(
            physics_spec=self.state.physics_spec,
            architecture=self.state.architecture,
            original_code=original_code
        )
        self.state.modern_design = design.to_dict()
        
        self.logger.log_info(f"Module: {design.module_name}")
        self.logger.log_info(f"Classes: {len(design.classes)}")
        self.logger.log_info(f"Functions: {len(design.functions)}")
        self.logger.log_info(f"Design decisions: {len(design.design_decisions)}")
        for i, decision in enumerate(design.design_decisions[:5]):
            self.logger.log_debug(f"  Decision {i+1}: {decision}")
        
        self.state.timestamps["modern_design_end"] = time.time()
        self.logger.log_artifact("modern_design", self.state.modern_design)
    
    def _run_code_generation(self):
        """Stage 5: Generate Python code."""
        self.state.stage = "code_generation"
        self.state.timestamps["code_gen_start"] = time.time()
        
        original_code = list(self.state.source_files.values())[0] if self.state.source_files else ""
        
        result = self.code_gen_agent.generate(
            design=self.state.modern_design,
            physics_spec=self.state.physics_spec,
            original_fortran=original_code
        )
        
        self.state.generated_code = result.code
        self.state.iterations["code_generation"] = result.iteration
        
        self.logger.log_info(f"Generation success: {result.success}")
        self.logger.log_info(f"Iterations used: {result.iteration}")
        self.logger.log_info(f"Code length: {len(result.code)} chars")
        if result.error_message:
            self.logger.log_info(f"Error: {result.error_message}")
        
        self.state.timestamps["code_gen_end"] = time.time()
        
        # Save generated code as artifact
        self.logger.log_artifact("generated_code.py", result.code, is_json=False)
        if result.metrics:
            self.logger.log_artifact("generation_metrics", result.metrics)
    
    def _run_equivalence_verification(self, reference_results: list[dict] = None):
        """Stage 6: Verify equivalence."""
        self.state.stage = "equivalence"
        self.state.timestamps["equivalence_start"] = time.time()
        
        if reference_results is None:
            self.logger.log_info("Generating test suite from physics spec...")
            reference_results = self.equivalence_agent.generate_test_suite(
                self.state.physics_spec
            )
            self.logger.log_info(f"Generated {len(reference_results)} test cases")
        else:
            self.logger.log_info(f"Using {len(reference_results)} provided reference results")
        
        report = self.equivalence_agent.verify_element_stiffness(
            self.state.generated_code,
            reference_results
        )
        
        self.state.equivalence_report = report.to_dict()
        
        self.logger.log_info(f"Equivalent: {report.equivalent}")
        self.logger.log_info(f"Tests passed: {report.test_cases_passed}/{report.test_cases_total}")
        if report.numerical_checks:
            self.logger.log_info(f"Max error: {report.numerical_checks.get('max_error_across_tests', 'N/A')}")
            self.logger.log_info(f"Within tolerance: {report.numerical_checks.get('within_tolerance', 'N/A')}")
        if report.failed_cases:
            self.logger.log_info(f"Failed cases:")
            for fc in report.failed_cases:
                self.logger.log_info(f"  - {fc.get('description', 'unknown')}: {fc.get('error', '')}")
        
        self.state.timestamps["equivalence_end"] = time.time()
        self.logger.log_artifact("equivalence_report", self.state.equivalence_report)
    
    def _check_quality_gate(self, stage: str) -> bool:
        """Check quality gate for a stage."""
        passed = True
        reason = ""
        
        if stage == "comprehension":
            if len(self.state.code_intents) == 0:
                passed = False
                reason = "No code intents extracted"
            elif self.state.code_intents[0].get("physical_domain") == "unknown" and \
                 self.state.code_intents[0].get("physical_operation") == "unknown":
                passed = False
                reason = "No useful information extracted"
        
        elif stage == "source_mining":
            # Source mining is best-effort: always pass (don't halt pipeline)
            # But log whether it found anything useful
            mining = self.state.source_mining
            has_conventions = bool(mining.get("critical_conventions") or mining.get("storage_warnings"))
            if has_conventions:
                self.logger.log_info("  Source mining: found useful conventions")
            else:
                self.logger.log_info("  Source mining: no critical conventions found (proceeding anyway)")
            passed = True
        
        elif stage == "architecture":
            if not self.state.architecture.get("call_graph"):
                passed = False
                reason = "Call graph is empty"
        
        elif stage == "physics_spec":
            spec = self.state.physics_spec
            if spec.get("domain") in ("unknown", "parse_error") and \
               spec.get("element_type") == "unknown" and \
               spec.get("num_nodes", 0) <= 0:
                passed = False
                reason = "No useful physics specification extracted"
        
        elif stage == "modern_design":
            # Pass if we got any useful design info (classes OR design decisions)
            has_classes = len(self.state.modern_design.get("classes", [])) > 0
            has_decisions = len(self.state.modern_design.get("design_decisions", [])) > 0
            has_module = self.state.modern_design.get("module_name", "parse_error") != "parse_error"
            passed = has_classes or (has_decisions and has_module)
            if not passed:
                reason = "No useful design information extracted"
        
        elif stage == "code_generation":
            if len(self.state.generated_code) < 100:
                passed = False
                reason = f"Generated code too short ({len(self.state.generated_code)} chars)"
        
        elif stage == "equivalence":
            passed = self.state.equivalence_report.get("equivalent", False)
            if not passed:
                failed = self.state.equivalence_report.get("test_cases_total", 0) - \
                         self.state.equivalence_report.get("test_cases_passed", 0)
                reason = f"{failed} test case(s) failed"
        
        self.state.quality_gates[stage] = passed
        self.logger.log_quality_gate(stage, passed, reason)
        
        return passed
    
    def _save_state(self):
        """Save the complete workflow state."""
        self.logger.log_artifact("workflow_state", self.state.to_dict())
        self.logger.log_info(f"Final state saved to artifacts/workflow_state.json")
