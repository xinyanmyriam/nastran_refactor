"""
Knowledge Base for Cross-Benchmark Learning.

Accumulates insights from successful repairs and provides them to subsequent
experiments. Implements the key finding from Level 3 experiments: once an
implicit convention is discovered (e.g., GMMATD row-major storage), it can
be reused across all related benchmarks at near-zero marginal cost.

Storage: JSON file persisted to disk, loaded at experiment start.

Design:
- Each "insight" is a structured record of what was learned from a successful run.
- Insights are indexed by: called_subroutines, category, and keywords.
- Query returns relevant insights ranked by similarity to the current benchmark.
"""

import json
import re
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional


@dataclass
class Insight:
    """A single knowledge insight from a successful (or informative) run."""
    
    # Identity
    id: str  # Unique ID (benchmark_model_timestamp)
    source_benchmark: str  # e.g., "B5"
    source_subroutine: str  # e.g., "KTRPLT"
    
    # What was learned
    category: str  # "storage_convention" | "algorithm" | "transform" | "indexing" | "general"
    description: str  # Human-readable description of the insight
    
    # Context: when does this insight apply?
    applies_to_subroutines: list[str] = field(default_factory=list)  # Called subs this is about
    applies_to_categories: list[str] = field(default_factory=list)  # e.g., ["plate_bending", "solid_3d"]
    keywords: list[str] = field(default_factory=list)  # Searchable keywords
    
    # The actual knowledge
    convention: str = ""  # e.g., "GMMATD stores output by rows"
    code_pattern: str = ""  # Example code showing correct handling
    failure_mode: str = ""  # What goes wrong if this isn't respected
    
    # Metadata
    confidence: float = 1.0  # 0-1, how reliable this insight is
    timestamp: float = 0.0
    model: str = ""
    method: str = ""
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, d: dict) -> "Insight":
        # Filter to valid fields only
        valid_fields = set(cls.__dataclass_fields__.keys())
        filtered = {k: v for k, v in d.items() if k in valid_fields}
        return cls(**filtered)


class KnowledgeBase:
    """Persistent knowledge store with similarity-based retrieval."""
    
    def __init__(self, storage_path: str = None):
        """
        Args:
            storage_path: Path to JSON file for persistence.
                Defaults to project_root/knowledge_base.json
        """
        if storage_path is None:
            storage_path = str(
                Path(__file__).parent / "knowledge_base.json"
            )
        self.storage_path = Path(storage_path)
        self.insights: list[Insight] = []
        self._load()
    
    def _load(self):
        """Load insights from disk."""
        if self.storage_path.exists():
            try:
                data = json.loads(
                    self.storage_path.read_text(encoding='utf-8')
                )
                self.insights = [
                    Insight.from_dict(d) for d in data.get("insights", [])
                ]
            except (json.JSONDecodeError, KeyError, TypeError):
                self.insights = []
    
    def _save(self):
        """Persist insights to disk."""
        data = {
            "version": "1.0",
            "last_updated": time.strftime("%Y-%m-%d %H:%M:%S"),
            "total_insights": len(self.insights),
            "insights": [ins.to_dict() for ins in self.insights],
        }
        self.storage_path.parent.mkdir(parents=True, exist_ok=True)
        self.storage_path.write_text(
            json.dumps(data, indent=2, ensure_ascii=False),
            encoding='utf-8'
        )
    
    def add_insight(self, insight: Insight) -> None:
        """Add a new insight and persist."""
        # Check for duplicates (same source + same convention)
        for existing in self.insights:
            if (existing.source_benchmark == insight.source_benchmark and
                existing.convention == insight.convention and
                existing.category == insight.category):
                # Update confidence if re-discovered
                existing.confidence = min(1.0, existing.confidence + 0.1)
                self._save()
                return
        
        if not insight.timestamp:
            insight.timestamp = time.time()
        self.insights.append(insight)
        self._save()
    
    def add_from_successful_run(
        self,
        benchmark_id: str,
        subroutine_name: str,
        model: str,
        method: str,
        category: str,
        called_subroutines: list[str],
        repair_history: list[dict] = None,
        final_code: str = "",
        mining_result: dict = None,
    ) -> list[Insight]:
        """
        Extract and store insights from a successful experiment run.
        
        Args:
            benchmark_id: e.g., "B5"
            subroutine_name: e.g., "KTRPLT"
            model: LLM model used
            method: experiment method
            category: benchmark category (e.g., "A4_plate_bending")
            called_subroutines: list of called routines
            repair_history: list of repair iterations with error info
            final_code: the successful C++ code
            mining_result: dict from source mining (if available)
        
        Returns:
            List of newly created insights.
        """
        new_insights = []
        timestamp = time.time()
        
        # Insight from mining results (critical conventions discovered)
        if mining_result:
            conventions = mining_result.get("critical_conventions", [])
            storage_warnings = mining_result.get("storage_warnings", [])
            
            for conv in conventions + storage_warnings:
                insight = Insight(
                    id=f"{benchmark_id}_{model}_{int(timestamp)}_{len(new_insights)}",
                    source_benchmark=benchmark_id,
                    source_subroutine=subroutine_name,
                    category="storage_convention" if "storage" in conv.lower() or "row" in conv.lower() else "general",
                    description=f"Convention from {subroutine_name}: {conv}",
                    applies_to_subroutines=called_subroutines,
                    applies_to_categories=[category],
                    keywords=self._extract_keywords(conv),
                    convention=conv,
                    failure_mode="Numerical mismatch due to incorrect matrix interpretation",
                    confidence=0.9,
                    timestamp=timestamp,
                    model=model,
                    method=method,
                )
                new_insights.append(insight)
        
        # Insight from repair history (what errors were fixed)
        if repair_history:
            for i, repair in enumerate(repair_history):
                error_type = repair.get("stage", "unknown")
                error_msg = repair.get("error", "")[:200]
                
                if error_type == "numerical" and "relative error" in error_msg.lower():
                    insight = Insight(
                        id=f"{benchmark_id}_{model}_{int(timestamp)}_repair{i}",
                        source_benchmark=benchmark_id,
                        source_subroutine=subroutine_name,
                        category="transform" if "transform" in error_msg.lower() else "indexing",
                        description=f"Numerical repair for {subroutine_name}: fixed via iteration {i+1}",
                        applies_to_subroutines=called_subroutines,
                        applies_to_categories=[category],
                        keywords=self._extract_keywords(error_msg),
                        failure_mode=error_msg,
                        confidence=0.7,
                        timestamp=timestamp,
                        model=model,
                        method=method,
                    )
                    new_insights.append(insight)
        
        # Insight from successful code patterns
        if final_code and len(final_code) > 100:
            # Extract key patterns from successful code
            patterns = self._extract_code_patterns(final_code, called_subroutines)
            for pattern_desc, pattern_code in patterns:
                insight = Insight(
                    id=f"{benchmark_id}_{model}_{int(timestamp)}_pattern{len(new_insights)}",
                    source_benchmark=benchmark_id,
                    source_subroutine=subroutine_name,
                    category="algorithm",
                    description=f"Successful pattern from {subroutine_name}: {pattern_desc}",
                    applies_to_subroutines=called_subroutines,
                    applies_to_categories=[category],
                    keywords=self._extract_keywords(pattern_desc),
                    code_pattern=pattern_code[:500],
                    confidence=0.8,
                    timestamp=timestamp,
                    model=model,
                    method=method,
                )
                new_insights.append(insight)
        
        # Store all new insights
        for ins in new_insights:
            self.add_insight(ins)
        
        return new_insights
    
    def query(
        self,
        called_subroutines: list[str] = None,
        category: str = None,
        keywords: list[str] = None,
        max_results: int = 5,
    ) -> list[Insight]:
        """
        Query the knowledge base for relevant insights.
        
        Scoring: each insight gets a relevance score based on overlap
        with the query parameters. Higher = more relevant.
        
        Args:
            called_subroutines: subroutines called by the current benchmark
            category: benchmark category
            keywords: additional search keywords
            max_results: maximum number of insights to return
        
        Returns:
            List of relevant insights, sorted by score descending.
        """
        if not self.insights:
            return []
        
        scored = []
        for ins in self.insights:
            score = 0.0
            
            # Subroutine overlap (highest signal)
            if called_subroutines:
                query_subs = set(s.upper() for s in called_subroutines)
                ins_subs = set(s.upper() for s in ins.applies_to_subroutines)
                overlap = len(query_subs & ins_subs)
                if overlap > 0:
                    score += 3.0 * overlap
            
            # Category match
            if category and category in ins.applies_to_categories:
                score += 2.0
            
            # Keyword overlap
            if keywords:
                query_kw = set(k.lower() for k in keywords)
                ins_kw = set(k.lower() for k in ins.keywords)
                kw_overlap = len(query_kw & ins_kw)
                if kw_overlap > 0:
                    score += 1.0 * kw_overlap
            
            # Confidence weight
            score *= ins.confidence
            
            if score > 0:
                scored.append((score, ins))
        
        # Sort by score descending
        scored.sort(key=lambda x: x[0], reverse=True)
        
        return [ins for _, ins in scored[:max_results]]
    
    def query_for_benchmark(
        self,
        benchmark_config: dict,
        called_subroutines: list[str] = None,
    ) -> list[Insight]:
        """
        Convenience method: query using a benchmark config dict.
        
        Args:
            benchmark_config: dict from benchmarks_config.json
            called_subroutines: optional override for called subroutines
        """
        category = benchmark_config.get("category", "")
        name = benchmark_config.get("name", "").upper()
        description = benchmark_config.get("description", "")
        
        # Extract keywords from description
        keywords = self._extract_keywords(description)
        keywords.append(name)
        
        return self.query(
            called_subroutines=called_subroutines,
            category=category,
            keywords=keywords,
        )
    
    def format_context(self, insights: list[Insight], max_length: int = 2000) -> str:
        """
        Format insights as context string for injection into prompts.
        
        Args:
            insights: list of relevant insights
            max_length: maximum character length of output
        
        Returns:
            Formatted context string.
        """
        if not insights:
            return ""
        
        lines = ["=== KNOWLEDGE BASE: Relevant Insights from Prior Experiments ==="]
        
        for i, ins in enumerate(insights, 1):
            lines.append(f"\n[{i}] From {ins.source_benchmark}/{ins.source_subroutine} (confidence: {ins.confidence:.1f}):")
            lines.append(f"    {ins.description}")
            
            if ins.convention:
                lines.append(f"    CONVENTION: {ins.convention}")
            
            if ins.failure_mode:
                lines.append(f"    If ignored: {ins.failure_mode}")
            
            if ins.code_pattern:
                # Show first few lines of code pattern
                pattern_lines = ins.code_pattern.strip().split('\n')[:5]
                lines.append(f"    Code pattern:")
                for pl in pattern_lines:
                    lines.append(f"      {pl}")
            
            # Check length
            current = '\n'.join(lines)
            if len(current) > max_length:
                lines.append("\n... (additional insights truncated)")
                break
        
        lines.append("\n=== END KNOWLEDGE BASE ===")
        return '\n'.join(lines)
    
    def _extract_keywords(self, text: str) -> list[str]:
        """Extract searchable keywords from text."""
        # Common FEM/NASTRAN keywords to look for
        fem_keywords = {
            'stiffness', 'mass', 'matrix', 'element', 'node', 'dof',
            'plate', 'shell', 'beam', 'rod', 'solid', 'membrane',
            'bending', 'torsion', 'axial', 'shear',
            'triangle', 'quad', 'tetra', 'hexa', 'wedge',
            'transform', 'rotation', 'transpose',
            'row', 'column', 'storage', 'packed', 'symmetric',
            'gauss', 'integration', 'jacobian',
            'gmmatd', 'gmmats', 'transd', 'ktrbsc',
        }
        
        words = set(re.findall(r'\b\w+\b', text.lower()))
        return list(words & fem_keywords)
    
    def _extract_code_patterns(
        self, code: str, called_subroutines: list[str]
    ) -> list[tuple[str, str]]:
        """
        Extract notable code patterns from successful implementations.
        
        Returns list of (description, code_snippet) tuples.
        """
        patterns = []
        
        # Look for matrix storage handling
        storage_patterns = [
            (r'(?:row[\s_]*major|RowMajor)', "Row-major storage handling"),
            (r'(?:col[\s_]*major|ColMajor)', "Column-major storage handling"),
            (r'\.transpose\(\)', "Explicit transpose operation"),
            (r'Map<.*RowMajor', "Eigen RowMajor map"),
        ]
        
        for pattern, desc in storage_patterns:
            match = re.search(pattern, code, re.IGNORECASE)
            if match:
                # Get surrounding context (5 lines)
                start = code.rfind('\n', 0, max(0, match.start() - 100))
                end = code.find('\n', match.end() + 100)
                if start == -1:
                    start = 0
                if end == -1:
                    end = len(code)
                snippet = code[start:end].strip()
                patterns.append((desc, snippet))
        
        # Look for condensation / sub-matrix assembly patterns
        if any(s in code for s in ['condense', 'static_condensation', 'submatrix']):
            match = re.search(r'(?:condense|static_condensation|submatrix)', code, re.IGNORECASE)
            if match:
                start = code.rfind('\n', 0, max(0, match.start() - 150))
                end = code.find('\n', match.end() + 150)
                snippet = code[max(0, start):min(len(code), end)].strip()
                patterns.append(("Static condensation / sub-matrix assembly", snippet))
        
        return patterns[:3]  # Limit to 3 patterns
    
    def get_stats(self) -> dict:
        """Get summary statistics of the knowledge base."""
        if not self.insights:
            return {"total": 0}
        
        categories = {}
        benchmarks = set()
        subroutines = set()
        
        for ins in self.insights:
            categories[ins.category] = categories.get(ins.category, 0) + 1
            benchmarks.add(ins.source_benchmark)
            subroutines.update(ins.applies_to_subroutines)
        
        return {
            "total": len(self.insights),
            "categories": categories,
            "source_benchmarks": sorted(benchmarks),
            "covered_subroutines": sorted(subroutines),
            "avg_confidence": sum(i.confidence for i in self.insights) / len(self.insights),
        }
    
    def clear(self):
        """Clear all insights (for fresh experiments)."""
        self.insights = []
        self._save()
    
    def seed_from_level3(self) -> None:
        """
        Seed the knowledge base with insights from Level 3 experiments.
        These are the critical discoveries that took 8+ rounds to find.
        """
        seeds = [
            Insight(
                id="seed_gmmatd_rowmajor",
                source_benchmark="B5",
                source_subroutine="KTRPLT",
                category="storage_convention",
                description="GMMATD stores output matrix by ROWS, not Fortran-default column-major",
                applies_to_subroutines=["GMMATD", "GMMATS"],
                applies_to_categories=["A4_plate_bending", "A5_solid_3d", "A6_special"],
                keywords=["gmmatd", "row", "storage", "matrix", "transpose"],
                convention="GMMATD output is stored row-by-row. When interpreting as a 2D matrix, "
                           "read consecutive elements as filling rows, not columns. This effectively "
                           "means the result appears transposed compared to standard Fortran column-major layout.",
                failure_mode="Matrix appears transposed → stiffness matrix is wrong → numerical mismatch",
                confidence=1.0,
                timestamp=time.time(),
                model="level3_human_ai",
                method="level3",
            ),
            Insight(
                id="seed_ktrbsc_condensation",
                source_benchmark="B5",
                source_subroutine="KTRPLT",
                category="algorithm",
                description="KTRPLT uses FMMS-55 sub-triangle condensation via KTRBSC, not standard DKT",
                applies_to_subroutines=["KTRBSC", "KTRPLT", "KQDPLT"],
                applies_to_categories=["A4_plate_bending"],
                keywords=["plate", "bending", "triangle", "condensation", "ktrbsc"],
                convention="Plate bending uses 3 sub-triangles (KTRBSC computes each), then "
                           "condenses internal DOFs to get the final 9x9 or 12x12 matrix. "
                           "Do NOT use standard DKT or Mindlin plate formulation.",
                code_pattern="// For each sub-triangle:\n"
                            "//   1. Compute 6x6 basic bending stiffness via KTRBSC\n"
                            "//   2. Assemble into 15x15 expanded matrix\n"
                            "//   3. Static condensation: K_final = K_ee - K_ei * K_ii^-1 * K_ie",
                failure_mode="Using textbook DKT gives completely different matrix → 0% match",
                confidence=1.0,
                timestamp=time.time(),
                model="level3_human_ai",
                method="level3",
            ),
            Insight(
                id="seed_fortran_colmajor_implicit_transpose",
                source_benchmark="B9",
                source_subroutine="KELBOW",
                category="transform",
                description="Passing Fortran 2D array to GMMATD causes implicit transpose due to column-major storage",
                applies_to_subroutines=["GMMATD", "GMMATS"],
                applies_to_categories=["A4_plate_bending", "A5_solid_3d", "A6_special"],
                keywords=["transpose", "column", "row", "gmmatd", "transform"],
                convention="When Fortran passes a 2D array A(M,N) to GMMATD which reads it "
                           "as a 1D array stored by rows, the effective matrix seen by GMMATD "
                           "is A^T. This implicit transpose must be accounted for in C++ translation.",
                failure_mode="Missing implicit transpose → coordinate transform wrong → wrong stiffness",
                confidence=1.0,
                timestamp=time.time(),
                model="level3_human_ai",
                method="level3",
            ),
            Insight(
                id="seed_ksolid_ktetra_decomposition",
                source_benchmark="B8",
                source_subroutine="KSOLID",
                category="algorithm",
                description="KSOLID decomposes wedge/hex into tetrahedra using KTETRA with averaging weights",
                applies_to_subroutines=["KTETRA", "KSOLID"],
                applies_to_categories=["A5_solid_3d"],
                keywords=["solid", "tetra", "wedge", "element", "stiffness"],
                convention="KSOLID does not use isoparametric hex formulation. Instead it "
                           "decomposes the element into tetrahedra (3 for wedge, 5 for hex), "
                           "computes each via KTETRA, and averages with specific weights.",
                failure_mode="Using isoparametric formulation gives different result",
                confidence=1.0,
                timestamp=time.time(),
                model="level3_human_ai",
                method="level3",
            ),
            Insight(
                id="seed_inverd_gauss_jordan",
                source_benchmark="B10",
                source_subroutine="INVERD",
                category="algorithm",
                description="INVERD is standard Gauss-Jordan elimination with partial pivoting - straightforward",
                applies_to_subroutines=["INVERD"],
                applies_to_categories=["A7_solver"],
                keywords=["matrix", "gauss", "solver", "inverse"],
                convention="INVERD implements standard Gauss-Jordan in-place inversion. "
                           "No special conventions. Can use Eigen's .inverse() directly.",
                confidence=1.0,
                timestamp=time.time(),
                model="level3_human_ai",
                method="level3",
            ),
        ]
        
        for seed in seeds:
            self.add_insight(seed)
