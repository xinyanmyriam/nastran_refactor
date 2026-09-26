"""
Source Annotation Mining Agent.

Resolves called subroutines from the NASTRAN-95 source tree, reads their
comments and documentation, and identifies implicit conventions that are
critical for correct refactoring:
- Storage conventions (row-major, column-major)
- Non-standard algorithms referenced
- Implicit type conversions / EQUIVALENCE tricks
- Data format expectations

This addresses the key failure mode discovered in Level 3 experiments:
LLMs fail because they cannot discover implicit conventions hidden in
callee documentation that are never explicitly stated in the caller.
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class SourceAnnotation:
    """Annotations extracted from a called subroutine's source."""
    subroutine_name: str
    source_file: str
    header_comments: str  # First block of comments (documentation)
    storage_conventions: list[str]  # e.g., "stored by rows", "column-major"
    algorithm_references: list[str]  # e.g., "FMMS-55", "Timoshenko"
    data_format_notes: list[str]  # e.g., "upper triangle only"
    key_parameters: dict[str, str]  # parameter name -> documented meaning
    warnings: list[str]  # any "NOTE", "WARNING", "CAUTION" comments


@dataclass
class MiningResult:
    """Complete result of source annotation mining for a benchmark."""
    primary_subroutine: str
    called_subroutines: list[str]
    resolved_sources: dict[str, str]  # sub_name -> file_path
    unresolved: list[str]  # subroutines whose source was not found
    annotations: list[SourceAnnotation]
    
    # LLM-synthesized summary of critical conventions
    critical_conventions: list[str] = field(default_factory=list)
    storage_warnings: list[str] = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "primary_subroutine": self.primary_subroutine,
            "called_subroutines": self.called_subroutines,
            "resolved_sources": self.resolved_sources,
            "unresolved": self.unresolved,
            "annotations": [
                {
                    "subroutine_name": a.subroutine_name,
                    "source_file": a.source_file,
                    "header_comments": a.header_comments[:500],
                    "storage_conventions": a.storage_conventions,
                    "algorithm_references": a.algorithm_references,
                    "data_format_notes": a.data_format_notes,
                    "key_parameters": a.key_parameters,
                    "warnings": a.warnings,
                }
                for a in self.annotations
            ],
            "critical_conventions": self.critical_conventions,
            "storage_warnings": self.storage_warnings,
        }
    
    def to_context_string(self) -> str:
        """Format as a context string to inject into downstream agents."""
        lines = ["=== SOURCE ANNOTATION MINING RESULTS ==="]
        
        if self.storage_warnings:
            lines.append("\n⚠️ CRITICAL STORAGE CONVENTIONS:")
            for w in self.storage_warnings:
                lines.append(f"  - {w}")
        
        if self.critical_conventions:
            lines.append("\n📋 CONVENTIONS FROM CALLED SUBROUTINES:")
            for c in self.critical_conventions:
                lines.append(f"  - {c}")
        
        for ann in self.annotations:
            if ann.storage_conventions or ann.algorithm_references or ann.warnings:
                lines.append(f"\n--- {ann.subroutine_name} ({ann.source_file}) ---")
                if ann.header_comments:
                    # Show first 3 lines of comments
                    comment_lines = ann.header_comments.strip().split('\n')[:3]
                    for cl in comment_lines:
                        lines.append(f"  {cl.strip()}")
                if ann.storage_conventions:
                    lines.append(f"  Storage: {'; '.join(ann.storage_conventions)}")
                if ann.algorithm_references:
                    lines.append(f"  Algorithm: {'; '.join(ann.algorithm_references)}")
                if ann.warnings:
                    for w in ann.warnings:
                        lines.append(f"  ⚠️ {w}")
        
        if self.unresolved:
            lines.append(f"\n❓ Unresolved callees: {', '.join(self.unresolved)}")
        
        lines.append("=== END MINING RESULTS ===")
        return '\n'.join(lines)


MINING_SYSTEM_PROMPT = """You are a Source Annotation Mining Agent specialized in extracting 
implicit conventions from legacy Fortran 77 subroutines in NASA NASTRAN-95.

Your task is to analyze comments and documentation from called subroutines and identify
CRITICAL IMPLICIT CONVENTIONS that a refactoring agent needs to know, especially:

1. **Storage conventions**: How matrices/arrays are stored (row-major, column-major,
   upper triangle, packed, banded). Look for keywords like "stored by rows", 
   "column major", "packed", "upper triangle", "lower triangle".

2. **Algorithm references**: Named algorithms or technical reports (e.g., "FMMS-55",
   "Timoshenko beam theory", "Mindlin plate"). These indicate non-standard implementations.

3. **Data format expectations**: How input/output arrays are organized. E.g., whether
   a matrix multiplication routine expects row-major or column-major input.

4. **Implicit transformations**: Operations that happen as a side effect of data layout.
   E.g., passing a Fortran column-major 2D array to a routine expecting row-major 
   effectively transposes the matrix.

NASTRAN-95 specific patterns to watch for:
- GMMATD/GMMATS: General matrix multiply routines with row-storage conventions
- TRANSD/TRANSS: Transpose routines with specific storage assumptions  
- FORMAA/FORMBB: Assembly routines with block-offset conventions
- KTRPLT/KTRBSC: Plate bending with sub-triangle condensation (FMMS-55)

Output your analysis as JSON:
{
    "critical_conventions": ["convention 1 that MUST be respected", ...],
    "storage_warnings": ["WARNING: GMMATD stores result by rows, not Fortran default column-major", ...],
    "algorithm_notes": ["Non-standard algorithm X referenced in subroutine Y", ...]
}
"""


# Common NASTRAN utility subroutines that we know exist in mis/
KNOWN_NASTRAN_UTILITIES = {
    'GMMATD', 'GMMATS', 'TRANSD', 'TRANSS', 'FORMAA', 'FORMBB',
    'MATADD', 'MATSUB', 'MATSCL', 'MATMPY', 'INVERS', 'INVERD',
    'DAXPY', 'DDOT', 'DNRM2', 'DSCAL', 'DCOPY',
    'SMA1A', 'SMA1B', 'SMA2A', 'SMA2B',
    'ANISOP', 'CROSS', 'DOT', 'LENGTH',
}


class SourceMiningAgent:
    """Agent that mines annotations from called subroutines' source code."""
    
    def __init__(self, llm_client, nastran_source_dir: str = None):
        """
        Args:
            llm_client: LLM client for semantic analysis.
            nastran_source_dir: Path to NASTRAN-95 mis/ directory.
        """
        self.llm_client = llm_client
        self.system_prompt = MINING_SYSTEM_PROMPT
        
        if nastran_source_dir:
            self.source_dir = Path(nastran_source_dir)
        else:
            # Default: relative to project root
            self.source_dir = Path(__file__).parent.parent / "nastran" / "NASTRAN-95" / "mis"
        
        # Build index of available source files
        self._source_index = self._build_source_index()
    
    def _build_source_index(self) -> dict[str, Path]:
        """Build mapping: subroutine_name (upper) -> file path."""
        index = {}
        if not self.source_dir.exists():
            return index
        
        for f in self.source_dir.glob("*.f"):
            # Index by filename (without extension, uppercased)
            name = f.stem.upper()
            index[name] = f
            
            # Also scan for SUBROUTINE declarations in the file
            try:
                content = f.read_text(encoding='utf-8', errors='replace')
                for match in re.finditer(
                    r'^\s*SUBROUTINE\s+(\w+)', content, re.MULTILINE | re.IGNORECASE
                ):
                    sub_name = match.group(1).upper()
                    if sub_name not in index:
                        index[sub_name] = f
            except Exception:
                pass
        
        return index
    
    def mine_annotations(
        self,
        primary_code: str,
        primary_filename: str,
        called_subroutines: list[str] = None,
    ) -> MiningResult:
        """
        Mine annotations from all called subroutines.
        
        Args:
            primary_code: The main Fortran source being refactored.
            primary_filename: Filename of the primary source.
            called_subroutines: List of called routine names (if already known).
                If None, will extract from primary_code via regex.
        
        Returns:
            MiningResult with all discovered annotations.
        """
        # Step 1: Identify called subroutines
        if called_subroutines is None:
            called_subroutines = self._extract_calls(primary_code)
        
        # Normalize to uppercase
        called_subroutines = [s.upper() for s in called_subroutines]
        
        # Step 2: Resolve source files
        resolved = {}
        unresolved = []
        for sub_name in called_subroutines:
            if sub_name in self._source_index:
                resolved[sub_name] = str(self._source_index[sub_name])
            else:
                unresolved.append(sub_name)
        
        # Step 3: Extract annotations from each resolved source
        annotations = []
        for sub_name, file_path in resolved.items():
            ann = self._extract_annotations(sub_name, Path(file_path))
            if ann:
                annotations.append(ann)
        
        # Step 4: Use LLM to synthesize critical conventions
        critical_conventions = []
        storage_warnings = []
        
        if annotations:
            synthesis = self._synthesize_conventions(
                primary_filename, annotations
            )
            critical_conventions = synthesis.get("critical_conventions", [])
            storage_warnings = synthesis.get("storage_warnings", [])
        
        # Derive primary subroutine name
        match = re.search(
            r'^\s*SUBROUTINE\s+(\w+)', primary_code, re.MULTILINE | re.IGNORECASE
        )
        primary_name = match.group(1).upper() if match else primary_filename
        
        return MiningResult(
            primary_subroutine=primary_name,
            called_subroutines=called_subroutines,
            resolved_sources=resolved,
            unresolved=unresolved,
            annotations=annotations,
            critical_conventions=critical_conventions,
            storage_warnings=storage_warnings,
        )
    
    def _extract_calls(self, code: str) -> list[str]:
        """Extract all CALL targets from Fortran source."""
        calls = []
        for match in re.finditer(r'^\s*CALL\s+(\w+)', code, re.MULTILINE | re.IGNORECASE):
            name = match.group(1).upper()
            if name not in calls:
                calls.append(name)
        return calls
    
    def _extract_annotations(self, sub_name: str, file_path: Path) -> Optional[SourceAnnotation]:
        """Extract annotations from a single source file."""
        try:
            content = file_path.read_text(encoding='utf-8', errors='replace')
        except Exception:
            return None
        
        # Find the subroutine in the file
        pattern = rf'^\s*SUBROUTINE\s+{sub_name}\b'
        match = re.search(pattern, content, re.MULTILINE | re.IGNORECASE)
        
        if not match:
            # If not found as subroutine, use entire file
            sub_start = 0
        else:
            sub_start = match.start()
        
        # Extract header comments (comments above/below SUBROUTINE line)
        header_comments = self._extract_header_comments(content, sub_start)
        
        # Pattern matching for storage conventions
        storage_conventions = self._find_storage_patterns(content)
        
        # Pattern matching for algorithm references
        algorithm_references = self._find_algorithm_references(content)
        
        # Pattern matching for data format notes
        data_format_notes = self._find_data_format_notes(content)
        
        # Extract parameter documentation
        key_parameters = self._extract_parameter_docs(content, sub_start)
        
        # Find warnings/notes
        warnings = self._find_warnings(content)
        
        return SourceAnnotation(
            subroutine_name=sub_name,
            source_file=file_path.name,
            header_comments=header_comments,
            storage_conventions=storage_conventions,
            algorithm_references=algorithm_references,
            data_format_notes=data_format_notes,
            key_parameters=key_parameters,
            warnings=warnings,
        )
    
    def _extract_header_comments(self, content: str, sub_start: int) -> str:
        """Extract the documentation comment block near the subroutine."""
        lines = content[:sub_start + 2000].split('\n')
        
        # Find the SUBROUTINE line index
        sub_line_idx = None
        for i, line in enumerate(lines):
            if re.match(r'^\s*SUBROUTINE\s+', line, re.IGNORECASE):
                sub_line_idx = i
                break
        
        if sub_line_idx is None:
            sub_line_idx = 0
        
        # Collect comment lines before and after the SUBROUTINE line
        comment_lines = []
        
        # Comments before (up to 30 lines)
        start = max(0, sub_line_idx - 30)
        for i in range(start, sub_line_idx):
            line = lines[i]
            if line.startswith('C') or line.startswith('c') or line.startswith('*'):
                comment_lines.append(line[1:].rstrip())
        
        # Comments after (up to 40 lines after SUBROUTINE)
        for i in range(sub_line_idx + 1, min(len(lines), sub_line_idx + 40)):
            line = lines[i]
            if line.startswith('C') or line.startswith('c') or line.startswith('*'):
                comment_lines.append(line[1:].rstrip())
            elif line.strip() and not line.startswith('C') and not line.startswith('c'):
                # Stop at first non-comment, non-blank line
                break
        
        return '\n'.join(comment_lines)
    
    def _find_storage_patterns(self, content: str) -> list[str]:
        """Find storage convention keywords in source."""
        patterns = [
            (r'(?i)stored?\s+by\s+rows?', 'stored by rows'),
            (r'(?i)row[\s-]*major', 'row-major storage'),
            (r'(?i)column[\s-]*major', 'column-major storage'),
            (r'(?i)upper\s+triangle', 'upper triangle storage'),
            (r'(?i)lower\s+triangle', 'lower triangle storage'),
            (r'(?i)packed\s+(?:format|storage|form)', 'packed storage'),
            (r'(?i)band(?:ed)?\s+(?:matrix|storage|form)', 'banded storage'),
            (r'(?i)symmetric\s+(?:matrix|storage)', 'symmetric storage'),
            (r'(?i)row\s+order', 'row order'),
            (r'(?i)output.*(?:row|column)', 'output storage convention'),
        ]
        found = []
        for pattern, label in patterns:
            if re.search(pattern, content):
                # Get the actual matching line for more context
                match = re.search(pattern, content)
                # Get surrounding context
                start = max(0, match.start() - 60)
                end = min(len(content), match.end() + 60)
                context = content[start:end].replace('\n', ' ').strip()
                found.append(f"{label}: ...{context}...")
        return found
    
    def _find_algorithm_references(self, content: str) -> list[str]:
        """Find references to named algorithms or technical reports."""
        patterns = [
            r'(?i)FMMS[\s-]*\d+',  # NASA technical memos
            r'(?i)NASA[\s-]*(?:TN|TR|TM|CR)[\s-]*\d+',  # NASA reports
            r'(?i)(?:Timoshenko|Kirchhoff|Mindlin|Reissner)',  # Named theories
            r'(?i)(?:Gauss|Newton|Newmark|Wilson)',  # Named methods
            r'(?i)(?:AIAA|ASME)[\s-]*\d+',  # Conference refs
            r'(?i)ref(?:erence)?\.?\s*\d+',  # Numbered references
        ]
        found = []
        for pattern in patterns:
            for match in re.finditer(pattern, content):
                # Get context around the match
                start = max(0, match.start() - 40)
                end = min(len(content), match.end() + 40)
                context = content[start:end].replace('\n', ' ').strip()
                if context not in found:
                    found.append(context)
        return found
    
    def _find_data_format_notes(self, content: str) -> list[str]:
        """Find notes about data formats in comments."""
        notes = []
        for line in content.split('\n'):
            if not (line.startswith('C') or line.startswith('c') or line.startswith('*')):
                continue
            comment = line[1:].strip().lower()
            # Look for data format indicators
            keywords = ['matrix', 'array', 'vector', 'stored', 'format',
                       'input', 'output', 'dimension', 'layout']
            if any(kw in comment for kw in keywords):
                if len(line.strip()) > 10:  # Skip trivial comments
                    notes.append(line[1:].strip())
        # Deduplicate and limit
        return list(dict.fromkeys(notes))[:10]
    
    def _extract_parameter_docs(self, content: str, sub_start: int) -> dict[str, str]:
        """Extract parameter documentation from comments near the SUBROUTINE line."""
        params = {}
        lines = content[sub_start:sub_start + 3000].split('\n')
        
        for line in lines:
            if not (line.startswith('C') or line.startswith('c') or line.startswith('*')):
                continue
            comment = line[1:].strip()
            # Look for "VARNAME - description" or "VARNAME = description" patterns
            match = re.match(r'^(\w{1,8})\s*[-=:]\s*(.+)$', comment)
            if match:
                name = match.group(1).upper()
                desc = match.group(2).strip()
                if len(desc) > 5:  # Skip trivial entries
                    params[name] = desc
        
        return params
    
    def _find_warnings(self, content: str) -> list[str]:
        """Find WARNING/NOTE/CAUTION comments."""
        warnings = []
        for line in content.split('\n'):
            if not (line.startswith('C') or line.startswith('c') or line.startswith('*')):
                continue
            comment = line[1:].strip()
            if re.match(r'(?i)(NOTE|WARNING|CAUTION|IMPORTANT|N\.B\.)', comment):
                warnings.append(comment)
        return warnings
    
    def _synthesize_conventions(
        self, primary_filename: str, annotations: list[SourceAnnotation]
    ) -> dict:
        """Use LLM to synthesize critical conventions from all annotations."""
        # Build context from annotations
        context_parts = []
        for ann in annotations:
            part = f"=== {ann.subroutine_name} ({ann.source_file}) ===\n"
            if ann.header_comments:
                part += f"Comments:\n{ann.header_comments[:800]}\n"
            if ann.storage_conventions:
                part += f"Storage: {ann.storage_conventions}\n"
            if ann.algorithm_references:
                part += f"Algorithms: {ann.algorithm_references}\n"
            if ann.warnings:
                part += f"Warnings: {ann.warnings}\n"
            context_parts.append(part)
        
        annotations_text = '\n'.join(context_parts)
        
        user_prompt = f"""I am refactoring the NASTRAN-95 subroutine in file '{primary_filename}'.
It calls the following subroutines, whose documentation I have extracted:

{annotations_text}

Based on these annotations, what are the CRITICAL CONVENTIONS that the refactoring
MUST respect to produce numerically equivalent results?

Focus on:
1. Storage conventions that affect how matrices are interpreted
2. Implicit transpositions due to Fortran column-major vs routine expectations
3. Non-standard algorithms that shouldn't be replaced with textbook formulas
4. Any data format expectations that differ from modern defaults

Output as JSON with:
{{
    "critical_conventions": ["convention that MUST be respected", ...],
    "storage_warnings": ["specific warning about storage", ...],
    "algorithm_notes": ["note about algorithm", ...]
}}"""

        response = self.llm_client.generate(
            system_prompt=self.system_prompt,
            user_prompt=user_prompt,
            agent_name="source_mining"
        )
        
        if not response:
            return {"critical_conventions": [], "storage_warnings": []}
        
        from .json_utils import robust_json_parse
        result = robust_json_parse(response, {
            "critical_conventions": [],
            "storage_warnings": [],
            "algorithm_notes": [],
        })
        
        # Merge algorithm_notes into critical_conventions
        if result.get("algorithm_notes"):
            result["critical_conventions"].extend(result["algorithm_notes"])
        
        return result
