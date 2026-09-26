"""
Architecture Recovery Agent.

Analyzes multiple Fortran subroutines to recover:
- Call graph (which subroutines call which)
- Data flow (COMMON blocks, parameter passing)
- Module boundaries (logical groupings)
- Dependency ordering for refactoring
"""

import json
import re
from dataclasses import dataclass, field


@dataclass
class ModuleArchitecture:
    """Recovered architecture of a set of related subroutines."""
    
    # Call graph: caller -> [callees]
    call_graph: dict[str, list[str]] = field(default_factory=dict)
    
    # Data flow: COMMON block name -> [variables]
    common_blocks: dict[str, list[str]] = field(default_factory=dict)
    
    # Which subroutines use which COMMON blocks
    subroutine_commons: dict[str, list[str]] = field(default_factory=dict)
    
    # Logical module groupings
    modules: dict[str, list[str]] = field(default_factory=dict)
    
    # Dependency order (topological sort for refactoring)
    refactoring_order: list[str] = field(default_factory=list)
    
    # Data coupling between subroutines
    data_coupling: dict[str, dict[str, list[str]]] = field(default_factory=dict)
    
    def to_dict(self) -> dict:
        return {
            "call_graph": self.call_graph,
            "common_blocks": self.common_blocks,
            "subroutine_commons": self.subroutine_commons,
            "modules": self.modules,
            "refactoring_order": self.refactoring_order,
            "data_coupling": self.data_coupling,
        }


ARCHITECTURE_SYSTEM_PROMPT = """You are an Architecture Recovery Agent specialized in analyzing 
legacy Fortran 77 code structure from NASA NASTRAN-95.

Your task is to analyze a collection of related Fortran subroutines and recover:
1. The call graph (which subroutines call which others)
2. Data flow through COMMON blocks and parameters
3. Logical module boundaries (which subroutines form a cohesive unit)
4. The optimal order for refactoring (bottom-up: leaf subroutines first)

NASTRAN-95 architecture context:
- The system uses DMAP (Direct Matrix Abstraction Program) for high-level flow control
- Core computation is in 'functional modules' (subroutines in mis/ folder)
- Data is shared via COMMON blocks and passed through the GINO I/O system
- Element-level computations are relatively self-contained

Output your analysis as JSON with this structure:
{
    "call_graph": {"caller": ["callee1", "callee2"]},
    "common_blocks": {"block_name": ["var1", "var2"]},
    "subroutine_commons": {"sub_name": ["block1", "block2"]},
    "modules": {"module_name": ["sub1", "sub2"]},
    "refactoring_order": ["leaf_sub", "mid_sub", "top_sub"],
    "data_coupling": {
        "sub1": {"sub2": ["shared_var1", "shared_var2"]}
    }
}
"""


class ArchitectureAgent:
    """Agent that recovers software architecture from Fortran code."""
    
    def __init__(self, llm_client):
        self.llm_client = llm_client
        self.system_prompt = ARCHITECTURE_SYSTEM_PROMPT
    
    def recover_architecture(
        self, subroutines: dict[str, str], context: str = ""
    ) -> ModuleArchitecture:
        """
        Analyze multiple subroutines and recover their architecture.
        
        Args:
            subroutines: Dict mapping filename -> source code.
            context: Additional context about the module being analyzed.
            
        Returns:
            ModuleArchitecture with recovered structure.
        """
        # First pass: static analysis (no LLM needed)
        static_info = self._static_analysis(subroutines)
        
        # Second pass: LLM-based semantic analysis
        code_summary = self._summarize_code(subroutines)
        
        user_prompt = f"""Analyze the architecture of the following NASTRAN-95 subroutines.

Context: {context}

Static analysis results:
{json.dumps(static_info, indent=2)}

Code summaries:
{code_summary}

Provide the full architecture recovery as JSON."""

        response = self.llm_client.generate(
            system_prompt=self.system_prompt,
            user_prompt=user_prompt
        )
        
        arch_dict = self._parse_response(response)
        
        # Merge static analysis with LLM analysis
        if not arch_dict.get("call_graph"):
            arch_dict["call_graph"] = static_info.get("call_graph", {})
        if not arch_dict.get("common_blocks"):
            arch_dict["common_blocks"] = static_info.get("common_blocks", {})
        
        # Filter to only valid ModuleArchitecture fields
        valid_fields = set(ModuleArchitecture.__dataclass_fields__.keys())
        filtered = {k: v for k, v in arch_dict.items() if k in valid_fields}
        return ModuleArchitecture(**filtered)
    
    def _static_analysis(self, subroutines: dict[str, str]) -> dict:
        """Perform static analysis without LLM (regex-based)."""
        call_graph = {}
        common_blocks = {}
        subroutine_commons = {}
        
        # Known subroutine names
        sub_names = set()
        for filename, code in subroutines.items():
            for match in re.finditer(
                r'^\s*SUBROUTINE\s+(\w+)', code, re.MULTILINE | re.IGNORECASE
            ):
                sub_names.add(match.group(1).upper())
        
        for filename, code in subroutines.items():
            # Find subroutine name
            sub_match = re.search(
                r'^\s*SUBROUTINE\s+(\w+)', code, re.MULTILINE | re.IGNORECASE
            )
            if not sub_match:
                continue
            sub_name = sub_match.group(1).upper()
            
            # Find CALL statements
            calls = []
            for match in re.finditer(
                r'^\s*CALL\s+(\w+)', code, re.MULTILINE | re.IGNORECASE
            ):
                called = match.group(1).upper()
                if called not in calls:
                    calls.append(called)
            call_graph[sub_name] = calls
            
            # Find COMMON blocks
            commons = []
            for match in re.finditer(
                r'COMMON\s*/(\w+)/', code, re.IGNORECASE
            ):
                block_name = match.group(1).upper()
                if block_name not in commons:
                    commons.append(block_name)
                # Extract variables in this COMMON block
                # (simplified - full parsing would need continuation lines)
                line_match = re.search(
                    rf'COMMON\s*/{block_name}/\s*(.+?)$',
                    code, re.MULTILINE | re.IGNORECASE
                )
                if line_match and block_name not in common_blocks:
                    vars_str = line_match.group(1).strip()
                    variables = [v.strip() for v in vars_str.split(',') if v.strip()]
                    common_blocks[block_name] = variables
            
            subroutine_commons[sub_name] = commons
        
        return {
            "call_graph": call_graph,
            "common_blocks": common_blocks,
            "subroutine_commons": subroutine_commons,
        }
    
    def _summarize_code(self, subroutines: dict[str, str]) -> str:
        """Create a brief summary of each subroutine for the LLM."""
        summaries = []
        for filename, code in subroutines.items():
            # Take first 50 lines as summary
            lines = code.split('\n')[:50]
            summaries.append(f"--- {filename} (first 50 lines) ---\n" + '\n'.join(lines))
        return '\n\n'.join(summaries)
    
    def _parse_response(self, response: str) -> dict:
        """Extract JSON from LLM response."""
        from .json_utils import robust_json_parse
        
        defaults = {
            "call_graph": {},
            "common_blocks": {},
            "subroutine_commons": {},
            "modules": {},
            "refactoring_order": [],
            "data_coupling": {},
        }
        
        parsed = robust_json_parse(response, defaults)
        for key, default_val in defaults.items():
            if key not in parsed:
                parsed[key] = default_val
        return parsed
