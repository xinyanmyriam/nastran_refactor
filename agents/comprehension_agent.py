"""
Code Comprehension Agent.

Reads Fortran 77 source code and extracts structured understanding of:
- Physical/mathematical intent
- Algorithm description
- Variable semantics
- Input/output specification
"""

import json
from dataclasses import dataclass, field
from typing import Optional

# Reuse the LLM client from the main project
import sys
sys.path.insert(0, str(__file__).rsplit("nastran_refactor", 1)[0])


@dataclass
class CodeIntent:
    """Structured representation of a subroutine's intent."""
    subroutine_name: str
    source_file: str
    
    # Physical intent
    physical_domain: str  # e.g., "structural mechanics", "heat transfer"
    physical_operation: str  # e.g., "stiffness matrix computation"
    governing_equation: str  # e.g., "F = EA/L * (u2 - u1)"
    
    # Algorithm description
    algorithm_steps: list[str] = field(default_factory=list)
    numerical_method: str = ""  # e.g., "direct stiffness method"
    
    # Variable semantics
    input_variables: dict[str, str] = field(default_factory=dict)  # name -> meaning
    output_variables: dict[str, str] = field(default_factory=dict)
    local_variables: dict[str, str] = field(default_factory=dict)
    
    # Dependencies
    called_subroutines: list[str] = field(default_factory=list)
    common_blocks_used: list[str] = field(default_factory=list)
    
    # Assumptions and constraints
    assumptions: list[str] = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "subroutine_name": self.subroutine_name,
            "source_file": self.source_file,
            "physical_domain": self.physical_domain,
            "physical_operation": self.physical_operation,
            "governing_equation": self.governing_equation,
            "algorithm_steps": self.algorithm_steps,
            "numerical_method": self.numerical_method,
            "input_variables": self.input_variables,
            "output_variables": self.output_variables,
            "local_variables": self.local_variables,
            "called_subroutines": self.called_subroutines,
            "common_blocks_used": self.common_blocks_used,
            "assumptions": self.assumptions,
        }


COMPREHENSION_SYSTEM_PROMPT = """You are a Code Comprehension Agent specialized in understanding 
legacy Fortran 77 finite element analysis code, specifically from NASA NASTRAN-95.

Your task is to read Fortran 77 source code and produce a structured understanding of:
1. What physical/mathematical operation the code performs
2. What algorithm it implements step by step
3. What each variable represents physically
4. What assumptions are made

Context about NASTRAN-95:
- Written in Fortran 77 with fixed-format (columns 1-72)
- Uses COMMON blocks for global data sharing
- Variable names are limited to 6 characters
- Uses EQUIVALENCE for memory reuse
- Part of a finite element analysis system for structural/thermal analysis
- The 'mis/' folder contains core computational subroutines

When analyzing code, pay special attention to:
- Matrix operations (stiffness, mass, load vectors)
- Element formulations (shape functions, Jacobians)
- Coordinate transformations
- Boundary condition handling
- Material property lookups

Output your analysis as a JSON object with the following structure:
{
    "subroutine_name": "...",
    "physical_domain": "structural mechanics | heat transfer | ...",
    "physical_operation": "what this code does physically",
    "governing_equation": "the mathematical equation being implemented",
    "algorithm_steps": ["step 1", "step 2", ...],
    "numerical_method": "method name",
    "input_variables": {"var_name": "physical meaning"},
    "output_variables": {"var_name": "physical meaning"},
    "local_variables": {"var_name": "physical meaning"},
    "called_subroutines": ["sub1", "sub2"],
    "common_blocks_used": ["block1", "block2"],
    "assumptions": ["assumption 1", "assumption 2"]
}
"""


class ComprehensionAgent:
    """Agent that understands Fortran 77 legacy code."""
    
    def __init__(self, llm_client):
        """
        Args:
            llm_client: LLM client with a `generate(system_prompt, user_prompt)` method.
        """
        self.llm_client = llm_client
        self.system_prompt = COMPREHENSION_SYSTEM_PROMPT
    
    def analyze_subroutine(self, fortran_code: str, source_file: str = "") -> CodeIntent:
        """
        Analyze a single Fortran subroutine and extract its intent.
        
        Args:
            fortran_code: The Fortran 77 source code text.
            source_file: Name of the source file (for reference).
            
        Returns:
            CodeIntent object with structured understanding.
        """
        user_prompt = f"""Analyze the following Fortran 77 subroutine from NASTRAN-95 
(file: {source_file}):

```fortran
{fortran_code}
```

Provide a structured analysis as JSON."""

        response = self.llm_client.generate(
            system_prompt=self.system_prompt,
            user_prompt=user_prompt
        )
        
        # Parse JSON from response
        intent_dict = self._parse_response(response)
        intent_dict["source_file"] = source_file
        
        # Filter to only valid CodeIntent fields
        valid_fields = set(CodeIntent.__dataclass_fields__.keys())
        filtered = {k: v for k, v in intent_dict.items() if k in valid_fields}
        return CodeIntent(**filtered)
    
    def analyze_multiple(self, subroutines: dict[str, str]) -> list[CodeIntent]:
        """
        Analyze multiple related subroutines.
        
        Args:
            subroutines: Dict mapping filename -> source code.
            
        Returns:
            List of CodeIntent objects.
        """
        results = []
        for filename, code in subroutines.items():
            intent = self.analyze_subroutine(code, filename)
            results.append(intent)
        return results
    
    def _parse_response(self, response: str) -> dict:
        """Extract JSON from LLM response, with robust error handling."""
        from .json_utils import robust_json_parse
        
        defaults = {
            "subroutine_name": "unknown",
            "physical_domain": "unknown",
            "physical_operation": "unknown",
            "governing_equation": "",
            "algorithm_steps": [],
            "numerical_method": "",
            "input_variables": {},
            "output_variables": {},
            "local_variables": {},
            "called_subroutines": [],
            "common_blocks_used": [],
            "assumptions": [],
        }
        
        parsed = robust_json_parse(response, defaults)
        
        # Ensure all required fields exist
        for key, default_val in defaults.items():
            if key not in parsed:
                parsed[key] = default_val
        
        return parsed
