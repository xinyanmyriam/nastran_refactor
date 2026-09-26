"""
Modern Design Agent.

Takes the physics specification and architecture information to produce
a modern Python design that:
- Preserves all physical/mathematical functionality
- Uses modern OOP patterns (classes, inheritance, composition)
- Follows NumPy/SciPy idioms for numerical computing
- Is readable, testable, and maintainable
"""

import json
from dataclasses import dataclass, field


@dataclass
class ModernDesign:
    """Design document for the modern Python implementation."""
    
    # Module structure
    module_name: str = ""
    class_hierarchy: dict[str, list[str]] = field(default_factory=dict)  # class -> [methods]
    
    # Class design
    classes: list[dict] = field(default_factory=list)
    # Each class: {"name", "base_class", "attributes", "methods", "docstring"}
    
    # Function signatures
    functions: list[dict] = field(default_factory=list)
    # Each function: {"name", "args", "returns", "docstring"}
    
    # Data structures
    data_structures: list[dict] = field(default_factory=list)
    # Replacement for COMMON blocks
    
    # Design decisions
    design_decisions: list[str] = field(default_factory=list)
    
    # Mapping from old to new
    fortran_to_python_map: dict[str, str] = field(default_factory=dict)
    
    # Dependencies
    python_dependencies: list[str] = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "module_name": self.module_name,
            "class_hierarchy": self.class_hierarchy,
            "classes": self.classes,
            "functions": self.functions,
            "data_structures": self.data_structures,
            "design_decisions": self.design_decisions,
            "fortran_to_python_map": self.fortran_to_python_map,
            "python_dependencies": self.python_dependencies,
        }


MODERN_DESIGN_SYSTEM_PROMPT = """You are a Modern Design Agent that transforms legacy Fortran 77 
finite element code designs into modern Python architectures.

Your task is to take a physics specification (extracted from NASTRAN-95 code) and produce
a modern Python design that:

1. PRESERVES all physical/mathematical functionality exactly
2. Uses modern Python patterns:
   - NumPy arrays instead of Fortran arrays
   - Classes instead of COMMON blocks
   - Type hints and docstrings
   - Clear variable names (not 6-char Fortran limits)
3. Is TESTABLE - each component can be independently verified
4. Is COMPOSABLE - elements can be assembled into larger systems

Design principles:
- Element classes should inherit from a base Element class
- Stiffness matrix computation should be a method returning numpy array
- Material properties should be encapsulated in Material classes
- Coordinate transformations should be explicit and reusable
- All physical quantities should have clear units in documentation

Fortran → Python mapping patterns:
- COMMON block → class attributes or dataclass
- SUBROUTINE with output args → method returning values
- EQUIVALENCE → separate variables (memory is cheap now)
- Fixed-format I/O → structured data (dict, dataclass)
- GOTO → structured control flow (if/else, loops)
- Implicit typing (I-N integer) → explicit type hints

Output as JSON:
{
    "module_name": "element_rod",
    "class_hierarchy": {"ClassName": ["method1", "method2"]},
    "classes": [
        {
            "name": "RodElement",
            "base_class": "Element",
            "attributes": [
                {"name": "node_ids", "type": "np.ndarray", "description": "..."},
                ...
            ],
            "methods": [
                {"name": "stiffness_matrix", "args": [], "returns": "np.ndarray", 
                 "description": "Compute 2x2 element stiffness matrix"},
                ...
            ],
            "docstring": "CROD element implementation..."
        }
    ],
    "functions": [...],
    "data_structures": [...],
    "design_decisions": ["Decision 1: ...", "Decision 2: ..."],
    "fortran_to_python_map": {"FORTRAN_VAR": "python_var"},
    "python_dependencies": ["numpy", "scipy"]
}
"""


class ModernDesignAgent:
    """Agent that designs modern Python architecture from physics specs."""
    
    def __init__(self, llm_client):
        self.llm_client = llm_client
        self.system_prompt = MODERN_DESIGN_SYSTEM_PROMPT
    
    def design(
        self,
        physics_spec: dict,
        architecture: dict = None,
        original_code: str = ""
    ) -> ModernDesign:
        """
        Produce a modern Python design from physics specification.
        
        Args:
            physics_spec: Dict from PhysicsSpecAgent.
            architecture: Dict from ArchitectureAgent (optional).
            original_code: Original Fortran code for reference (optional).
            
        Returns:
            ModernDesign object.
        """
        user_prompt = f"""Design a modern Python implementation for the following 
finite element component.

Physics Specification:
{json.dumps(physics_spec, indent=2)}

"""
        if architecture:
            user_prompt += f"""Architecture Context:
{json.dumps(architecture, indent=2)}

"""
        if original_code:
            # Include first 100 lines of original for reference
            lines = original_code.split('\n')[:100]
            user_prompt += f"""Original Fortran Code (first 100 lines, for reference):
```fortran
{chr(10).join(lines)}
```

"""
        user_prompt += "Produce the modern Python design as JSON."

        response = self.llm_client.generate(
            system_prompt=self.system_prompt,
            user_prompt=user_prompt
        )
        
        design_dict = self._parse_response(response)
        # Filter to only valid ModernDesign fields
        valid_fields = set(ModernDesign.__dataclass_fields__.keys())
        filtered = {k: v for k, v in design_dict.items() if k in valid_fields}
        return ModernDesign(**filtered)
    
    def _parse_response(self, response: str) -> dict:
        """Extract JSON from LLM response."""
        from .json_utils import robust_json_parse
        
        defaults = ModernDesign().to_dict()
        parsed = robust_json_parse(response, defaults)
        for key, default_val in defaults.items():
            if key not in parsed:
                parsed[key] = default_val
        return parsed
