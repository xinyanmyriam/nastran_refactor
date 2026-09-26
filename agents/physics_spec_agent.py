"""
Physics Specification Agent.

Takes the code intent from the Comprehension Agent and produces a formal,
structured physics specification that can drive modern code generation.

This is analogous to the Physics Model Agent in the original paper,
but its input comes from legacy code understanding rather than user requirements.
"""

import json
from dataclasses import dataclass, field


@dataclass
class PhysicsSpecification:
    """Formal physics specification extracted from legacy code."""
    
    # Problem classification
    domain: str  # "structural_mechanics", "heat_transfer", etc.
    element_type: str  # "rod", "bar", "tria3", "quad4", etc.
    analysis_type: str  # "static", "dynamic", "buckling", etc.
    dimensionality: str  # "1D", "2D", "3D"
    
    # Governing equations
    governing_equation: str  # LaTeX or text form
    weak_form: str  # Weak/variational form if applicable
    
    # Element formulation
    num_nodes: int = 0
    dof_per_node: int = 0
    total_dof: int = 0
    shape_functions: str = ""  # Description or formula
    
    # Material model
    constitutive_law: str = ""  # e.g., "linear elastic isotropic"
    material_parameters: list[str] = field(default_factory=list)  # ["E", "nu", "rho"]
    
    # Stiffness matrix
    stiffness_matrix_formula: str = ""  # Symbolic or algorithmic description
    integration_scheme: str = ""  # "exact", "1-point Gauss", "2x2 Gauss"
    
    # Load vector
    load_types: list[str] = field(default_factory=list)  # ["concentrated", "distributed"]
    
    # Coordinate system
    local_coordinates: str = ""  # Description of local element coordinate system
    transformation: str = ""  # How to transform to global
    
    # Stress/strain recovery
    stress_recovery_method: str = ""
    output_quantities: list[str] = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "domain": self.domain,
            "element_type": self.element_type,
            "analysis_type": self.analysis_type,
            "dimensionality": self.dimensionality,
            "governing_equation": self.governing_equation,
            "weak_form": self.weak_form,
            "num_nodes": self.num_nodes,
            "dof_per_node": self.dof_per_node,
            "total_dof": self.total_dof,
            "shape_functions": self.shape_functions,
            "constitutive_law": self.constitutive_law,
            "material_parameters": self.material_parameters,
            "stiffness_matrix_formula": self.stiffness_matrix_formula,
            "integration_scheme": self.integration_scheme,
            "load_types": self.load_types,
            "local_coordinates": self.local_coordinates,
            "transformation": self.transformation,
            "stress_recovery_method": self.stress_recovery_method,
            "output_quantities": self.output_quantities,
        }


PHYSICS_SPEC_SYSTEM_PROMPT = """You are a Physics Specification Agent for finite element analysis.

Your task is to take a code intent description (extracted from legacy Fortran code) and produce
a formal, complete physics specification that captures ALL the physical and mathematical 
content of the original code.

This specification will be used to generate modern Python code that is behaviorally equivalent
to the original Fortran implementation.

Requirements for the specification:
1. Must be COMPLETE - capture every physical assumption and mathematical operation
2. Must be PRECISE - use standard notation and unambiguous descriptions
3. Must be VERIFIABLE - include enough detail to write equivalence tests
4. Must preserve any numerical tricks or special cases from the original code

Output as JSON with this structure:
{
    "domain": "structural_mechanics | heat_transfer | ...",
    "element_type": "rod | bar | tria3 | quad4 | ...",
    "analysis_type": "static | dynamic | ...",
    "dimensionality": "1D | 2D | 3D",
    "governing_equation": "equation in text/LaTeX",
    "weak_form": "variational form",
    "num_nodes": integer,
    "dof_per_node": integer,
    "total_dof": integer,
    "shape_functions": "description or formula",
    "constitutive_law": "material model description",
    "material_parameters": ["param1", "param2"],
    "stiffness_matrix_formula": "K = ...",
    "integration_scheme": "method description",
    "load_types": ["type1", "type2"],
    "local_coordinates": "coordinate system description",
    "transformation": "transformation method",
    "stress_recovery_method": "method description",
    "output_quantities": ["displacement", "stress", ...]
}
"""


class PhysicsSpecAgent:
    """Agent that produces formal physics specifications from code intent."""
    
    def __init__(self, llm_client):
        self.llm_client = llm_client
        self.system_prompt = PHYSICS_SPEC_SYSTEM_PROMPT
    
    def generate_specification(self, code_intent: dict, context: str = "") -> PhysicsSpecification:
        """
        Generate a physics specification from code intent.
        
        Args:
            code_intent: Dict from ComprehensionAgent (CodeIntent.to_dict()).
            context: Additional context (e.g., NASTRAN documentation excerpts).
            
        Returns:
            PhysicsSpecification object.
        """
        user_prompt = f"""Based on the following code intent extracted from NASTRAN-95 legacy code,
produce a complete physics specification.

Code Intent:
{json.dumps(code_intent, indent=2)}

Additional Context:
{context if context else "None"}

Produce the physics specification as JSON."""

        response = self.llm_client.generate(
            system_prompt=self.system_prompt,
            user_prompt=user_prompt
        )
        
        spec_dict = self._parse_response(response)
        # Filter to only valid PhysicsSpecification fields
        valid_fields = set(PhysicsSpecification.__dataclass_fields__.keys())
        filtered = {k: v for k, v in spec_dict.items() if k in valid_fields}
        return PhysicsSpecification(**filtered)
    
    def _parse_response(self, response: str) -> dict:
        """Extract JSON from LLM response."""
        from .json_utils import robust_json_parse
        
        defaults = PhysicsSpecification(
            domain="unknown", element_type="unknown",
            analysis_type="unknown", dimensionality="unknown",
            governing_equation="", weak_form=""
        ).to_dict()
        
        parsed = robust_json_parse(response, defaults)
        for key, default_val in defaults.items():
            if key not in parsed:
                parsed[key] = default_val
        return parsed
