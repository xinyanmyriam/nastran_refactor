"""
Equivalence Verification Agent.

Verifies that the refactored Python code produces results equivalent to
the original NASTRAN Fortran code. This is the key innovation for the
refactoring scenario (vs. the original paper's physics-aware verification).

Verification levels:
1. Structural equivalence: same inputs/outputs, same matrix dimensions
2. Numerical equivalence: same results within tolerance
3. Behavioral equivalence: same behavior across a test suite
"""

import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np


@dataclass
class EquivalenceReport:
    """Report on equivalence between original and refactored code."""
    
    # Overall result
    equivalent: bool = False
    
    # Structural checks
    structural_checks: dict = field(default_factory=dict)
    # e.g., {"matrix_dimensions": True, "output_count": True}
    
    # Numerical checks
    numerical_checks: dict = field(default_factory=dict)
    # e.g., {"max_error": 1e-15, "l2_error": 1e-16, "within_tolerance": True}
    
    # Behavioral checks (across test suite)
    test_cases_passed: int = 0
    test_cases_total: int = 0
    failed_cases: list[dict] = field(default_factory=list)
    
    # Detailed comparison
    comparison_details: list[dict] = field(default_factory=list)
    
    def to_dict(self) -> dict:
        return {
            "equivalent": self.equivalent,
            "structural_checks": self.structural_checks,
            "numerical_checks": self.numerical_checks,
            "test_cases_passed": self.test_cases_passed,
            "test_cases_total": self.test_cases_total,
            "failed_cases": self.failed_cases,
            "comparison_details": self.comparison_details,
        }


EQUIVALENCE_SYSTEM_PROMPT = """You are an Equivalence Verification Agent that checks whether
refactored Python code produces results identical to the original NASTRAN Fortran code.

Your verification approach:
1. STRUCTURAL: Check that the Python code has the same interface (inputs/outputs)
2. NUMERICAL: Compare computed values (stiffness matrices, displacements, stresses)
3. BEHAVIORAL: Run multiple test cases and compare results

For finite element code, key equivalence checks include:
- Element stiffness matrix: K_python ≈ K_fortran (within machine precision)
- Stress recovery: σ_python ≈ σ_fortran
- Boundary condition handling: same DOF elimination
- Assembly: global matrix matches

Tolerance guidelines:
- For direct computation (no iteration): expect machine precision (< 1e-12)
- For iterative solvers: expect solver tolerance (< 1e-8)
- For different algorithms (same physics): expect engineering accuracy (< 1e-6)

Generate test cases that cover:
- Simple cases with analytical solutions
- Edge cases (zero length, zero area, etc.)
- Typical engineering values
- Extreme values (very stiff, very flexible)
"""


class EquivalenceAgent:
    """Agent that verifies equivalence between original and refactored code."""
    
    def __init__(self, llm_client=None, tolerance: float = 1e-6):
        """
        Args:
            llm_client: Optional LLM client for generating test cases.
            tolerance: Relative numerical tolerance for equivalence.
        """
        self.llm_client = llm_client
        self.tolerance = tolerance
    
    def verify_element_stiffness(
        self,
        python_code: str,
        reference_results: list[dict],
    ) -> EquivalenceReport:
        """
        Verify that Python element code produces correct stiffness matrices.
        
        Args:
            python_code: The refactored Python code.
            reference_results: List of test cases with expected results.
                Each dict: {"inputs": {...}, "expected_K": [[...]], "description": "..."}
                
        Returns:
            EquivalenceReport with detailed comparison.
        """
        report = EquivalenceReport()
        report.test_cases_total = len(reference_results)
        
        for i, test_case in enumerate(reference_results):
            result = self._run_single_test(python_code, test_case)
            
            if result["passed"]:
                report.test_cases_passed += 1
            else:
                report.failed_cases.append({
                    "index": i,
                    "description": test_case.get("description", f"Test {i}"),
                    "error": result.get("error", "Unknown"),
                    "max_diff": result.get("max_diff", None),
                })
            
            report.comparison_details.append(result)
        
        # Overall equivalence
        report.equivalent = (report.test_cases_passed == report.test_cases_total)
        
        # Aggregate numerical checks
        if report.comparison_details:
            max_errors = [
                d.get("max_diff", 0) for d in report.comparison_details 
                if d.get("max_diff") is not None
            ]
            if max_errors:
                report.numerical_checks = {
                    "max_error_across_tests": max(max_errors),
                    "mean_error_across_tests": sum(max_errors) / len(max_errors),
                    "tolerance": self.tolerance,
                    "within_tolerance": max(max_errors) < self.tolerance,
                }
        
        return report
    
    def verify_with_nastran_output(
        self,
        python_code: str,
        nastran_input_file: str,
        nastran_output_file: str,
    ) -> EquivalenceReport:
        """
        Verify Python code against actual NASTRAN input/output files.
        
        Args:
            python_code: The refactored Python code.
            nastran_input_file: Path to NASTRAN .inp file.
            nastran_output_file: Path to NASTRAN .out file with expected results.
            
        Returns:
            EquivalenceReport.
        """
        # Parse NASTRAN input to extract problem definition
        problem = self._parse_nastran_input(nastran_input_file)
        
        # Parse NASTRAN output to extract reference results
        reference = self._parse_nastran_output(nastran_output_file)
        
        # Run Python code with same inputs
        python_result = self._run_python_solver(python_code, problem)
        
        # Compare
        report = self._compare_results(python_result, reference)
        return report
    
    def generate_test_suite(
        self, physics_spec: dict, num_cases: int = 10
    ) -> list[dict]:
        """
        Use LLM to generate a comprehensive test suite for equivalence testing.
        
        Args:
            physics_spec: Physics specification of the element.
            num_cases: Number of test cases to generate.
            
        Returns:
            List of test case dicts.
        """
        if not self.llm_client:
            return self._default_test_suite(physics_spec)
        
        prompt = f"""Generate {num_cases} test cases for verifying a 
{physics_spec.get('element_type', 'unknown')} finite element implementation.

Physics specification:
{json.dumps(physics_spec, indent=2)}

Each test case should include:
- Input parameters (geometry, material properties, loads)
- Expected stiffness matrix (computed analytically)
- Expected displacements (if applicable)
- Description of what the test verifies

Output as JSON array of test cases."""

        response = self.llm_client.generate(
            system_prompt=EQUIVALENCE_SYSTEM_PROMPT,
            user_prompt=prompt
        )
        
        try:
            if "```json" in response:
                json_str = response.split("```json")[1].split("```")[0].strip()
            elif "```" in response:
                json_str = response.split("```")[1].split("```")[0].strip()
            else:
                json_str = response.strip()
            return json.loads(json_str)
        except (json.JSONDecodeError, IndexError):
            return self._default_test_suite(physics_spec)
    
    def _run_single_test(self, python_code: str, test_case: dict) -> dict:
        """Run a single equivalence test using separate files to avoid conflicts."""
        import tempfile
        import os
        
        # Write generated code to a temp module file
        temp_dir = tempfile.mkdtemp()
        code_file = Path(temp_dir) / "generated_element.py"
        code_file.write_text(python_code, encoding='utf-8')
        
        # Write test harness as a separate script
        harness_file = Path(temp_dir) / "test_harness.py"
        harness_code = self._build_test_harness(test_case)
        harness_file.write_text(harness_code, encoding='utf-8')
        
        try:
            proc = subprocess.run(
                [sys.executable, str(harness_file)],
                capture_output=True, text=True, timeout=30,
                cwd=temp_dir
            )
            
            if proc.returncode != 0:
                return {"passed": False, "error": proc.stderr[:1000]}
            
            if "__RESULT__" in proc.stdout:
                result_str = proc.stdout.split("__RESULT__")[1].strip()
                return json.loads(result_str)
            else:
                return {"passed": False, "error": "No result output", "stdout": proc.stdout[:500]}
                
        except subprocess.TimeoutExpired:
            return {"passed": False, "error": "Timeout"}
        except Exception as e:
            return {"passed": False, "error": str(e)}
        finally:
            # Cleanup
            import shutil
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    def _build_test_harness(self, test_case: dict) -> str:
        """Build a test harness script that imports and tests the generated code."""
        inputs_json = json.dumps(test_case.get("inputs", {}))
        expected_K_json = json.dumps(test_case.get("expected_K", []))
        tolerance = self.tolerance
        
        return f'''import numpy as np
import json
import sys
import importlib.util

# Load the generated element module
spec = importlib.util.spec_from_file_location("generated_element", "generated_element.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

# Test case parameters
inputs = {inputs_json}
expected_K = np.array({expected_K_json})

L = inputs.get("L", 1.0)
E = inputs.get("E", 1.0)
A = inputs.get("A", 1.0)
G = E / (2 * (1 + 0.3))
theta = inputs.get("theta", 0.0)
phi = inputs.get("phi", 0.0)

cx = np.cos(phi) * np.cos(theta)
cy = np.cos(phi) * np.sin(theta)
cz = np.sin(phi)
coord_a = np.array([0.0, 0.0, 0.0])
coord_b = L * np.array([cx, cy, cz])
node_coords = np.array([coord_a, coord_b])

def try_create():
    """Try multiple interface patterns to create the element."""
    # Pattern 1: create_element factory
    if hasattr(mod, "create_element"):
        try:
            return mod.create_element(**inputs)
        except Exception:
            pass

    # Find element class (any class with stiffness_matrix method)
    elem_cls = None
    for name in dir(mod):
        obj = getattr(mod, name)
        if isinstance(obj, type) and hasattr(obj, "stiffness_matrix") and name != "Element":
            elem_cls = obj
            break
    
    if elem_cls is None:
        raise RuntimeError("No element class found with stiffness_matrix method")
    
    # Get Material class if exists
    mat_cls = getattr(mod, "Material", None)
    prop_cls = getattr(mod, "ElementProperties", None)
    
    # Try various constructor patterns
    attempts = []
    
    # Pattern A: (node_ids, coord_a, coord_b, material, properties, cs_ids)
    if mat_cls and prop_cls:
        try:
            mat = mat_cls(youngs_modulus=E, shear_modulus=G)
            prop = prop_cls(area=A, polar_moment=0.0)
            return elem_cls(node_ids=(1,2), coord_a=coord_a, coord_b=coord_b, material=mat, properties=prop)
        except Exception as e:
            attempts.append(f"A1: {{e}}")
        try:
            mat = mat_cls(youngs_modulus=E, shear_modulus=G)
            prop = prop_cls(area=A, polar_moment=0.0)
            return elem_cls((1,2), coord_a, coord_b, mat, prop)
        except Exception as e:
            attempts.append(f"A2: {{e}}")

    # Pattern B: (node_ids, material, area, polar_moment, node_coords) - Qwen3 style
    if mat_cls:
        try:
            mat = mat_cls(E, G)
            return elem_cls(np.array([1,2]), mat, A, 0.0, node_coords)
        except Exception as e:
            attempts.append(f"B1: {{e}}")
        try:
            mat = mat_cls(E=E, G=G)
            return elem_cls(node_ids=np.array([1,2]), material=mat, area=A, polar_moment=0.0, node_coords=node_coords)
        except Exception as e:
            attempts.append(f"B2: {{e}}")

    # Pattern C: (x1,y1,z1,x2,y2,z2,area,J,E,G) - direct coords style
    try:
        return elem_cls(coord_a[0],coord_a[1],coord_a[2], coord_b[0],coord_b[1],coord_b[2], A, 0.0, E, G)
    except Exception as e:
        attempts.append(f"C1: {{e}}")
    try:
        return elem_cls(x1=coord_a[0],y1=coord_a[1],z1=coord_a[2], x2=coord_b[0],y2=coord_b[1],z2=coord_b[2], area=A, torsion_constant=0.0, young_modulus=E, shear_modulus=G)
    except Exception as e:
        attempts.append(f"C2: {{e}}")

    # Pattern D: (node_coords, material) simple
    if mat_cls:
        try:
            mat = mat_cls(E=E, A=A)
            return elem_cls(node_coords, mat)
        except Exception as e:
            attempts.append(f"D1: {{e}}")

    # Pattern E: (node_coords, E, A)
    try:
        return elem_cls(node_coords=node_coords, E=E, A=A)
    except Exception as e:
        attempts.append(f"E1: {{e}}")

    raise RuntimeError(f"All patterns failed: {{'; '.join(attempts[:6])}}")

try:
    element = try_create()
    K_computed = element.stiffness_matrix()
    
    # Handle 12x12 -> 6x6 extraction
    if K_computed.shape == (12, 12) and expected_K.shape == (6, 6):
        idx = [0, 1, 2, 6, 7, 8]
        K_computed = K_computed[np.ix_(idx, idx)]
    
    if K_computed.shape != expected_K.shape:
        result = {{"passed": False, "error": f"Shape mismatch: {{K_computed.shape}} vs {{expected_K.shape}}", "K_shape": list(K_computed.shape)}}
    else:
        max_diff = float(np.max(np.abs(K_computed - expected_K)))
        max_val = max(float(np.max(np.abs(expected_K))), 1e-15)
        rel_diff = max_diff / max_val
        is_symmetric = bool(np.allclose(K_computed, K_computed.T, atol=1e-12))
        result = {{"passed": rel_diff < {tolerance} or max_diff < 1e-10, "max_diff": max_diff, "rel_diff": rel_diff, "is_symmetric": is_symmetric, "K_shape": list(K_computed.shape)}}
except Exception as e:
    result = {{"passed": False, "error": str(e)}}

print("__RESULT__")
print(json.dumps(result, default=str))
'''

    def _default_test_suite(self, physics_spec: dict) -> list[dict]:
        """Generate default test cases based on element type."""
        element_type = str(physics_spec.get("element_type", "rod")).lower()
        
        if "rod" in element_type:
            return [
                {
                    "description": "Unit rod element along x-axis",
                    "inputs": {"L": 1.0, "E": 1.0, "A": 1.0, "theta": 0.0, "phi": 0.0},
                    "expected_K": [
                        [1, 0, 0, -1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [-1, 0, 0, 1, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0],
                    ],
                },
            ]
        else:
            return []
    
    def _parse_nastran_input(self, filepath: str) -> dict:
        problem = {"nodes": [], "elements": [], "materials": [], "loads": [], "constraints": []}
        return problem
    
    def _parse_nastran_output(self, filepath: str) -> dict:
        reference = {"displacements": {}, "stresses": {}, "forces": {}}
        return reference
    
    def _run_python_solver(self, python_code: str, problem: dict) -> dict:
        return {}
    
    def _compare_results(self, computed: dict, reference: dict) -> EquivalenceReport:
        report = EquivalenceReport()
        return report
