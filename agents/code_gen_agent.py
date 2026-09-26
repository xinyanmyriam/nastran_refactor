"""
Code Generation Agent.

Generates modern Python code from the design specification.
Includes iterative repair mechanism (same as the original paper).
"""

import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class GenerationResult:
    """Result of code generation attempt."""
    success: bool
    code: str
    iteration: int
    error_message: str = ""
    metrics: dict = field(default_factory=dict)


CODE_GEN_SYSTEM_PROMPT = """You are a Code Generation Agent that produces Python implementations
of finite element components based on a design specification.

Requirements for generated code:
1. SELF-CONTAINED: The code must run with only numpy and scipy as dependencies
2. CORRECT: Must implement the exact physics and numerics from the specification
3. COMPLETE: Must include all methods and classes from the design
4. DOCUMENTED: Must have docstrings with physical meaning of parameters
5. TESTABLE: Must include a `verify()` function that runs a basic self-test

Code structure requirements:
- Use numpy for all array operations
- Use scipy.sparse for large matrices (optional for small elements)
- Include type hints
- Include a main block that demonstrates usage
- The verify() function should:
  - Create a simple test case with known analytical solution
  - Compute the element stiffness matrix
  - Verify dimensions and symmetry
  - Return a dict with test results

Example output structure:
```python
import numpy as np
from typing import Optional

class RodElement:
    \"\"\"CROD element - axial rod with tension/compression only.\"\"\"
    
    def __init__(self, node_coords, E, A):
        ...
    
    def stiffness_matrix(self) -> np.ndarray:
        ...
    
    def stress(self, displacements: np.ndarray) -> float:
        ...

def verify() -> dict:
    \"\"\"Self-verification with known analytical solution.\"\"\"
    ...
    return {"passed": True, "details": {...}}

if __name__ == "__main__":
    results = verify()
    print(json.dumps(results, indent=2))
```

IMPORTANT: Generate COMPLETE, RUNNABLE code. Do not use placeholders or '...'.
"""

ERROR_FEEDBACK_TEMPLATE = """The previous code generation attempt failed.

Previous code:
```python
{code}
```

Error:
{error}

Please fix the error and regenerate the complete code. 
Ensure the code is self-contained and runnable with only numpy/scipy."""

METRIC_FEEDBACK_TEMPLATE = """The previous code generation attempt produced incorrect results.

Previous code:
```python
{code}
```

Test results:
{metrics}

Expected behavior:
{expected}

Please fix the numerical/physical errors and regenerate the complete code."""


class CodeGenAgent:
    """Agent that generates Python code with iterative repair."""
    
    def __init__(self, llm_client, max_iterations: int = 5, timeout: int = 60):
        self.llm_client = llm_client
        self.max_iterations = max_iterations
        self.timeout = timeout
        self.system_prompt = CODE_GEN_SYSTEM_PROMPT
    
    def generate(
        self,
        design: dict,
        physics_spec: dict,
        original_fortran: str = "",
        expected_behavior: str = ""
    ) -> GenerationResult:
        """
        Generate Python code with iterative repair.
        
        Args:
            design: Dict from ModernDesignAgent.
            physics_spec: Dict from PhysicsSpecAgent.
            original_fortran: Original Fortran code for reference.
            expected_behavior: Description of expected test results.
            
        Returns:
            GenerationResult with generated code and metrics.
        """
        user_prompt = self._build_initial_prompt(
            design, physics_spec, original_fortran
        )
        
        for iteration in range(1, self.max_iterations + 1):
            # Generate code
            response = self.llm_client.generate(
                system_prompt=self.system_prompt,
                user_prompt=user_prompt
            )
            code = self._extract_code(response)
            
            # Execute in sandbox
            exec_result = self._execute_code(code)
            
            if exec_result["success"]:
                # Check if verify() passed
                if exec_result.get("metrics", {}).get("passed", False):
                    return GenerationResult(
                        success=True,
                        code=code,
                        iteration=iteration,
                        metrics=exec_result.get("metrics", {})
                    )
                else:
                    # Code runs but verification failed
                    user_prompt = METRIC_FEEDBACK_TEMPLATE.format(
                        code=code,
                        metrics=json.dumps(exec_result.get("metrics", {}), indent=2),
                        expected=expected_behavior or "Stiffness matrix should be symmetric and positive semi-definite"
                    )
            else:
                # Code failed to execute
                user_prompt = ERROR_FEEDBACK_TEMPLATE.format(
                    code=code,
                    error=exec_result.get("error", "Unknown error")
                )
        
        # Max iterations reached
        return GenerationResult(
            success=False,
            code=code,
            iteration=self.max_iterations,
            error_message="Max iterations reached without passing verification"
        )
    
    def _build_initial_prompt(
        self, design: dict, physics_spec: dict, original_fortran: str
    ) -> str:
        """Build the initial generation prompt."""
        prompt = f"""Generate a complete Python implementation based on the following specification.

Design Specification:
{json.dumps(design, indent=2)}

Physics Specification:
{json.dumps(physics_spec, indent=2)}
"""
        if original_fortran:
            lines = original_fortran.split('\n')[:80]
            prompt += f"""
Original Fortran Code (for reference - replicate its behavior exactly):
```fortran
{chr(10).join(lines)}
```
"""
        prompt += """
Generate the complete, runnable Python code including the verify() function.
The verify() function must return a dict with at least {"passed": bool, "details": {...}}.
"""
        return prompt
    
    def _extract_code(self, response: str) -> str:
        """Extract Python code from LLM response."""
        if "```python" in response:
            code = response.split("```python")[1].split("```")[0].strip()
        elif "```" in response:
            code = response.split("```")[1].split("```")[0].strip()
        else:
            code = response.strip()
        return code
    
    def _execute_code(self, code: str) -> dict:
        """Execute code in sandbox and return results."""
        # Write to temp file
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.py', delete=False, encoding='utf-8'
        ) as f:
            # Wrap code to capture verify() output
            wrapped_code = code + """

# === Auto-appended verification runner ===
import json as _json
import sys as _sys

try:
    _result = verify()
    print("__RESULT_START__")
    print(_json.dumps(_result, default=str))
    print("__RESULT_END__")
except Exception as _e:
    print("__RESULT_START__")
    print(_json.dumps({"passed": False, "error": str(_e)}))
    print("__RESULT_END__")
"""
            f.write(wrapped_code)
            temp_path = f.name
        
        try:
            result = subprocess.run(
                [sys.executable, temp_path],
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=tempfile.gettempdir()
            )
            
            if result.returncode != 0:
                return {
                    "success": False,
                    "error": result.stderr[:2000]
                }
            
            # Extract metrics from output
            output = result.stdout
            if "__RESULT_START__" in output and "__RESULT_END__" in output:
                metrics_str = output.split("__RESULT_START__")[1].split("__RESULT_END__")[0].strip()
                try:
                    metrics = json.loads(metrics_str)
                except json.JSONDecodeError:
                    metrics = {"passed": False, "error": "Could not parse verify() output"}
            else:
                metrics = {"passed": False, "error": "No verification output found"}
            
            return {
                "success": True,
                "metrics": metrics,
                "stdout": output
            }
            
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "error": f"Execution timed out after {self.timeout} seconds"
            }
        except Exception as e:
            return {
                "success": False,
                "error": str(e)
            }
        finally:
            Path(temp_path).unlink(missing_ok=True)
