"""
Multi-Language Refactoring Experiments v2 — with Numerical Verification.

Key difference from v1: the agent loop now COMPILES and RUNS generated code,
compares output against analytical reference, and feeds back compile errors
or numerical deviations to the LLM for iterative repair.

Flow per run:
  1. Generate code (direct prompt or multi-agent)
  2. Compile (MSVC+Eigen for C++, cargo+nalgebra for Rust)
  3. If compile fails → feed error to LLM → repair → goto 2
  4. Run executable, capture JSON output
  5. Parse stiffness matrix, compare against analytical reference
  6. If numerical mismatch → feed deviation to LLM → repair → goto 2
  7. Max 3 repair iterations

Usage:
    python run_multilang_v2.py --lang cpp --element krod --model deepseek-v4-flash
    python run_multilang_v2.py --lang rust --element all --model deepseek-v4-flash
    python run_multilang_v2.py --lang all --element all --repeat 3
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).parent
sys.path.insert(0, str(PROJECT_ROOT))

from experiment_logger import ExperimentLogger, create_experiment_id
from llm_client import LLMClient
from agents.orchestrator import RefactoringOrchestrator

# ============================================================
# Tool paths
# ============================================================

MSVC_CL = r"C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC\14.36.32532\bin\Hostx64\x64\cl.exe"
MSVC_INCLUDE = r"C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC\14.36.32532\include"
MSVC_LIB = r"C:\Program Files\Microsoft Visual Studio\2022\Professional\VC\Tools\MSVC\14.36.32532\lib\x64"
EIGEN_INCLUDE = r"F:\eigen\eigen-3.4.0"

CARGO_HOME = r"F:\cargo"
RUSTUP_HOME = r"F:\rustup"
CARGO_BIN = os.path.join(CARGO_HOME, "bin")

# Fixed Rust build directory (reuses cached dependencies)
RUST_BUILD_DIR = PROJECT_ROOT / "rust_build"

# Windows SDK
def find_windows_sdk():
    sdk_base = r"C:\Program Files (x86)\Windows Kits\10"
    include_base = os.path.join(sdk_base, "Include")
    if not os.path.exists(include_base):
        return None, None
    versions = sorted(os.listdir(include_base), reverse=True)
    for v in versions:
        if os.path.exists(os.path.join(include_base, v, "ucrt")):
            return os.path.join(include_base, v), os.path.join(sdk_base, "Lib", v)
    return None, None

SDK_INCLUDE, SDK_LIB = find_windows_sdk()

# ============================================================
# Element configurations
# ============================================================

ELEMENTS = {
    "krod": {
        "name": "CROD",
        "source": "nastran/NASTRAN-95/mis/krod.f",
        "description": "Rod element (axial + torsion), 2 nodes, 6 DOF/node",
        "matrix_size": 12,
        "test_params": {
            "E": 200.0e9, "A": 0.01, "G": 76.923e9, "J": 5.0e-6,
            "node_a": [0.0, 0.0, 0.0], "node_b": [2.0, 0.0, 0.0],
        },
    },
    "kbar": {
        "name": "CBAR",
        "source": "nastran/NASTRAN-95/mis/kbar.f",
        "description": "Beam element (bending + torsion + axial), 2 nodes, 12x12 stiffness",
        "matrix_size": 12,
        "test_params": {
            "E": 200.0e9, "G": 76.923e9, "A": 0.01, "Iy": 8.333e-6, "Iz": 8.333e-6,
            "J": 1.667e-5, "L": 2.0,
            "node_a": [0.0, 0.0, 0.0], "node_b": [2.0, 0.0, 0.0],
        },
    },
    "ktrmem": {
        "name": "CTRIA3",
        "source": "nastran/NASTRAN-95/mis/ktrmem.f",
        "description": "Triangular membrane (CST), 3 nodes, 2 DOF/node, 6x6 stiffness",
        "matrix_size": 6,
        "test_params": {
            "E": 2.1e11, "nu": 0.3, "t": 0.01,
            "nodes": [[0.0, 0.0], [2.0, 0.0], [1.0, 1.5]],
        },
    },
    "ktriqd": {
        "name": "CQUAD4",
        "source": ["nastran/NASTRAN-95/mis/ktriqd.f", "nastran/NASTRAN-95/mis/ktrmem.f"],
        "description": "Quad membrane (2-triangle decomposition), 4 nodes, 8x8 stiffness",
        "matrix_size": 8,
        "test_params": {
            "E": 200.0e9, "nu": 0.3, "t": 0.01,
            "nodes": [[0.0, 0.0], [2.0, 0.0], [2.0, 1.5], [0.0, 1.5]],
        },
    },
}

# ============================================================
# Analytical Reference Computation
# ============================================================

def compute_reference(element_key):
    """Compute the analytical reference stiffness matrix."""
    p = ELEMENTS[element_key]["test_params"]

    if element_key == "krod":
        na, nb = np.array(p["node_a"]), np.array(p["node_b"])
        L = np.linalg.norm(nb - na)
        n = (nb - na) / L
        nnT = np.outer(n, n)
        ea_l = p["E"] * p["A"] / L
        gj_l = p["G"] * p["J"] / L
        block = np.zeros((6, 6))
        block[0:3, 0:3] = ea_l * nnT
        block[3:6, 3:6] = gj_l * nnT
        K = np.zeros((12, 12))
        K[0:6, 0:6] = block
        K[0:6, 6:12] = -block
        K[6:12, 0:6] = -block
        K[6:12, 6:12] = block
        return K

    elif element_key == "kbar":
        E, G, A = p["E"], p["G"], p["A"]
        Iy, Iz, J, L = p["Iy"], p["Iz"], p["J"], p["L"]
        K = np.zeros((12, 12))
        ea_l = E * A / L
        K[0,0] = K[6,6] = ea_l; K[0,6] = K[6,0] = -ea_l
        gj_l = G * J / L
        K[3,3] = K[9,9] = gj_l; K[3,9] = K[9,3] = -gj_l
        # Bending xz (Iz): DOF 1,5,7,11
        c = E * Iz / L**3
        K[1,1]=12*c; K[1,5]=6*c*L; K[1,7]=-12*c; K[1,11]=6*c*L
        K[5,1]=6*c*L; K[5,5]=4*c*L**2; K[5,7]=-6*c*L; K[5,11]=2*c*L**2
        K[7,1]=-12*c; K[7,5]=-6*c*L; K[7,7]=12*c; K[7,11]=-6*c*L
        K[11,1]=6*c*L; K[11,5]=2*c*L**2; K[11,7]=-6*c*L; K[11,11]=4*c*L**2
        # Bending xy (Iy): DOF 2,4,8,10
        c2 = E * Iy / L**3
        K[2,2]=12*c2; K[2,4]=-6*c2*L; K[2,8]=-12*c2; K[2,10]=-6*c2*L
        K[4,2]=-6*c2*L; K[4,4]=4*c2*L**2; K[4,8]=6*c2*L; K[4,10]=2*c2*L**2
        K[8,2]=-12*c2; K[8,4]=6*c2*L; K[8,8]=12*c2; K[8,10]=6*c2*L
        K[10,2]=-6*c2*L; K[10,4]=2*c2*L**2; K[10,8]=6*c2*L; K[10,10]=4*c2*L**2
        return K

    elif element_key == "ktrmem":
        E, nu, t = p["E"], p["nu"], p["t"]
        nodes = p["nodes"]
        x1,y1 = nodes[0]; x2,y2 = nodes[1]; x3,y3 = nodes[2]
        A = 0.5 * abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1))
        b1,b2,b3 = y2-y3, y3-y1, y1-y2
        c1,c2,c3 = x3-x2, x1-x3, x2-x1
        B = (1/(2*A)) * np.array([[b1,0,b2,0,b3,0],[0,c1,0,c2,0,c3],[c1,b1,c2,b2,c3,b3]])
        D = (E/(1-nu**2)) * np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
        return t * A * (B.T @ D @ B)

    elif element_key == "ktriqd":
        E, nu, t = p["E"], p["nu"], p["t"]
        nodes = p["nodes"]
        def cst_k(n0, n1, n2):
            x1,y1=n0; x2,y2=n1; x3,y3=n2
            A=0.5*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
            b1,b2,b3=y2-y3,y3-y1,y1-y2
            c1,c2,c3=x3-x2,x1-x3,x2-x1
            B=(1/(2*A))*np.array([[b1,0,b2,0,b3,0],[0,c1,0,c2,0,c3],[c1,b1,c2,b2,c3,b3]])
            D=(E/(1-nu**2))*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
            return t*A*(B.T@D@B)
        K1 = cst_k(nodes[0], nodes[1], nodes[2])
        K2 = cst_k(nodes[0], nodes[2], nodes[3])
        K = np.zeros((8, 8))
        idx1 = [0,1,2,3,4,5]
        for i,ii in enumerate(idx1):
            for j,jj in enumerate(idx1):
                K[ii,jj] += K1[i,j]
        idx2 = [0,1,4,5,6,7]
        for i,ii in enumerate(idx2):
            for j,jj in enumerate(idx2):
                K[ii,jj] += K2[i,j]
        return K

    return None

# ============================================================
# Compilation & Execution
# ============================================================

# Windows fatal exit codes that abort a process without writing to stderr.
# Without this table a crash was reported as the bare string "runtime: ", the
# repair prompt therefore said nothing actionable, and all five repair
# iterations produced byte-identical code. Observed on a B4 artifact whose
# exit code was 0xC0000409.
_WIN_FATAL = {
    0xC0000005: ("ACCESS_VIOLATION",
                 "dereferenced or wrote through an invalid pointer, or indexed "
                 "past the end of an array"),
    0xC0000094: ("INTEGER_DIVIDE_BY_ZERO", "integer division by zero"),
    0xC0000095: ("INTEGER_OVERFLOW", "integer overflow"),
    0xC00000FD: ("STACK_OVERFLOW",
                 "stack exhausted, usually unbounded recursion or a very large "
                 "stack-allocated array"),
    0xC0000409: ("STACK_BUFFER_OVERRUN",
                 "MSVC fail-fast: a stack buffer was written out of bounds, or "
                 "a security check tripped. With Eigen this is typically a "
                 "write past the end of a fixed-size Matrix or a raw array"),
    0xC0000374: ("HEAP_CORRUPTION", "the heap was corrupted by an out-of-bounds "
                 "write"),
    0xC000041D: ("FATAL_USER_CALLBACK_EXCEPTION",
                 "an exception escaped a callback"),
    0xC0000417: ("INVALID_CRT_PARAMETER",
                 "a CRT function was called with an invalid parameter, e.g. an "
                 "out-of-range index passed to a checked accessor"),
    3: ("ABORT", "abort() was called, e.g. an uncaught exception or a failed "
        "assert"),
}


def _describe_crash(returncode: int, stderr: str, stdout: str) -> str:
    """
    Turn a nonzero exit into something a repair prompt can act on.

    A process killed by a Windows fatal exception writes nothing to stderr, so
    passing stderr through verbatim yielded an empty diagnostic.
    """
    parts = []
    unsigned = returncode & 0xFFFFFFFF
    name = _WIN_FATAL.get(unsigned) or _WIN_FATAL.get(returncode)
    if name:
        parts.append(f"process terminated by {name[0]} "
                     f"(exit 0x{unsigned:08X}): {name[1]}")
    else:
        parts.append(f"process exited with code {returncode} "
                     f"(0x{unsigned:08X})")
    if stderr and stderr.strip():
        parts.append("stderr:\n" + stderr[:800])
    else:
        parts.append("stderr was empty, so the failure is a hard crash rather "
                     "than a reported error")
    if stdout and stdout.strip():
        tail = stdout.strip().splitlines()[-12:]
        parts.append("last output before the crash:\n" + "\n".join(tail)[:600])
    else:
        parts.append("nothing was written to stdout, so the crash happened "
                     "before the first print")
    return "\n".join(parts)[:1600]


def compile_and_run_cpp(source_code: str, work_dir: Path) -> dict:
    """Compile and run C++ code. Returns {success, stdout, error, stage}."""
    src_file = work_dir / "main.cpp"
    exe_file = work_dir / "main.exe"
    src_file.write_text(source_code, encoding='utf-8')

    include_dirs = [EIGEN_INCLUDE, MSVC_INCLUDE]
    lib_dirs = [MSVC_LIB]
    if SDK_INCLUDE:
        include_dirs += [os.path.join(SDK_INCLUDE,"ucrt"), os.path.join(SDK_INCLUDE,"um"), os.path.join(SDK_INCLUDE,"shared")]
    if SDK_LIB:
        lib_dirs += [os.path.join(SDK_LIB,"ucrt","x64"), os.path.join(SDK_LIB,"um","x64")]

    cmd = [MSVC_CL, "/EHsc", "/std:c++17", "/O2", "/nologo", f"/Fe:{exe_file}", str(src_file)]
    for inc in include_dirs:
        cmd.append(f"/I{inc}")
    cmd.append("/link")
    for lib in lib_dirs:
        cmd.append(f"/LIBPATH:{lib}")

    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=60, cwd=str(work_dir))
        if r.returncode != 0:
            return {"success": False, "stage": "compile", "error": (r.stderr + r.stdout)[:2000]}
    except subprocess.TimeoutExpired:
        return {"success": False, "stage": "compile", "error": "Compilation timed out (60s)"}
    except Exception as e:
        return {"success": False, "stage": "compile", "error": str(e)}

    try:
        r = subprocess.run([str(exe_file)], capture_output=True, text=True, timeout=30, cwd=str(work_dir))
        if r.returncode != 0:
            return {"success": False, "stage": "runtime",
                    "error": _describe_crash(r.returncode, r.stderr, r.stdout),
                    "stdout": r.stdout}
        return {"success": True, "stage": "run", "stdout": r.stdout, "error": ""}
    except subprocess.TimeoutExpired:
        return {"success": False, "stage": "runtime", "error": "Execution timed out"}
    except Exception as e:
        return {"success": False, "stage": "runtime", "error": str(e)}


def compile_and_run_rust(source_code: str) -> dict:
    """Compile and run Rust code using fixed build dir. Returns {success, stdout, error, stage}."""
    cargo_dir = RUST_BUILD_DIR
    src_dir = cargo_dir / "src"
    src_dir.mkdir(parents=True, exist_ok=True)

    cargo_toml = '[package]\nname = "nastran_verify"\nversion = "0.1.0"\nedition = "2021"\n\n[dependencies]\nnalgebra = "0.33"\nserde_json = "1.0"\n'
    (cargo_dir / "Cargo.toml").write_text(cargo_toml, encoding='utf-8')
    (src_dir / "main.rs").write_text(source_code, encoding='utf-8')

    env = os.environ.copy()
    env["RUSTUP_HOME"] = RUSTUP_HOME
    env["CARGO_HOME"] = CARGO_HOME
    env["PATH"] = CARGO_BIN + ";" + env.get("PATH", "")

    try:
        r = subprocess.run([os.path.join(CARGO_BIN,"cargo.exe"), "build", "--release"],
                          capture_output=True, text=True, timeout=300, cwd=str(cargo_dir), env=env)
        if r.returncode != 0:
            return {"success": False, "stage": "compile", "error": r.stderr[:2000]}
    except subprocess.TimeoutExpired:
        return {"success": False, "stage": "compile", "error": "Compilation timed out (300s)"}
    except Exception as e:
        return {"success": False, "stage": "compile", "error": str(e)}

    exe_path = cargo_dir / "target" / "release" / "nastran_verify.exe"
    if not exe_path.exists():
        return {"success": False, "stage": "compile", "error": "Executable not found after build"}

    try:
        r = subprocess.run([str(exe_path)], capture_output=True, text=True, timeout=30, cwd=str(cargo_dir), env=env)
        if r.returncode != 0:
            return {"success": False, "stage": "runtime", "error": r.stderr[:1000], "stdout": r.stdout}
        return {"success": True, "stage": "run", "stdout": r.stdout, "error": ""}
    except subprocess.TimeoutExpired:
        return {"success": False, "stage": "runtime", "error": "Execution timed out"}
    except Exception as e:
        return {"success": False, "stage": "runtime", "error": str(e)}

# ============================================================
# Output Parsing & Numerical Comparison
# ============================================================

def parse_matrix_from_output(stdout: str, size: int) -> np.ndarray:
    """Parse stiffness matrix from program JSON output."""
    # Strategy 1: find JSON with stiffness_matrix key
    try:
        brace_start = stdout.find('{')
        if brace_start >= 0:
            depth = 0
            for i in range(brace_start, len(stdout)):
                if stdout[i] == '{': depth += 1
                elif stdout[i] == '}':
                    depth -= 1
                    if depth == 0:
                        data = json.loads(stdout[brace_start:i+1])
                        for key in ['stiffness_matrix','matrix','K','stiffness_matrix_12x12','stiffness_matrix_6x6']:
                            if key in data:
                                return np.array(data[key], dtype=float)
                        break
    except (json.JSONDecodeError, ValueError):
        pass

    # Strategy 2: find [[...],[...]] pattern
    try:
        arr_match = re.search(r'\[\s*\[[\d\s.,eE+\-]+\](?:\s*,?\s*\[[\d\s.,eE+\-]+\])*\s*\]', stdout)
        if arr_match:
            return np.array(json.loads(arr_match.group()), dtype=float)
    except:
        pass

    # Strategy 3: parse rows of numbers
    try:
        lines = stdout.strip().split('\n')
        rows = []
        for line in lines:
            cleaned = line.strip().strip('[]').replace(',', ' ')
            nums = re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?', cleaned)
            if len(nums) >= size:
                rows.append([float(x) for x in nums[:size]])
        if len(rows) >= size:
            return np.array(rows[:size], dtype=float)
    except:
        pass

    return None


def numerical_check(computed: np.ndarray, reference: np.ndarray, tol=1e-4) -> dict:
    """Compare computed vs reference matrix. Returns {passed, max_rel_diff, details}."""
    if computed is None:
        return {"passed": False, "error": "Could not parse matrix from output"}

    # Handle size mismatch (12x12 output vs 6x6 reference or vice versa)
    if computed.shape != reference.shape:
        if computed.shape == (12, 12) and reference.shape == (6, 6):
            idx = [0,1,2,6,7,8]
            computed = computed[np.ix_(idx, idx)]
        elif computed.shape[0] != reference.shape[0]:
            return {"passed": False, "error": f"Shape mismatch: got {computed.shape}, expected {reference.shape}"}

    max_ref = max(np.max(np.abs(reference)), 1e-15)
    abs_diff = np.max(np.abs(computed - reference))
    rel_diff = abs_diff / max_ref

    passed = rel_diff < tol or abs_diff < 1e-10
    return {
        "passed": passed,
        "max_abs_diff": float(abs_diff),
        "max_rel_diff": float(rel_diff),
        "is_symmetric": bool(np.allclose(computed, computed.T, atol=1e-10)),
    }

# ============================================================
# Prompt Templates (with mandatory test case)
# ============================================================

def get_test_code_requirement(lang: str, element_key: str) -> str:
    """Generate the mandatory test case code requirement for the prompt."""
    p = ELEMENTS[element_key]["test_params"]
    size = ELEMENTS[element_key]["matrix_size"]

    if element_key == "krod":
        if lang == "cpp":
            return f"""
Your main() MUST create a rod element with these EXACT parameters and print the result as JSON:
- node_a = ({p['node_a'][0]}, {p['node_a'][1]}, {p['node_a'][2]})
- node_b = ({p['node_b'][0]}, {p['node_b'][1]}, {p['node_b'][2]})
- E = {p['E']}, A = {p['A']}, G = {p['G']}, J = {p['J']}
- Expected matrix size: {size}x{size}

Output format (MUST be valid JSON to stdout):
{{"test": "CROD", "passed": true/false, "stiffness_matrix": [[row0], [row1], ...]}}
"""
        else:  # rust
            return f"""
Your main() MUST create a rod element with these EXACT parameters and print the result as JSON:
- node_a = [{p['node_a'][0]}, {p['node_a'][1]}, {p['node_a'][2]}]
- node_b = [{p['node_b'][0]}, {p['node_b'][1]}, {p['node_b'][2]}]
- E = {p['E']}, A = {p['A']}, G = {p['G']}, J = {p['J']}
- Expected matrix size: {size}x{size}

Output format (MUST be valid JSON to stdout):
{{"test": "CROD", "passed": true, "stiffness_matrix": [[row0], [row1], ...]}}
"""

    elif element_key == "kbar":
        if lang == "cpp":
            return f"""
Your main() MUST create a beam element along x-axis with these EXACT parameters:
- node_a = (0,0,0), node_b = ({p['L']},0,0), L = {p['L']}
- E = {p['E']}, G = {p['G']}, A = {p['A']}
- Iy = {p['Iy']}, Iz = {p['Iz']}, J = {p['J']}
- Expected matrix size: 12x12

Output format (MUST be valid JSON to stdout):
{{"test": "CBAR", "stiffness_matrix": [[row0], ..., [row11]]}}
"""
        else:
            return f"""
Your main() MUST create a beam element along x-axis with these EXACT parameters:
- node_a = [0,0,0], node_b = [{p['L']},0,0], L = {p['L']}
- E = {p['E']}, G = {p['G']}, A = {p['A']}
- Iy = {p['Iy']}, Iz = {p['Iz']}, J = {p['J']}
- Expected matrix size: 12x12

Output format (MUST be valid JSON to stdout):
{{"test": "CBAR", "stiffness_matrix": [[row0], ..., [row11]]}}
"""

    elif element_key == "ktrmem":
        nodes = p["nodes"]
        if lang == "cpp":
            return f"""
Your main() MUST create a CST triangle element with these EXACT parameters:
- node1 = ({nodes[0][0]}, {nodes[0][1]}), node2 = ({nodes[1][0]}, {nodes[1][1]}), node3 = ({nodes[2][0]}, {nodes[2][1]})
- E = {p['E']}, nu = {p['nu']}, thickness = {p['t']}
- Expected matrix size: 6x6

Output format (MUST be valid JSON to stdout):
{{"test": "CTRIA3", "stiffness_matrix": [[row0], ..., [row5]]}}
"""
        else:
            return f"""
Your main() MUST create a CST triangle element with these EXACT parameters:
- node1 = [{nodes[0][0]}, {nodes[0][1]}], node2 = [{nodes[1][0]}, {nodes[1][1]}], node3 = [{nodes[2][0]}, {nodes[2][1]}]
- E = {p['E']}, nu = {p['nu']}, thickness = {p['t']}
- Expected matrix size: 6x6

Output format (MUST be valid JSON to stdout):
{{"test": "CTRIA3", "stiffness_matrix": [[row0], ..., [row5]]}}
"""

    elif element_key == "ktriqd":
        nodes = p["nodes"]
        if lang == "cpp":
            return f"""
Your main() MUST create a quad membrane element (decomposed into 2 CST triangles) with:
- node1=({nodes[0][0]},{nodes[0][1]}), node2=({nodes[1][0]},{nodes[1][1]}), node3=({nodes[2][0]},{nodes[2][1]}), node4=({nodes[3][0]},{nodes[3][1]})
- E = {p['E']}, nu = {p['nu']}, thickness = {p['t']}
- Decompose into triangles (1-2-3) and (1-3-4), assemble into 8x8

Output format (MUST be valid JSON to stdout):
{{"test": "CQUAD4", "stiffness_matrix": [[row0], ..., [row7]]}}
"""
        else:
            return f"""
Your main() MUST create a quad membrane element (decomposed into 2 CST triangles) with:
- node1=[{nodes[0][0]},{nodes[0][1]}], node2=[{nodes[1][0]},{nodes[1][1]}], node3=[{nodes[2][0]},{nodes[2][1]}], node4=[{nodes[3][0]},{nodes[3][1]}]
- E = {p['E']}, nu = {p['nu']}, thickness = {p['t']}
- Decompose into triangles (1-2-3) and (1-3-4), assemble into 8x8

Output format (MUST be valid JSON to stdout):
{{"test": "CQUAD4", "stiffness_matrix": [[row0], ..., [row7]]}}
"""
    return ""

def build_generation_prompt(lang: str, element_key: str, fortran_code: str) -> str:
    """Build the initial code generation prompt."""
    elem = ELEMENTS[element_key]
    test_req = get_test_code_requirement(lang, element_key)

    if lang == "cpp":
        return f"""Translate the following NASTRAN-95 Fortran 77 subroutine to modern C++17.

Element: {elem['description']}

Requirements:
1. Use Eigen library for matrix operations (#include <Eigen/Dense>).
2. Create a class with a method that returns the stiffness matrix as Eigen::MatrixXd.
3. Include a main() function that runs the test case below and prints JSON output.
4. The code must compile with MSVC (cl.exe /std:c++17 /EHsc).
5. Use double for all floating point.
6. Include ALL necessary headers (#include <iostream>, <iomanip>, <Eigen/Dense>).
7. Do NOT use any external libraries other than Eigen and the C++ standard library.

{test_req}

Fortran source:
```fortran
{fortran_code}
```

Generate COMPLETE C++ code (single file). The JSON output must contain the full stiffness matrix as a 2D array.
"""
    else:  # rust
        return f"""Translate the following NASTRAN-95 Fortran 77 subroutine to modern Rust.

Element: {elem['description']}

Requirements:
1. Use nalgebra crate for matrix operations (nalgebra = "0.33" in Cargo.toml).
2. Create a struct with a method that returns the stiffness matrix as DMatrix<f64>.
3. Include a main() function that runs the test case below and prints JSON output.
4. The code must compile with: cargo build --release
5. Use f64 for all floating point.
6. You may use serde_json = "1.0" for JSON output, or manual println!.
7. Do NOT use any crates other than nalgebra and serde_json.

{test_req}

Fortran source:
```fortran
{fortran_code}
```

Generate COMPLETE Rust code (single file, main.rs). The JSON output must contain the full stiffness matrix as a 2D array.
"""

def build_repair_prompt(lang: str, element_key: str, code: str, error_info: dict) -> str:
    """Build a repair prompt based on the specific failure."""
    stage = error_info.get("stage", "unknown")
    error = error_info.get("error", "")
    test_req = get_test_code_requirement(lang, element_key)

    if stage == "compile":
        return f"""The following {lang.upper()} code failed to compile. Fix ALL compilation errors.

Compiler error:
```
{error[:1500]}
```

Current code:
```{lang if lang != 'cpp' else 'cpp'}
{code}
```

Requirements:
- Fix the compilation errors
- Keep the same algorithm and test case
- {'Use Eigen (#include <Eigen/Dense>) for matrices' if lang == 'cpp' else 'Use nalgebra (DMatrix<f64>) for matrices'}
{test_req}

Generate the COMPLETE fixed code (entire file).
"""

    elif stage == "runtime":
        return f"""The following {lang.upper()} code compiled but crashed at runtime.

Runtime error:
```
{error[:1000]}
```

Current code:
```{lang if lang != 'cpp' else 'cpp'}
{code}
```

Fix the runtime error and generate the COMPLETE fixed code.
{test_req}
"""

    elif stage == "numerical":
        num_info = error_info.get("numerical", {})
        ref_matrix = error_info.get("reference_preview", "")
        return f"""The following {lang.upper()} code compiles and runs, but produces INCORRECT numerical results.

Numerical deviation:
- Max relative error: {num_info.get('max_rel_diff', 'N/A')}
- Max absolute error: {num_info.get('max_abs_diff', 'N/A')}
- Expected matrix (first 3 rows): {ref_matrix}

The stiffness matrix computation is wrong. Please check:
- Direction cosines / coordinate transformations
- Matrix assembly (signs, indices)
- Physical formulas (EA/L for axial, GJ/L for torsion, 12EI/L^3 for bending, etc.)

Current code:
```{lang if lang != 'cpp' else 'cpp'}
{code}
```

{test_req}

Generate the COMPLETE fixed code with correct physics.
"""

    else:
        return f"""The following {lang.upper()} code has issues. Please fix it.

Error: {error[:500]}

Current code:
```{lang if lang != 'cpp' else 'cpp'}
{code}
```

{test_req}

Generate the COMPLETE fixed code.
"""

# ============================================================
# Code extraction helper
# ============================================================

def extract_code(response: str, lang: str) -> str:
    """Extract code block from LLM response."""
    markers = {
        "rust": ["```rust", "```rs"],
        "cpp": ["```cpp", "```c++", "```c"],
    }
    for marker in markers.get(lang, ["```"]):
        if marker in response:
            parts = response.split(marker, 1)
            if len(parts) > 1:
                return parts[1].split("```")[0].strip()
    # Fallback: try generic code block
    if "```" in response:
        parts = response.split("```", 2)
        if len(parts) >= 3:
            # Skip the language tag on first line
            code = parts[1]
            first_newline = code.find('\n')
            if first_newline > 0 and first_newline < 20:
                code = code[first_newline+1:]
            return code.strip()
    return response.strip()

# ============================================================
# Main experiment runner (with compile-run-verify loop)
# ============================================================

MAX_REPAIR_ITERATIONS = 3


def run_single_experiment(model_name: str, lang: str, element_key: str,
                          fortran_code: str, method: str, logger, run_idx: int) -> dict:
    """
    Run a single experiment with integrated numerical verification.
    
    Returns dict with: passed, compile_ok, run_ok, numerical_ok, iterations, etc.
    """
    logger.log_info(f"--- {lang.upper()} | {ELEMENTS[element_key]['name']} | {method} | Run {run_idx} ---")
    logger.log_stage_start(f"{method}_{lang}_{element_key}_run{run_idx}")

    llm = LLMClient(model=model_name, logger=logger)
    reference = compute_reference(element_key)
    mat_size = ELEMENTS[element_key]["matrix_size"]

    # Step 1: Generate initial code
    if method == "direct":
        prompt = build_generation_prompt(lang, element_key, fortran_code)
        response = llm.generate(
            system_prompt=f"You are an expert in Fortran, {'C++17 with Eigen' if lang == 'cpp' else 'Rust with nalgebra'}, and finite element analysis. Generate complete, compilable code.",
            user_prompt=prompt,
            agent_name=f"direct_{lang}_{element_key}"
        )
        code = extract_code(response, lang)
    else:
        # Multi-agent: use orchestrator for comprehension/design, then generate
        orchestrator = RefactoringOrchestrator(
            llm_client=llm, logger=logger,
            output_dir=str(PROJECT_ROOT / "results" / logger.experiment_id / f"run_{run_idx}")
        )
        state = orchestrator.run(
            source_files={f"{element_key}.f": fortran_code},
            context=ELEMENTS[element_key]["description"] + f" Target: {lang.upper()}.",
            reference_results=[],
        )
        if not state.modern_design:
            logger.log_info("  No design produced, falling back to direct prompt")
            prompt = build_generation_prompt(lang, element_key, fortran_code)
        else:
            prompt = build_generation_prompt(lang, element_key, fortran_code)
            prompt += f"\n\nAdditional design context from analysis:\n{json.dumps(state.modern_design, indent=2, default=str)[:3000]}"

        response = llm.generate(
            system_prompt=f"You are an expert in {'C++17 with Eigen' if lang == 'cpp' else 'Rust with nalgebra'} and FEA.",
            user_prompt=prompt,
            agent_name=f"codegen_{lang}_{element_key}"
        )
        code = extract_code(response, lang)

    if len(code) < 50:
        logger.log_info("  ERROR: Generated code too short")
        logger.log_stage_end(f"{method}_{lang}_{element_key}_run{run_idx}", False)
        return {"passed": False, "compile_ok": False, "run_ok": False, "numerical_ok": False,
                "iterations": 0, "error": "Code too short"}

    # Step 2: Compile-Run-Verify loop with repair
    for iteration in range(MAX_REPAIR_ITERATIONS + 1):
        logger.log_info(f"  Iteration {iteration}: compiling and running...")
        logger.log_artifact(f"code_{lang}_{element_key}_run{run_idx}_iter{iteration}{'.' + ('rs' if lang == 'rust' else 'cpp')}", code, is_json=False)

        # Compile and run
        if lang == "cpp":
            with tempfile.TemporaryDirectory(prefix="nastran_cpp_") as tmpdir:
                result = compile_and_run_cpp(code, Path(tmpdir))
        else:
            result = compile_and_run_rust(code)

        if not result["success"]:
            stage = result["stage"]
            error = result["error"]
            logger.log_info(f"  FAIL at {stage}: {error[:200]}")

            if iteration >= MAX_REPAIR_ITERATIONS:
                logger.log_info(f"  Max repair iterations reached. Final failure at {stage}.")
                logger.log_stage_end(f"{method}_{lang}_{element_key}_run{run_idx}", False)
                return {"passed": False, "compile_ok": stage != "compile",
                        "run_ok": False, "numerical_ok": False,
                        "iterations": iteration, "error": f"{stage}: {error[:200]}"}

            # Repair
            logger.log_info(f"  Requesting repair (iteration {iteration+1})...")
            repair_prompt = build_repair_prompt(lang, element_key, code, result)
            response = llm.generate(
                system_prompt=f"Fix the {lang.upper()} code. Return the COMPLETE fixed file.",
                user_prompt=repair_prompt,
                agent_name=f"repair_{lang}_{element_key}_iter{iteration+1}"
            )
            new_code = extract_code(response, lang)
            if len(new_code) > 50:
                code = new_code
            continue

        # Compiled and ran successfully — now check numerical output
        stdout = result.get("stdout", "")
        logger.log_info(f"  Compiled & ran OK. Output length: {len(stdout)}")

        computed = parse_matrix_from_output(stdout, mat_size)
        num_result = numerical_check(computed, reference)

        if num_result["passed"]:
            logger.log_info(f"  ✓ NUMERICAL PASS (rel_diff={num_result.get('max_rel_diff', 0):.2e})")
            logger.log_artifact(f"final_{lang}_{element_key}_run{run_idx}{'.' + ('rs' if lang == 'rust' else 'cpp')}", code, is_json=False)
            logger.log_artifact(f"verify_{lang}_{element_key}_run{run_idx}", num_result)
            logger.log_stage_end(f"{method}_{lang}_{element_key}_run{run_idx}", True)
            return {"passed": True, "compile_ok": True, "run_ok": True, "numerical_ok": True,
                    "iterations": iteration, "max_rel_diff": num_result.get("max_rel_diff", 0)}

        # Numerical mismatch
        logger.log_info(f"  NUMERICAL FAIL: rel_diff={num_result.get('max_rel_diff', 'N/A')}")

        if iteration >= MAX_REPAIR_ITERATIONS:
            logger.log_info(f"  Max repair iterations reached. Numerical verification failed.")
            logger.log_artifact(f"final_{lang}_{element_key}_run{run_idx}{'.' + ('rs' if lang == 'rust' else 'cpp')}", code, is_json=False)
            logger.log_artifact(f"verify_{lang}_{element_key}_run{run_idx}", num_result)
            logger.log_stage_end(f"{method}_{lang}_{element_key}_run{run_idx}", False)
            return {"passed": False, "compile_ok": True, "run_ok": True, "numerical_ok": False,
                    "iterations": iteration, "max_rel_diff": num_result.get("max_rel_diff", 0)}

        # Repair with numerical feedback
        ref_preview = str(reference[:3, :3].tolist()) if reference is not None else "N/A"
        error_info = {
            "stage": "numerical",
            "error": f"Relative error {num_result.get('max_rel_diff', 'N/A')}",
            "numerical": num_result,
            "reference_preview": ref_preview,
        }
        logger.log_info(f"  Requesting numerical repair (iteration {iteration+1})...")
        repair_prompt = build_repair_prompt(lang, element_key, code, error_info)
        response = llm.generate(
            system_prompt=f"Fix the numerical error in the {lang.upper()} FEA code. The stiffness matrix values are wrong.",
            user_prompt=repair_prompt,
            agent_name=f"repair_num_{lang}_{element_key}_iter{iteration+1}"
        )
        new_code = extract_code(response, lang)
        if len(new_code) > 50:
            code = new_code

    # Should not reach here
    logger.log_stage_end(f"{method}_{lang}_{element_key}_run{run_idx}", False)
    return {"passed": False, "compile_ok": False, "run_ok": False, "numerical_ok": False,
            "iterations": MAX_REPAIR_ITERATIONS, "error": "Unexpected exit"}

# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Multi-Language Refactoring v2 (with numerical verification)")
    parser.add_argument("--lang", choices=["rust", "cpp", "all"], default="all")
    parser.add_argument("--element", choices=list(ELEMENTS.keys()) + ["all"], default="all")
    parser.add_argument("--model", default="deepseek-v4-flash")
    parser.add_argument("--method", choices=["direct", "multiagent", "both"], default="both")
    parser.add_argument("--repeat", type=int, default=3)
    args = parser.parse_args()

    languages = ["rust", "cpp"] if args.lang == "all" else [args.lang]
    elements = list(ELEMENTS.keys()) if args.element == "all" else [args.element]
    methods = ["direct", "multiagent"] if args.method == "both" else [args.method]

    print(f"{'='*70}")
    print(f"Multi-Language NASTRAN Refactoring v2 (Numerical Verification)")
    print(f"{'='*70}")
    print(f"Languages: {languages}")
    print(f"Elements: {elements}")
    print(f"Methods: {methods}")
    print(f"Model: {args.model}")
    print(f"Repeats: {args.repeat}")
    print(f"Max repair iterations: {MAX_REPAIR_ITERATIONS}")
    print(f"Time: {datetime.now().isoformat()}")
    print(f"{'='*70}\n")

    all_results = []

    for lang in languages:
        for element_key in elements:
            # Load Fortran source
            elem_config = ELEMENTS[element_key]
            source_paths = elem_config.get("source", "")
            if isinstance(source_paths, str):
                source_paths = [source_paths]

            fortran_code = ""
            for sp in source_paths:
                fpath = PROJECT_ROOT / sp
                if fpath.exists():
                    fortran_code += fpath.read_text(encoding='utf-8', errors='replace') + "\n\n"

            if not fortran_code:
                print(f"  SKIP: {element_key} - source not found")
                continue

            for method in methods:
                exp_id = create_experiment_id(f"{lang}_{method}_{args.model}", f"multilang_v2_{element_key}")
                logger = ExperimentLogger(exp_id, base_dir=str(PROJECT_ROOT / "results"))
                logger.log_info(f"Language: {lang}, Element: {element_key}, Method: {method}, Model: {args.model}")
                logger.log_artifact("source_fortran.f", fortran_code, is_json=False)

                results = []
                for run_idx in range(1, args.repeat + 1):
                    r = run_single_experiment(args.model, lang, element_key, fortran_code, method, logger, run_idx)
                    r["lang"] = lang
                    r["element"] = element_key
                    r["method"] = method
                    r["run"] = run_idx
                    results.append(r)

                    status = "✓" if r["passed"] else "✗"
                    detail = f"iter={r.get('iterations',0)}"
                    if r.get("max_rel_diff") is not None:
                        detail += f" rel_diff={r['max_rel_diff']:.2e}"
                    print(f"  {status} {lang:4s} | {elem_config['name']:6s} | {method:10s} | run {run_idx}: {detail}")

                successes = sum(1 for r in results if r["passed"])
                summary = f"{lang} {element_key} {method}: {successes}/{args.repeat}"
                logger.log_info(f"Summary: {summary}")
                logger.log_artifact("results_summary", {"results": results, "success_rate": successes/args.repeat})
                logger.finalize(successes > 0, summary)
                all_results.extend(results)

    # Final summary
    print(f"\n{'='*70}")
    print("FINAL SUMMARY (Numerical Verification)")
    print(f"{'='*70}")
    print(f"{'Lang':<6} {'Element':<8} {'Method':<12} {'Compile':<9} {'Run':<5} {'Numerical':<10} {'Pass':<6}")
    print("-" * 60)

    from itertools import groupby
    sorted_results = sorted(all_results, key=lambda x: (x["lang"], x["element"], x["method"]))
    for (lang, elem, method), group in groupby(sorted_results, key=lambda x: (x["lang"], x["element"], x["method"])):
        items = list(group)
        n = len(items)
        c = sum(1 for r in items if r.get("compile_ok"))
        ru = sum(1 for r in items if r.get("run_ok"))
        nu = sum(1 for r in items if r.get("numerical_ok"))
        p = sum(1 for r in items if r.get("passed"))
        print(f"{lang:<6} {elem:<8} {method:<12} {c}/{n:<5} {ru}/{n:<3} {nu}/{n:<7} {p}/{n}")

    total = len(all_results)
    total_pass = sum(1 for r in all_results if r.get("passed"))
    total_compile = sum(1 for r in all_results if r.get("compile_ok"))
    total_run = sum(1 for r in all_results if r.get("run_ok"))
    total_num = sum(1 for r in all_results if r.get("numerical_ok"))

    print(f"\n{'─'*60}")
    print(f"TOTAL: {total} experiments")
    print(f"  Compile: {total_compile}/{total} ({total_compile/total*100:.1f}%)" if total else "")
    print(f"  Run:     {total_run}/{total} ({total_run/total*100:.1f}%)" if total else "")
    print(f"  Numerical: {total_num}/{total} ({total_num/total*100:.1f}%)" if total else "")
    print(f"  PASS:    {total_pass}/{total} ({total_pass/total*100:.1f}%)" if total else "")

    # Save
    summary_file = PROJECT_ROOT / "results" / "multilang_v2_summary.json"
    summary_file.write_text(json.dumps({
        "timestamp": datetime.now().isoformat(),
        "config": {"languages": languages, "elements": elements, "methods": methods,
                   "model": args.model, "repeat": args.repeat, "max_repair": MAX_REPAIR_ITERATIONS},
        "results": all_results,
        "overall": {"total": total, "compile": total_compile, "run": total_run,
                    "numerical": total_num, "pass": total_pass},
    }, indent=2, default=str), encoding='utf-8')
    print(f"\nSaved to: {summary_file}")


if __name__ == "__main__":
    main()
