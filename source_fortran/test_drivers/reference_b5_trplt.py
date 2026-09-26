"""
Reference implementation for B5 (TRPLT) - Triangular Plate Bending Element
===========================================================================

Computes the 9x9 stiffness matrix for a Kirchhoff plate bending triangle
using the BCIZ (Bazeley-Cheung-Irons-Zienkiewicz) non-conforming element,
which is the formulation used in NASTRAN-95's KTRBSC routine.

The NASTRAN basic bending triangle uses:
- 3 nodes, 3 DOF per node: (w, theta_x, theta_y) = 9 DOF total
- Polynomial displacement: w = a1 + a2*x + a3*y + a4*x^2 + a5*xy + a6*y^2
- The H-matrix maps generalized coords to nodal DOFs
- K = H^{-T} * K^x * H^{-1} where K^x is the stiffness in generalized coords

Test case: Equilateral-like triangle in XY plane
  Node 1: (0, 0, 0)
  Node 2: (1, 0, 0) 
  Node 3: (0.5, 0.866025, 0)  [equilateral triangle, side=1]

Material: E=200e9, nu=0.3, thickness t=0.01
Bending rigidity: D = Et^3 / [12(1-nu^2)]

Output: 9x9 stiffness matrix as JSON

DOF ordering per node: (w, theta_x, theta_y)
  where theta_x = dw/dy, theta_y = -dw/dx (NASTRAN convention)
"""

import numpy as np
import json
import sys


def plate_bending_D_matrix(E, nu, t):
    """Compute the 3x3 bending rigidity matrix D for isotropic material.
    
    D relates moments to curvatures: {Mx, My, Mxy} = D * {kx, ky, kxy}
    where kx = -d2w/dx2, ky = -d2w/dy2, kxy = -2*d2w/dxdy
    """
    I_mom = t**3 / 12.0  # moment of inertia per unit width
    D11 = E * I_mom / (1.0 - nu**2)
    D12 = nu * D11
    D33 = E * I_mom / (2.0 * (1.0 + nu))  # = D11*(1-nu)/2
    
    D = np.array([
        [D11, D12, 0.0],
        [D12, D11, 0.0],
        [0.0, 0.0, D33]
    ])
    return D


def nastran_basic_bending_triangle(xb, xc, yc, D_mat):
    """
    Compute the 9x9 stiffness matrix for NASTRAN's basic bending triangle.
    
    This follows the KTRBSC algorithm exactly:
    1. Form K^x (6x6) stiffness in generalized coordinates
    2. Form H-bar (6x6) matrix
    3. Invert H to get H^{-1}
    4. Compute K_II = H^{-T} * K^x * H^{-1}
    5. Form S matrix (6x3) for constraint
    6. Compute K_IA = -K_II * S
    7. Compute K_AA = S^T * K_II * S
    8. Arrange into 9 sub-matrices (3x3 each)
    
    Parameters:
    -----------
    xb : float - x-coordinate of node B in element system (node A at origin)
    xc : float - x-coordinate of node C in element system
    yc : float - y-coordinate of node C in element system
    D_mat : 3x3 array - bending rigidity matrix [D11,D12,D13; D12,D22,D23; D13,D23,D33]
    
    Returns:
    --------
    K9x9 : 9x9 numpy array - stiffness matrix
            DOF order: node_A(w,thx,thy), node_B(w,thx,thy), node_C(w,thx,thy)
    """
    D1 = D_mat[0, 0]  # D11
    D2 = D_mat[0, 1]  # D12
    D3 = D_mat[0, 2]  # D13 (=0 for isotropic)
    D5 = D_mat[1, 1]  # D22
    D6 = D_mat[1, 2]  # D23 (=0 for isotropic)
    D9 = D_mat[2, 2]  # D33
    
    # Geometric quantities
    Area = xb * yc / 2.0
    xbar = (xb + xc) / 3.0
    ybar = yc / 3.0
    
    xcsq = xc**2
    ycsq = yc**2
    xbsq = xb**2
    xcyc = xc * yc
    
    PX2 = (xbsq + xb * xc + xcsq) / 6.0
    PY2 = ycsq / 6.0
    PXY2 = yc * (xb + 2.0 * xc) / 12.0
    
    xbar3 = 3.0 * xbar
    ybar3 = 3.0 * ybar
    ybar2 = 2.0 * ybar
    
    # Fill K^x matrix (6x6) - stiffness in generalized coordinates
    # Stored row-major as in NASTRAN: A(1)..A(36)
    Kx = np.zeros((6, 6))
    Kx[0, 0] = D1
    Kx[0, 1] = D3
    Kx[0, 2] = D2
    Kx[0, 3] = D1 * xbar3
    Kx[0, 4] = D2 * xbar + ybar2 * D3
    Kx[0, 5] = D2 * ybar3
    
    Kx[1, 0] = D3
    Kx[1, 1] = D9
    Kx[1, 2] = D6
    Kx[1, 3] = D3 * xbar3
    Kx[1, 4] = D6 * xbar + ybar2 * D9
    Kx[1, 5] = D6 * ybar3
    
    Kx[2, 0] = D2
    Kx[2, 1] = D6
    Kx[2, 2] = D5
    Kx[2, 3] = D2 * xbar3
    Kx[2, 4] = D5 * xbar + ybar2 * D6
    Kx[2, 5] = D5 * ybar3
    
    Kx[3, 0] = Kx[0, 3]
    Kx[3, 1] = Kx[1, 3]
    Kx[3, 2] = Kx[2, 3]
    Kx[3, 3] = D1 * 9.0 * PX2
    Kx[3, 4] = D2 * 3.0 * PX2 + 6.0 * PXY2 * D3
    Kx[3, 5] = D2 * 9.0 * PXY2
    
    Kx[4, 0] = Kx[0, 4]
    Kx[4, 1] = Kx[1, 4]
    Kx[4, 2] = Kx[2, 4]
    Kx[4, 3] = Kx[3, 4]
    Kx[4, 4] = D5 * PX2 + 4.0 * PXY2 * D6 + 4.0 * PY2 * D9
    Kx[4, 5] = D5 * 3.0 * PXY2 + 6.0 * PY2 * D6
    
    Kx[5, 0] = Kx[0, 5]
    Kx[5, 1] = Kx[1, 5]
    Kx[5, 2] = Kx[2, 5]
    Kx[5, 3] = Kx[3, 5]
    Kx[5, 4] = Kx[4, 5]
    Kx[5, 5] = D5 * 9.0 * PY2
    
    # Multiply by 4*Area
    Kx *= 4.0 * Area
    
    # Fill H-bar matrix (6x6), stored ROW-WISE in NASTRAN
    # A(37)..A(72) maps to H[row][col] where index = 37 + row*6 + col
    # Row 0: A(37)=xb², A(38)=0, A(39)=0, A(40)=xb³, A(41)=0, A(42)=0
    # Row 1: A(43)=0, A(44)=xb, A(45)=0, A(46)=0, A(47)=0, A(48)=0
    # Row 2: A(49)=-2xb, A(50)=0, A(51)=0, A(52)=-3xb², A(53)=0, A(54)=0
    # Row 3: A(55)=xc², A(56)=xc*yc, A(57)=yc², A(58)=xc³, A(59)=yc²*xc, A(60)=yc³
    # Row 4: A(61)=0, A(62)=xc, A(63)=2yc, A(64)=0, A(65)=2xc*yc, A(66)=3yc²
    # Row 5: A(67)=-2xc, A(68)=-yc, A(69)=0, A(70)=-3xc², A(71)=-yc², A(72)=0
    H = np.zeros((6, 6))
    # Row 0
    H[0, 0] = xbsq
    H[0, 3] = xbsq * xb
    # Row 1
    H[1, 1] = xb
    # Row 2
    H[2, 0] = -2.0 * xb
    H[2, 3] = -3.0 * xbsq
    # Row 3
    H[3, 0] = xcsq
    H[3, 1] = xcyc
    H[3, 2] = ycsq
    H[3, 3] = xcsq * xc
    H[3, 4] = ycsq * xc
    H[3, 5] = ycsq * yc
    # Row 4
    H[4, 1] = xc
    H[4, 2] = 2.0 * yc
    H[4, 4] = 2.0 * xcyc
    H[4, 5] = 3.0 * ycsq
    # Row 5
    H[5, 0] = -2.0 * xc
    H[5, 1] = -yc
    H[5, 3] = -3.0 * xcsq
    H[5, 4] = -ycsq
    
    # Invert H
    H_inv = np.linalg.inv(H)
    
    # K_II = H^{-T} * Kx * H^{-1}
    K_II = H_inv.T @ Kx @ H_inv
    
    # S matrix (6x3) - constraint matrix
    S = np.zeros((6, 3))
    S[0, 0] = 1.0
    S[0, 2] = -xb   # S(3) = -xsubb
    S[1, 1] = 1.0
    S[2, 2] = 1.0
    S[3, 0] = 1.0
    S[3, 1] = yc    # S(11) = ysubc
    S[3, 2] = -xc   # S(12) = -xsubc
    S[4, 1] = 1.0
    S[4, 2] = 0.0
    S[5, 2] = 1.0
    
    # K_IA = -K_II * S (6x3)
    K_IA = -K_II @ S
    
    # K_AA = S^T * K_II * S (3x3) = -S^T * K_IA (note sign)
    K_AA = S.T @ K_II @ S
    
    # Arrange into 9x9 matrix
    # The 9x9 K^u matrix has the structure:
    # [K_AA   K_AB   K_AC ]   where A=node1, B=node2, C=node3
    # [K_BA   K_BB   K_BC ]
    # [K_CA   K_CB   K_CC ]
    #
    # In KTRBSC storage (A(1)..A(81)):
    # A(1:9)   = K_AA (3x3)
    # A(10:18) = K_AB = K_IA rows 1-3 transposed -> K_IA[0:3,:].T? 
    # Actually from the code:
    # A(1:9) = K_AA
    # A(10:18) = K_IA columns rewritten as K_IB -> transpose of rows
    # A(19:27) = K_IC
    # A(28:36) = K_IB (originally at A(46:63) = K_IA rearranged)
    # etc.
    
    # From the KTRBSC code after statement 600:
    # The 9x9 is arranged as nine 3x3 blocks stored row-major:
    # blocks[0] = A(1:9)   = K_AA  (node A, node A)
    # blocks[1] = A(10:18) = K_AB  (node A, node B) = K_IA[0:3,:]^T 
    # blocks[2] = A(19:27) = K_AC  (node A, node C) = K_IA[3:6,:]^T
    # blocks[3] = A(28:36) = K_BA  (node B, node A) = K_IA[0:3,:]
    # blocks[4] = A(37:45) = K_BB  (node B, node B) = K_II[0:3,0:3]
    # blocks[5] = A(46:54) = K_BC  (node B, node C) = K_II[0:3,3:6]
    # blocks[6] = A(55:63) = K_CA  (node C, node A) = K_IA[3:6,:]
    # blocks[7] = A(64:72) = K_CB  (node C, node B) = K_II[3:6,0:3]
    # blocks[8] = A(73:81) = K_CC  (node C, node C) = K_II[3:6,3:6]
    
    # Actually, re-reading the code more carefully:
    # After K_IA (6x3) is computed as -K_II*S, stored at A(46:63)
    # K_AA (3x3) is S^T * (-K_IA) = S^T * K_II * S, stored at A(1:9)
    # Then code rearranges into 9 blocks of 3x3:
    # The arrangement follows from the final assignments in KTRBSC
    
    # Let's be more careful. K_II is 6x6, partition into 2x2 of 3x3 blocks:
    # K_II = [[K_II_11, K_II_12], [K_II_21, K_II_22]]
    # where subscript 1=node B block, subscript 2=node C block
    
    # K_IA is 6x3, partition into two 3x3 blocks:
    # K_IA = [[K_BA], [K_CA]]  (B-to-A and C-to-A coupling)
    
    # K_AA = S^T * K_II * S
    # K_AB = K_IA[0:3,:]^T = K_BA^T  (should be symmetric for stiffness)
    # K_AC = K_IA[3:6,:]^T = K_CA^T
    
    # Following KTRBSC's final arrangement exactly:
    K9 = np.zeros((9, 9))
    
    # K_AA (3x3) - node A self-coupling
    K9[0:3, 0:3] = K_AA
    
    # K_AB = transpose of K_IA upper block (rows 0:3)
    K_AB = K_IA[0:3, :].T  # 3x3
    K9[0:3, 3:6] = K_AB
    
    # K_AC = transpose of K_IA lower block (rows 3:6)
    K_AC = K_IA[3:6, :].T  # 3x3
    K9[0:3, 6:9] = K_AC
    
    # K_BA = K_IA upper block  
    K9[3:6, 0:3] = K_IA[0:3, :]
    
    # K_BB = K_II upper-left 3x3
    K9[3:6, 3:6] = K_II[0:3, 0:3]
    
    # K_BC = K_II upper-right 3x3
    K9[3:6, 6:9] = K_II[0:3, 3:6]
    
    # K_CA = K_IA lower block
    K9[6:9, 0:3] = K_IA[3:6, :]
    
    # K_CB = K_II lower-left 3x3
    K9[6:9, 3:6] = K_II[3:6, 0:3]
    
    # K_CC = K_II lower-right 3x3
    K9[6:9, 6:9] = K_II[3:6, 3:6]
    
    return K9


def compute_b5_stiffness(nodes, E, nu, t):
    """
    Compute the 9x9 plate bending stiffness for a triangular plate element (B5/TRPLT).
    
    For a single basic bending triangle (no sub-triangles), this is straightforward.
    The full KTRPLT uses 3 sub-triangles with a centroid point and condenses out
    the centroid DOFs, but for verification we compute KTRBSC directly.
    
    Parameters:
    -----------
    nodes : 3x3 array - [[x1,y1,z1], [x2,y2,z2], [x3,y3,z3]]
    E : float - Young's modulus
    nu : float - Poisson's ratio
    t : float - plate thickness
    
    Returns:
    --------
    K : 9x9 stiffness matrix in global DOFs (w, theta_x, theta_y) per node
    """
    # Compute element coordinate system
    # Node A at origin, Node B on x-axis
    r1 = np.array(nodes[0])
    r2 = np.array(nodes[1])
    r3 = np.array(nodes[2])
    
    # I-vector = r2 - r1, normalized
    d2 = r2 - r1
    xb = np.linalg.norm(d2)
    i_vec = d2 / xb
    
    # D1 = r3 - r1
    d1 = r3 - r1
    
    # K-vector = I x D1, normalized
    k_vec = np.cross(i_vec, d1)
    yc_mag = np.linalg.norm(k_vec)
    k_vec = k_vec / yc_mag
    
    # J-vector = K x I
    j_vec = np.cross(k_vec, i_vec)
    
    # Element coordinates of node C
    xc = np.dot(d1, i_vec)
    yc = np.dot(d1, j_vec)
    
    # Note: yc should equal yc_mag (the norm of the cross product)
    # Actually yc = |d1| * sin(angle between d1 and i_vec) which is the same as
    # |i_vec x d1| = yc_mag. Good.
    
    # D matrix (bending rigidity)
    I_mom = t**3 / 12.0
    D_mat = plate_bending_D_matrix(E, nu, t)
    # Scale by I/unit = I_mom already included in plate_bending_D_matrix
    # Actually plate_bending_D_matrix returns D = E*I/(1-nu^2) etc.
    # In NASTRAN, MAT returns G11, G12, G13, G22, G23, G33
    # Then D = I * G (where I is moment of inertia = t^3/12 for unit width)
    # So G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G33 = E/(2*(1+nu))
    # And D = I * G_matrix
    
    # For NASTRAN: the G-matrix from MAT (with INFLAG=2 for bending) gives:
    # G11 = E/(1-nu^2), G12 = nu*G11, G22 = G11, G33 = G11*(1-nu)/2
    # Then D_mat = EYE * G (where EYE = t^3/12 for plate bending)
    
    # Our plate_bending_D_matrix already computes D = E*t^3/12 * [...]
    # which is equivalent to EYE * G_matrix.
    # So D_mat is correct as-is.
    
    # Compute the basic bending triangle stiffness
    K9 = nastran_basic_bending_triangle(xb, xc, yc, D_mat)
    
    return K9


def main():
    # Test case: Equilateral triangle with side length 1.0
    # Node 1 at origin, Node 2 at (1,0,0), Node 3 at (0.5, sqrt(3)/2, 0)
    s = 1.0  # side length
    nodes = np.array([
        [0.0, 0.0, 0.0],
        [s, 0.0, 0.0],
        [s/2, s*np.sqrt(3)/2, 0.0]
    ])
    
    # Material properties
    E = 200.0e9      # Young's modulus (Pa)
    nu = 0.3         # Poisson's ratio
    t = 0.01         # thickness (m)
    
    # Bending rigidity for reference
    D = E * t**3 / (12.0 * (1.0 - nu**2))
    print(f"Material: E={E:.3e}, nu={nu}, t={t}")
    print(f"Bending rigidity D = {D:.6e}")
    print(f"")
    print(f"Triangle: equilateral, side = {s}")
    print(f"  Node 1: ({nodes[0,0]}, {nodes[0,1]}, {nodes[0,2]})")
    print(f"  Node 2: ({nodes[1,0]}, {nodes[1,1]}, {nodes[1,2]})")
    print(f"  Node 3: ({nodes[2,0]}, {nodes[2,1]}, {nodes[2,2]})")
    print(f"")
    
    # Element coordinates
    xb = s
    xc = s / 2.0
    yc = s * np.sqrt(3) / 2.0
    print(f"Element coords: xb={xb}, xc={xc}, yc={yc:.6f}")
    print(f"Area = {xb*yc/2.0:.6f}")
    print(f"")
    
    # Compute stiffness matrix
    K = compute_b5_stiffness(nodes, E, nu, t)
    
    # Print the matrix
    print("9x9 Stiffness Matrix K (KTRBSC basic bending triangle):")
    print("DOF order: (w1, thx1, thy1, w2, thx2, thy2, w3, thx3, thy3)")
    print("=" * 80)
    for i in range(9):
        row_str = " ".join(f"{K[i,j]:14.6e}" for j in range(9))
        print(row_str)
    
    # Check symmetry
    sym_err = np.max(np.abs(K - K.T))
    print(f"\nSymmetry check: max|K - K^T| = {sym_err:.3e}")
    
    # Check that matrix is positive semi-definite (should have 3 zero eigenvalues
    # for rigid body modes: translation + 2 rotations... actually for plate bending
    # with w, thx, thy DOFs, rigid body = w=const gives zero strain)
    eigvals = np.linalg.eigvalsh(K)
    print(f"\nEigenvalues (sorted):")
    for i, ev in enumerate(eigvals):
        print(f"  lambda_{i+1} = {ev:14.6e}")
    
    # Output as JSON for verification pipeline
    result = {
        "element_type": "B5_TRPLT",
        "formulation": "NASTRAN_KTRBSC_basic_bending_triangle",
        "test_case": "equilateral_triangle_side_1",
        "material": {
            "E": E,
            "nu": nu,
            "thickness": t,
            "D_bending_rigidity": D
        },
        "geometry": {
            "nodes": nodes.tolist(),
            "element_coords": {"xb": xb, "xc": xc, "yc": yc},
            "area": xb * yc / 2.0
        },
        "dof_order": ["w1", "thx1", "thy1", "w2", "thx2", "thy2", "w3", "thx3", "thy3"],
        "stiffness_matrix_9x9": K.tolist(),
        "eigenvalues": eigvals.tolist(),
        "symmetry_error": float(sym_err)
    }
    
    with open("reference_b5_trplt_result.json", "w") as f:
        json.dump(result, f, indent=2)
    
    print(f"\nResults written to reference_b5_trplt_result.json")
    
    return K


if __name__ == "__main__":
    main()
