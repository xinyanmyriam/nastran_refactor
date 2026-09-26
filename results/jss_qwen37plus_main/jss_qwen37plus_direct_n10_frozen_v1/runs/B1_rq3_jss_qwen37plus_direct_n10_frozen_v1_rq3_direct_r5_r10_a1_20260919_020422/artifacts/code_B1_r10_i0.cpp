#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Using double precision throughout
using Real = double;
using Matrix6d = Eigen::Matrix<Real, 6, 6>;
using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix<Real, 3, 3>;

// Mock MAT subroutine: sets E and G for given material ID
// In real NASTRAN, this would look up in MAT1 table
void MAT(int /*matid*/, Real& E, Real& G, Real& /*nu*/, Real& /*rho*/, Real& /*alpha*/) {
    // For test case: E = 200e9, G = 76.923e9 (implied nu ≈ 0.3)
    E = 200.0e9;
    G = 76.923e9;
}

// Mock TRANSD: returns identity for basic coordinate system (CID=0)
// For non-basic, we'd transform, but test uses basic coords only
void TRANSD(int cid, Matrix3d& T) {
    if (cid == 0) {
        T.setIdentity();
    } else {
        // For test case, we assume basic coordinates only
        T.setIdentity();
    }
}

// Mock GMMATD: matrix multiplication C = A * B (with transpose flags)
// Here: GMMATD(A, arows, acols, at, B, brows, bcols, bt, C)
// at/bt: 0 = no transpose, 1 = transpose
void GMMATD(const Matrix3d& A, int /*arows*/, int /*acols*/, int at,
            const Matrix3d& B, int /*brows*/, int /*bcols*/, int bt,
            Matrix3d& C) {
    Matrix3d A_use = (at == 0) ? A : A.transpose();
    Matrix3d B_use = (bt == 0) ? B : B.transpose();
    C = A_use * B_use;
}

// Mock SMA1B: inserts KE into global stiffness matrix
// We'll just accumulate local K into our result
void SMA1B(Matrix6d& KE, int /*nonpvt*/, int /*npvt*/, int /*ifkgg*/, Real /*scale*/) {
    // We'll handle assembly in main logic; this is just a placeholder
}

// Main rod element stiffness computation
Matrix6d compute_crod_stiffness(
    const Vector3d& nodeA,  // coordinates of node A (basic coord)
    const Vector3d& nodeB,  // coordinates of node B (basic coord)
    Real E,                 // Young's modulus
    Real A,                 // cross-sectional area
    Real G,                 // shear modulus
    Real J                   // polar moment of inertia
) {
    // Compute length
    Vector3d vec = nodeB - nodeA;
    Real L = vec.norm();
    
    // Handle zero-length
    if (L == 0.0) {
        throw std::runtime_error("Zero-length rod element");
    }
    
    // Direction cosines (normalized vector)
    Vector3d n = vec / L;
    
    // Axial and torsional stiffness coefficients
    Real k_axial = E * A / L;   // axial stiffness
    Real k_torsion = G * J / L; // torsional stiffness
    
    // Build the 6x6 local stiffness matrix in basic coordinates
    // DOF order: ux, uy, uz, rx, ry, rz
    Matrix6d K = Matrix6d::Zero();
    
    // Axial part (affects ux,uy,uz at both ends)
    // [ n_x^2   n_x*n_y   n_x*n_z   -n_x^2   -n_x*n_y   -n_x*n_z ]
    // [ n_x*n_y n_y^2     n_y*n_z   -n_x*n_y -n_y^2     -n_y*n_z ]
    // [ n_x*n_z n_y*n_z   n_z^2     -n_x*n_z -n_y*n_z   -n_z^2   ]
    // [ -n_x^2  -n_x*n_y  -n_x*n_z  n_x^2    n_x*n_y    n_x*n_z  ]
    // [ -n_x*n_y -n_y^2   -n_y*n_z  n_x*n_y  n_y^2      n_y*n_z  ]
    // [ -n_x*n_z -n_y*n_z  -n_z^2    n_x*n_z  n_y*n_z    n_z^2    ]
    
    // Torsional part (affects rx,ry,rz at both ends)
    // Only diagonal terms for rx-rx, ry-ry, rz-rz coupling with sign
    // [ 0 0 0 0 0 0 ]
    // [ 0 0 0 0 0 0 ]
    // [ 0 0 0 0 0 0 ]
    // [ 0 0 0 k_torsion*nx^2 k_torsion*nx*ny k_torsion*nx*nz ]
    // [ 0 0 0 k_torsion*nx*ny k_torsion*ny^2 k_torsion*ny*nz ]
    // [ 0 0 0 k_torsion*nx*nz k_torsion*ny*nz k_torsion*nz^2 ]
    // and symmetric negative block at bottom-left
    
    // But standard rod element has:
    // Axial: connects translational DOFs (1-3 and 7-9 in 12x12, but here 6x6 local)
    // Torsional: connects rotational DOFs (4-6 and 10-12 in 12x12, but here 4-6 and 1-3? Wait...)
    
    // Actually, standard CROD has:
    // - Axial stiffness between ux_A and ux_B (and coupled uy, uz via direction cosines)
    // - Torsional stiffness between rx_A and rx_B, ry_A and ry_B, rz_A and rz_B, also via direction cosines
    
    // So full 6x6 local stiffness (in basic coordinates) is:
    // [ k_axial * n*n^T     -k_axial * n*n^T ]
    // [ -k_axial * n*n^T     k_axial * n*n^T ]
    // plus
    // [ 0                    0                ]
    // [ 0    k_torsion * n*n^T               ]
    
    // But wait — torsion couples rotations, not translations.
    // Correct CROD stiffness (per NASTRAN theory) is:
    // K = [ k_axial * n*n^T         0                  -k_axial * n*n^T         0               ]
    //     [ 0                  k_torsion * n*n^T       0                  -k_torsion * n*n^T ]
    //     [ -k_axial * n*n^T     0                   k_axial * n*n^T          0               ]
    //     [ 0                  -k_torsion * n*n^T    0                   k_torsion * n*n^T ]
    
    // However, the Fortran code builds a 6x6 as:
    // Rows 1-3: axial terms (DSCL * N matrix)
    // Rows 4-6: torsional terms (DSCR * N matrix) but placed at rows 22-24 etc. — that's 12x12 indexing
    
    // Since we need 6x6 for one "end pair", and the problem asks for 6x6 stiffness matrix,
    // and the test expects K(1,1) = E*A/L = 1e9, that's the axial term for ux-ux.
    
    // The standard local stiffness for CROD (with x-axis aligned) is:
    // [ k_axial  0    0    0    0    0 ]
    // [ 0        0    0    0    0    0 ]
    // [ 0        0    0    0    0    0 ]
    // [ 0        0    0    k_torsion 0  0 ]
    // [ 0        0    0    0         0  0 ]
    // [ 0        0    0    0         0  0 ]
    // and then transformed.
    
    // But the Fortran builds the full directional version.
    // From Fortran: it computes N = n * n^T (3x3), then places DSCL*N at top-left, 
    // -DSCL*N at top-right, -DSCL*N at bottom-left, DSCL*N at bottom-right for axial,
    // and similarly DSCR*N for rotational blocks.
    
    // So the 6x6 local stiffness matrix (in basic coordinates) is:
    //
    // [ DSCL*N     0         -DSCL*N     0       ]
    // [ 0          DSCR*N    0          -DSCR*N  ]
    // [ -DSCL*N    0          DSCL*N     0       ]
    // [ 0         -DSCR*N     0          DSCR*N   ]
    //
    // But arranged as 6x6 with DOF order: ux,uy,uz,rx,ry,rz for node A, then same for node B.
    // So indices: 
    // 0,1,2 = ux,uy,uz of A
    // 3,4,5 = rx,ry,rz of A
    // 6,7,8 = ux,uy,uz of B
    // 9,10,11 = rx,ry,rz of B
    //
    // Therefore, the 6x6 we return is actually the upper-left 6x6 of the 12x12? 
    // No — the problem says "6x6 stiffness matrix" and "Output: 12x12 stiffness matrix (6 DOF per node)".
    // But then it says "Print the 6x6 stiffness matrix as JSON".
    //
    // Re-reading: "Output specification: Print the 6x6 stiffness matrix as JSON"
    // and "Expected K(1,1) = E*A/L = 1e9 N/m"
    //
    // In standard FEM, the element stiffness is 12x12. But the problem says "6x6 stiffness matrix",
    // and the Fortran computes two 6x6 matrices: K(NPVT,NPVT) and K(NPVT,J). 
    // So likely they want the 6x6 "coupling" matrix between the two nodes, i.e., the off-diagonal 6x6 block.
    //
    // However, the test says "Expected K(1,1) = 1e9", which is the axial stiffness, and in the 12x12,
    // K(1,1) is the first diagonal term (ux_A - ux_A), which is +k_axial * n_x^2, but if aligned with x-axis, n_x=1, so = k_axial.
    //
    // So they want the full 12x12? But output spec says "6x6 stiffness matrix".
    //
    // Let's check the Fortran: it computes a 6x6 matrix KE(36) and calls SMA1B to insert into global matrix.
    // And the comment says "TWO 6X6 MATRICES K(NPVT,NPVT) AND K(NPVT,J)" — meaning the diagonal block (6x6 for node NPVT) and off-diagonal (6x6 for NPVT-J).
    //
    // But the problem says: "CROD - Rod element (axial + torsion), 2 nodes, 3 DOF/node, 6x6 stiffness matrix"
    // That's ambiguous. However, standard practice is to return the 12x12 element stiffness.
    // But the output spec says "6x6 stiffness matrix" and the JSON has 6 rows.
    //
    // Looking again at the problem statement: "Output specification: Print the 6x6 stiffness matrix as JSON"
    // and the test case expects K(1,1)=1e9. In a 6x6, (1,1) would be the first DOF of the first node vs first DOF — which is ux_A vs ux_A, so diagonal term.
    //
    // However, a 6x6 cannot represent a 2-node element with 6 DOF per node — that requires 12x12.
    //
    // I think there's a misunderstanding. Let me re-read: "6x6 stiffness matrix" — but then it says "Output: 12x12 stiffness matrix (6 DOF per node)".
    // That's contradictory.
    //
    // The key is the Fortran: it computes two 6x6 matrices: K(NPVT,NPVT) and K(NPVT,J). 
    // So for node NPVT (say node 1), K(NPVT,NPVT) is the 6x6 diagonal block (DOFs of node 1 coupled to themselves),
    // and K(NPVT,J) is the 6x6 off-diagonal block (DOFs of node 1 coupled to DOFs of node 2).
    //
    // But the problem says "CROD ... 6x6 stiffness matrix". 
    // And the test expects K(1,1)=1e9, which matches the diagonal term of the diagonal block.
    //
    // However, the diagonal block for a rod has:
    // - For ux_A: k_axial * n_x^2 (but also coupled to uy_A, uz_A)
    // - For rx_A: k_torsion * n_x^2 (coupled to ry_A, rz_A)
    //
    // So the 6x6 diagonal block is:
    // [ k_axial*n_x^2  k_axial*n_x*n_y  k_axial*n_x*n_z  0  0  0 ]
    // [ ... ]
    // [ 0 0 0 k_torsion*n_x^2 k_torsion*n_x*n_y k_torsion*n_x*n_z ]
    // etc.
    //
    // But the problem says "Expected K(1,1) = E*A/L = 1e9", which is exactly k_axial when n_x=1.
    //
    // Given the test case: Node A (0,0,0), Node B (2,0,0), so n = (1,0,0).
    // So n_x=1, n_y=0, n_z=0.
    //
    // Then the 6x6 diagonal block (for node A) is:
    // [ k_axial  0  0  0  0  0 ]
    // [ 0        0  0  0  0  0 ]
    // [ 0        0  0  0  0  0 ]
    // [ 0        0  0  k_torsion 0  0 ]
    // [ 0        0  0  0         0  0 ]
    // [ 0        0  0  0         0  0 ]
    //
    // But that doesn't match standard — the diagonal block should have positive terms only on diagonal for the self-coupling.
    //
    // Alternatively, perhaps they want the 6x6 "element stiffness" in local coordinates (aligned with rod), which is diagonal:
    // [ k_axial, 0, 0, 0, 0, 0,
    //   0, 0, 0, 0, 0, 0,
    //   0, 0, 0, 0, 0, 0,
    //   0, 0, 0, k_torsion, 0, 0,
    //   0, 0, 0, 0, 0, 0,
    //   0, 0, 0, 0, 0, 0 ]
    // but that's 6x6 with zeros — not right.
    //
    // Let me check standard CROD theory: the element stiffness matrix is 12x12. The upper-left 6x6 block is:
    // [ k_axial*n_x^2    k_axial*n_x*n_y    k_axial*n_x*n_z    0    0    0 ]
    // [ k_axial*n_x*n_y  k_axial*n_y^2      k_axial*n_y*n_z    0    0    0 ]
    // [ k_axial*n_x*n_z  k_axial*n_y*n_z    k_axial*n_z^2      0    0    0 ]
    // [ 0                0                  0                  k_torsion*n_x^2  k_torsion*n_x*n_y  k_torsion*n_x*n_z ]
    // [ 0                0                  0                  k_torsion*n_x*n_y  k_torsion*n_y^2    k_torsion*n_y*n_z ]
    // [ 0                0                  0                  k_torsion*n_x*n_z  k_torsion*n_y*n_z  k_torsion*n_z^2 ]
    //
    // And the upper-right 6x6 block is the negative of that.
    //
    // So the 6x6 matrix they want is likely the upper-left 6x6 block (the "K(NPVT,NPVT)" from Fortran).
    //
    // In the test case, n=(1,0,0), so:
    // Row 0: [k_axial, 0, 0, 0, 0, 0]
    // Row 1: [0, 0, 0, 0, 0, 0]
    // Row 2: [0, 0, 0, 0, 0, 0]
    // Row 3: [0, 0, 0, k_torsion, 0, 0]
    // Row 4: [0, 0, 0, 0, 0, 0]
    // Row 5: [0, 0, 0, 0, 0, 0]
    //
    // But that gives K(1,1)=k_axial=1e9, which matches.
    //
    // So we'll compute the 6x6 diagonal block for node A (NPVT).
    
    // Build N = n * n^T (3x3)
    Matrix3d N = n * n.transpose();
    
    // Diagonal block for node A: 
    // Top-left 3x3: k_axial * N
    // Bottom-right 3x3: k_torsion * N
    Matrix6d K_diag = Matrix6d::Zero();
    K_diag.topLeftCorner(3,3) = k_axial * N;
    K_diag.bottomRightCorner(3,3) = k_torsion * N;
    
    return K_diag;
}

// Helper to format double in scientific notation with fixed precision
std::string format_double(Real x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros after decimal
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos+1);
    }
    return s;
}

int main() {
    try {
        // Test case data
        Vector3d nodeA(0.0, 0.0, 0.0);
        Vector3d nodeB(2.0, 0.0, 0.0);
        Real E = 200.0e9;      // Pa
        Real A = 0.01;         // m^2
        Real G = 76.923e9;     // Pa
        Real J = 5.0e-6;       // m^4
        
        // Compute stiffness matrix
        Matrix6d K = compute_crod_stiffness(nodeA, nodeB, E, A, G, J);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 6; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << "[";
            for (int j = 0; j < 6; ++j) {
                if (j > 0) std::cout << ",";
                std::cout << format_double(K(i,j));
            }
            std::cout << "]";
        }
        std::cout << "]}" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}