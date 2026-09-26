#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>

// For JSON output
using json = nlohmann::json;

// Mock MAT subroutine - returns E and G for given material ID
// In real NASTRAN, this would look up material properties
void MAT(int /*matidc*/, int /*matflg*/, double eltemp, double& E, double& G, double& /*nu*/, double& /*rho*/, double& /*alpha*/) {
    // For our test case: steel-like material
    // E = 200e9 Pa, G = 76.923e9 Pa (so nu = 0.3)
    E = 200.0e9;
    G = 76.923e9;
}

// Mock TRANSD - returns identity matrix for basic coordinate system (ID=0)
// For non-basic systems, we'd need the transformation, but our test uses basic
void TRANSD(int coord_id, Eigen::Matrix<double, 3, 3>& T) {
    if (coord_id == 0) {
        // Basic coordinate system: identity matrix
        T.setIdentity();
    } else {
        // For simplicity in test case, use identity (all nodes in basic coords)
        T.setIdentity();
    }
}

// Matrix multiplication: C = A * B, where A is m x k, B is k x n, C is m x n
void GMMATD(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B, Eigen::MatrixXd& C) {
    C = A * B;
}

// SMA1B - assemble stiffness contribution to global matrix
// For our test, we just compute the local stiffness matrix
void SMA1B(const Eigen::MatrixXd& ke_local, int /*npvt*/, int /*nonpvt*/, 
           int /*ifkgg*/, double /*dampc*/) {
    // This is a placeholder - in real code it would assemble into global matrix
}

int main() {
    // Test case parameters
    const double nodeA[3] = {0.0, 0.0, 0.0};  // Node A coordinates
    const double nodeB[3] = {2.0, 0.0, 0.0};  // Node B coordinates
    
    const double E = 200.0e9;      // Young's modulus (Pa)
    const double A = 0.01;         // Cross-sectional area (m^2)
    const double G = 76.923e9;     // Shear modulus (Pa)
    const double J = 5.0e-6;       // Polar moment of inertia (m^4)
    
    // Compute length
    double dx = nodeB[0] - nodeA[0];
    double dy = nodeB[1] - nodeA[1];
    double dz = nodeB[2] - nodeA[2];
    double L = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Direction cosines (normalized direction vector)
    double xn[3] = {dx/L, dy/L, dz/L};
    
    // Compute stiffness coefficients
    double dscl = A * E / L;  // Axial stiffness coefficient
    double dscr = J * G / L;  // Torsional stiffness coefficient
    
    // Build the 3x3 "N" matrix (outer product of direction vector)
    // N = [xn[0]^2, xn[0]*xn[1], xn[0]*xn[2]]
    //     [xn[0]*xn[1], xn[1]^2, xn[1]*xn[2]]
    //     [xn[0]*xn[2], xn[1]*xn[2], xn[2]^2]
    Eigen::Matrix3d N;
    N << xn[0]*xn[0], xn[0]*xn[1], xn[0]*xn[2],
         xn[0]*xn[1], xn[1]*xn[1], xn[1]*xn[2],
         xn[0]*xn[2], xn[1]*xn[2], xn[2]*xn[2];
    
    // Since both nodes are in basic coordinate system (coord_id = 0),
    // no coordinate transformation needed - T = I, so T'NT = N
    
    // Build the 6x6 local stiffness matrix for CROD element
    // DOF order: ux, uy, uz, rx, ry, rz for each node (12 total, but we build 6x6 block)
    // The stiffness matrix has the form:
    // [ Kaa  Kab ]
    // [ Kba  Kbb ]
    // where each block is 6x6, but due to symmetry and rod assumptions:
    // Kaa = [ dscl*N   0     ]
    //       [   0    dscr*N ]
    // Kab = [-dscl*N   0     ]
    //       [   0   -dscr*N ]
    // Kba = Kab^T, Kbb = Kaa
    
    // Initialize 12x12 stiffness matrix (but we only need the 6x6 blocks)
    // For output, we need the 6x6 matrix that represents the element stiffness
    // in local coordinates (which is what the Fortran computes for K(NPVT,NPVT) and K(NPVT,J))
    
    // The Fortran code computes two 6x6 matrices: K(NPVT,NPVT) and K(NPVT,J)
    // But for the JSON output requirement, we need the full 6x6 stiffness matrix
    // that relates the 6 DOFs of node A to the 6 DOFs of node B.
    // However, the problem states "Output: 12x12 stiffness matrix (6 DOF per node)"
    // and the JSON should contain the 6x6 matrix.
    // Looking at the Fortran, it builds KE(36) which is a 6x6 matrix.
    // From the assembly logic, the 6x6 matrix computed is for the coupling between
    // the pivot node (NPVT) and the non-pivot node (J), i.e., the off-diagonal block.
    
    // But the problem says "Expected K(1,1) = E*A/L = 1e9 N/m", which is the axial stiffness
    // and appears in the diagonal of the 6x6 matrix for the axial DOFs.
    
    // Actually, the standard rod element stiffness in local coordinates is:
    // [  dscl  0     0     0     0     0   ]
    // [  0     0     0     0     0     0   ]
    // [  0     0     0     0     0     0   ]
    // [  0     0     0    dscr   0     0   ]
    // [  0     0     0     0     0     0   ]
    // [  0     0     0     0     0     0   ]
    // But this is not right for a 6x6 representation.
    
    // Let me re-read the Fortran: It builds KE(36) as a 6x6 matrix with:
    // KE(1), KE(2), KE(3), KE(7), KE(8), KE(9), KE(13), KE(14), KE(15) = DSCL * D(K2+i)
    // which are the first 3x3 block (axial part) at positions (1,1), (1,2), (1,3), (2,1), etc.
    // and KE(22), KE(23), KE(24), KE(28), KE(29), KE(30), KE(34), KE(35), KE(36) = DSCR * D(K2+i)
    // which are the second 3x3 block (torsional part) at positions (4,4), (4,5), (4,6), (5,4), etc.
    
    // So the 6x6 matrix has:
    // rows 1-3, cols 1-3: dscl * N
    // rows 4-6, cols 4-6: dscr * N
    // and the off-diagonal blocks are zero for a simple rod.
    
    // But the Fortran also does:
    // KE(1) = DSCL * D(K2)   -> (1,1)
    // KE(2) = DSCL * D(K2+1) -> (1,2)  
    // KE(3) = DSCL * D(K2+2) -> (1,3)
    // KE(7) = DSCL * D(K2+3) -> (2,1) since 7 = 1 + 6*1
    // KE(8) = DSCL * D(K2+4) -> (2,2)
    // KE(9) = DSCL * D(K2+5) -> (2,3)
    // KE(13) = DSCL * D(K2+6) -> (3,1) since 13 = 1 + 6*2
    // KE(14) = DSCL * D(K2+7) -> (3,2)
    // KE(15) = DSCL * D(K2+8) -> (3,3)
    // So the first 3x3 block is placed in rows 1-3, cols 1-3.
    
    // Similarly for torsion:
    // KE(22) = DSCR * D(K2) -> (4,4) since 22 = 1 + 6*3 + 1? Let's check indexing:
    // In Fortran, KE is 1-indexed, so KE(22) is row 4, col 4? 
    // Actually, KE is stored column-wise: KE(1) = (1,1), KE(2) = (2,1), ..., KE(6) = (6,1), KE(7) = (1,2), etc.
    // So KE(22) = 22nd element = row = (22-1)%6 + 1 = 5, col = (22-1)/6 + 1 = 4 -> (5,4)? 
    // But the Fortran comment says "KE(22) = DSCR * D(K2)" and then "KE(22) ... KE(36)" for the torsional part.
    // Looking at the pattern: KE(22), KE(23), KE(24) are consecutive, then KE(28), KE(29), KE(30), then KE(34), KE(35), KE(36).
    // This suggests they are placing in positions (4,4), (4,5), (4,6), (5,4), (5,5), (5,6), (6,4), (6,5), (6,6).
    // Because: 
    // KE(22) = (4,4) -> index = (4-1) + (4-1)*6 = 3 + 18 = 21 -> but Fortran is 1-indexed so 22
    // Yes! In Fortran column-major: element (i,j) is at position (i-1) + (j-1)*6 + 1
    // So (4,4): (4-1) + (4-1)*6 + 1 = 3 + 18 + 1 = 22 ✓
    // (4,5): 3 + 24 + 1 = 28 ✓
    // (4,6): 3 + 30 + 1 = 34 ✓
    // So the torsional part goes in rows 4-6, cols 4-6.
    
    // Therefore, the 6x6 matrix is:
    // [ dscl*N   0    ]
    // [   0    dscr*N ]
    // where N is 3x3, so the full matrix is 6x6.
    
    // Build the 6x6 local stiffness matrix
    Eigen::Matrix<double, 6, 6> K;
    K.setZero();
    
    // Axial part: rows 0-2, cols 0-2 (0-indexed)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            K(i, j) = dscl * N(i, j);
        }
    }
    
    // Torsional part: rows 3-5, cols 3-5
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            K(i+3, j+3) = dscr * N(i, j);
        }
    }
    
    // Verify expected value: K(1,1) = dscl * xn[0]^2
    // Since xn[0] = dx/L = 2.0/2.0 = 1.0, so K(0,0) = dscl * 1.0 = E*A/L = 200e9 * 0.01 / 2.0 = 1e9
    // So K(0,0) should be 1e9
    
    // Output as JSON
    json j;
    std::vector<std::vector<double>> stiffness_matrix;
    
    for (int i = 0; i < 6; ++i) {
        std::vector<double> row;
        for (int j = 0; j < 6; ++j) {
            row.push_back(K(i, j));
        }
        stiffness_matrix.push_back(row);
    }
    
    j["stiffness_matrix"] = stiffness_matrix;
    
    // Set precision for scientific notation
    std::cout << std::setprecision(15);
    std::cout << j.dump() << std::endl;
    
    return 0;
}