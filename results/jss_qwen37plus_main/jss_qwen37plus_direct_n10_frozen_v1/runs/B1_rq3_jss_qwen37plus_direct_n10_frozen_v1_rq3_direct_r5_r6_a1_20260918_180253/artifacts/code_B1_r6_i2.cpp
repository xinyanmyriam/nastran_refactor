#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

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
    
    // Build the 12x12 local stiffness matrix for CROD element
    // DOF order: ux, uy, uz, rx, ry, rz for node A (indices 0-5), then same for node B (indices 6-11)
    // The stiffness matrix has the form:
    // [ Kaa  Kab ]
    // [ Kba  Kbb ]
    // where each block is 6x6.
    //
    // For a rod element:
    // - Axial part (affects translational DOFs 0-2 and 6-8): 
    //   Kaa_axial =  dscl * [N  0; 0  0] (6x6)
    //   Kbb_axial =  dscl * [N  0; 0  0] (6x6)  
    //   Kab_axial = -dscl * [N  0; 0  0] (6x6)
    //   Kba_axial = Kab_axial^T
    //
    // - Torsional part (affects rotational DOFs 3-5 and 9-11):
    //   Kaa_torsion =  dscr * [0  0; 0  N] (6x6)
    //   Kbb_torsion =  dscr * [0  0; 0  N] (6x6)
    //   Kab_torsion = -dscr * [0  0; 0  N] (6x6)
    //   Kba_torsion = Kab_torsion^T
    //
    // So total:
    // Kaa = Kaa_axial + Kaa_torsion
    // Kbb = Kbb_axial + Kbb_torsion
    // Kab = Kab_axial + Kab_torsion
    // Kba = Kba_axial + Kba_torsion
    
    // Initialize 12x12 stiffness matrix
    Eigen::Matrix<double, 12, 12> K;
    K.setZero();
    
    // Build 6x6 blocks
    Eigen::Matrix<double, 6, 6> Kaa, Kbb, Kab, Kba;
    Kaa.setZero();
    Kbb.setZero();
    Kab.setZero();
    Kba.setZero();
    
    // Axial part: affects DOFs 0-2 (ux,uy,uz) for each node
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Kaa(i, j) =  dscl * N(i, j);
            Kbb(i, j) =  dscl * N(i, j);
            Kab(i, j) = -dscl * N(i, j);
            Kba(i, j) = -dscl * N(i, j); // since symmetric, Kba = Kab^T, and N is symmetric
        }
    }
    
    // Torsional part: affects DOFs 3-5 (rx,ry,rz) for each node
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Kaa(i+3, j+3) =  dscr * N(i, j);
            Kbb(i+3, j+3) =  dscr * N(i, j);
            Kab(i+3, j+3) = -dscr * N(i, j);
            Kba(i+3, j+3) = -dscr * N(i, j);
        }
    }
    
    // Assemble 12x12 matrix
    // Top-left 6x6: Kaa
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K(i, j) = Kaa(i, j);
        }
    }
    
    // Top-right 6x6: Kab
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K(i, j+6) = Kab(i, j);
        }
    }
    
    // Bottom-left 6x6: Kba
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K(i+6, j) = Kba(i, j);
        }
    }
    
    // Bottom-right 6x6: Kbb
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K(i+6, j+6) = Kbb(i, j);
        }
    }
    
    // Output as JSON string manually
    std::ostringstream oss;
    oss << std::setprecision(15);
    oss << "{\n  \"stiffness_matrix\": [\n";

    for (int i = 0; i < 12; ++i) {
        oss << "    [";
        for (int j = 0; j < 12; ++j) {
            oss << K(i, j);
            if (j < 11) oss << ", ";
        }
        oss << "]";
        if (i < 11) oss << ",";
        oss << "\n";
    }

    oss << "  ]\n}";

    std::cout << oss.str() << std::endl;
    
    return 0;
}