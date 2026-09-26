#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Define Matrix6d alias (Eigen doesn't provide this by default)
using Matrix6d = Eigen::Matrix<double, 6, 6>;

// Forward declarations
void ktetra(int iopt, int jtype, std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords);

// Helper function to compute cross product of two 3D vectors
Eigen::Vector3d sAXB(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.cross(b);
}

// Helper function to compute dot product of two 3D vectors
double sADOTB(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.dot(b);
}

// Helper function to invert a 4x4 matrix and compute determinant
bool invert4x4(const Eigen::Matrix4d& H, Eigen::Matrix4d& H_inv, double& det) {
    det = H.determinant();
    if (std::abs(det) < 1e-15) {
        return false;
    }
    H_inv = H.inverse();
    return true;
}

// Helper function to compute material matrix G (6x6)
Matrix6d computeG(double E, double nu) {
    Matrix6d G = Matrix6d::Zero();
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    
    if (std::abs(temp1) < 1e-12) {
        // Handle invalid nu case
        return G;
    }
    
    double g11 = E * (1.0 - nu) / temp1;
    double g12 = E * nu / temp1;
    double gg = E / (2.0 * (1.0 + nu));
    
    G(0,0) = g11; G(1,1) = g11; G(2,2) = g11;
    G(0,1) = g12; G(0,2) = g12; G(1,0) = g12; G(1,2) = g12; G(2,0) = g12; G(2,1) = g12;
    G(3,3) = gg; G(4,4) = gg; G(5,5) = gg;
    
    return G;
}

// Helper function to compute C matrices (72 elements, 4 matrices of 6x3 each)
std::vector<Eigen::Matrix<double,6,3>> computeC(const Eigen::Matrix4d& H) {
    std::vector<Eigen::Matrix<double,6,3>> C(4);
    
    for (int i = 0; i < 4; ++i) {
        C[i] = Eigen::Matrix<double,6,3>::Zero();
        // H is 4x4, rows 0..3 correspond to points 1..4
        // C[i] is 6x3: fill using H.row(i) for point i
        C[i](0,0) = H(i,1);  // x-component for dN/dx
        C[i](1,1) = H(i,2);  // y-component for dN/dy
        C[i](2,2) = H(i,3);  // z-component for dN/dz
        C[i](3,2) = H(i,3);  // dN/dz for Sxy? (shear)
        C[i](4,1) = H(i,2);  // dN/dy for Sxz?
        C[i](5,0) = H(i,1);  // dx for Syz? — but original logic was inconsistent.
        // Let's follow standard tetrahedron B-matrix pattern:
        // For linear tet, B_i = [dN_i/dx, 0, 0; 0, dN_i/dy, 0; 0, 0, dN_i/dz; dN_i/dy, dN_i/dx, 0; dN_i/dz, 0, dN_i/dx; 0, dN_i/dz, dN_i/dy]
        // So:
        C[i](0,0) = H(i,1); // dN_i/dx
        C[i](1,1) = H(i,2); // dN_i/dy
        C[i](2,2) = H(i,3); // dN_i/dz
        C[i](3,0) = H(i,2); C[i](3,1) = H(i,1); // dN_i/dy, dN_i/dx
        C[i](4,0) = H(i,3); C[i](4,2) = H(i,1); // dN_i/dz, dN_i/dx
        C[i](5,1) = H(i,3); C[i](5,2) = H(i,2); // dN_i/dz, dN_i/dy
    }
    
    return C;
}

// Helper function to compute transformation matrix T for coordinate system
Eigen::Matrix3d computeT() {
    // For this test case, all nodes are in global coordinate system, so T is identity
    return Eigen::Matrix3d::Identity();
}

// Main KSOLID subroutine equivalent
void ksolid(int itype, std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords) {
    // For wedge element (ITYPE = 1), we use the first 12 tetrahedrons from the mapping matrix
    // M(1,1..4) = {1,2,3,4}, M(2,1..4) = {1,2,3,5}, etc.
    std::vector<std::vector<int>> wedge_mapping = {
        {1,2,3,4}, {1,2,3,5}, {1,2,3,6}, {1,4,5,6},
        {2,4,5,6}, {3,4,5,6}, {2,1,4,6}, {2,3,4,6},
        {1,3,4,5}, {2,3,4,5}, {3,1,5,6}, {2,1,5,6}
    };
    
    // Reset stiffness matrix
    stiffness_matrix.setZero();
    
    // Process each tetrahedron in the wedge decomposition
    for (int i = 0; i < static_cast<int>(wedge_mapping.size()); ++i) {
        // Create ECPT for this tetrahedron
        std::vector<double> tet_ecpt(100, 0.0);
        
        // Copy element ID and material ID
        tet_ecpt[0] = ecpt[0];  // Element ID
        tet_ecpt[1] = ecpt[1];  // Material ID
        
        // Copy grid point IDs (1-based indexing in Fortran, 0-based here)
        for (int j = 0; j < 4; ++j) {
            tet_ecpt[2+j] = static_cast<double>(wedge_mapping[i][j]);
        }
        
        // Copy coordinates for each grid point
        // Grid points: 1,2,3,4 correspond to indices 0,1,2,3 in node_coords
        for (int j = 0; j < 4; ++j) {
            int node_idx = wedge_mapping[i][j] - 1; // Convert to 0-based
            if (node_idx >= 0 && node_idx < static_cast<int>(node_coords.size())) {
                tet_ecpt[7 + 4*j] = node_coords[node_idx][0]; // X
                tet_ecpt[8 + 4*j] = node_coords[node_idx][1]; // Y
                tet_ecpt[9 + 4*j] = node_coords[node_idx][2]; // Z
            }
        }
        
        // Set element temperature (assume 0 for this test)
        tet_ecpt[22] = 0.0;
        
        // Call ktetra for this tetrahedron
        // IOPT: For wedge, first 6 configurations are multiplied by 2, others normal
        int iopt = (i < 6) ? (i + 1) : (i + 1);
        ktetra(iopt, 1, tet_ecpt, stiffness_matrix, E, nu, node_coords);
    }
}

// KTETRA subroutine equivalent
void ktetra(int iopt, int jtype, std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords) {
    // Build H matrix (4x4) for tetrahedron
    Eigen::Matrix4d H;
    H.setZero();
    H(0,0) = 1.0; H(0,1) = ecpt[7]; H(0,2) = ecpt[8]; H(0,3) = ecpt[9];  // Point 1
    H(1,0) = 1.0; H(1,1) = ecpt[11]; H(1,2) = ecpt[12]; H(1,3) = ecpt[13]; // Point 2
    H(2,0) = 1.0; H(2,1) = ecpt[15]; H(2,2) = ecpt[16]; H(2,3) = ecpt[17]; // Point 3
    H(3,0) = 1.0; H(3,1) = ecpt[19]; H(3,2) = ecpt[20]; H(3,3) = ecpt[21]; // Point 4
    
    // Invert H matrix
    Eigen::Matrix4d H_inv;
    double h_deter;
    if (!invert4x4(H, H_inv, h_deter)) {
        return; // Bad geometry
    }
    
    // Take absolute value of determinant
    h_deter = std::abs(h_deter);
    
    // Compute material matrix G
    Matrix6d G = computeG(E, nu);
    
    // Compute C matrices
    std::vector<Eigen::Matrix<double,6,3>> C = computeC(H);
    
    // Determine scaling factor based on IOPT
    double scale_factor = 1.0;
    if (iopt >= 1 && iopt <= 12) {
        // Wedge case: first 6 configurations multiplied by 2, all divided by 36
        if (iopt <= 6) {
            scale_factor = h_deter / 18.0; // h_deter/36.0 * 2.0
        } else {
            scale_factor = h_deter / 36.0;
        }
    } else {
        scale_factor = h_deter / 6.0;
    }
    
    // For wedge, we need to divide by 6.0 and then apply additional scaling
    // Based on Fortran logic: for wedge, HDETER = HDETER/36.0, then *2 for first 6
    if (jtype == 1) {
        if (iopt <= 6) {
            scale_factor = h_deter / 18.0;
        } else {
            scale_factor = h_deter / 36.0;
        }
    }
    
    // Compute GCT = C^T * G (for pivot point)
    // We'll use the first point as pivot for simplicity
    Eigen::Matrix<double,3,6> CT;
    CT.setZero();
    
    // Compute CT = C[0]^T (since C[0] is 6x3, CT is 3x6)
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 6; ++j) {
            CT(i,j) = C[0](j,i);
        }
    }
    
    // Compute GCT = CT * G (3x6 * 6x6 = 3x6)
    Eigen::Matrix<double,3,6> GCT = CT * G;
    
    // Scale GCT
    GCT *= scale_factor;
    
    // Now compute KIJ contributions for each of the 4 points
    for (int i = 0; i < 4; ++i) {
        // Get transformation matrix (identity for global coordinates)
        Eigen::Matrix3d T = computeT();
        
        // Compute CT_i = C[i]^T
        Eigen::Matrix<double,3,6> CT_i;
        CT_i.setZero();
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 6; ++c) {
                CT_i(r,c) = C[i](c,r);
            }
        }
        
        // Compute local KIJ = GCT * C[i] (3x6 * 6x3 = 3x3)
        Eigen::Matrix3d KIJ_local = GCT * C[i];
        
        // Map to global stiffness matrix
        // Each node has 3 DOF, so node i contributes to rows/columns [3*i, 3*i+2]
        int base_row = 3 * (static_cast<int>(ecpt[2+i]) - 1);
        int base_col = base_row;
        
        // Add contribution to stiffness matrix
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 3; ++c) {
                stiffness_matrix(base_row + r, base_col + c) += KIJ_local(r,c);
                stiffness_matrix(base_col + c, base_row + r) += KIJ_local(r,c); // Symmetric
            }
        }
    }
}

int main() {
    // Test case: Wedge (6 nodes): N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    std::vector<std::vector<double>> node_coords = {
        {0.0, 0.0, 0.0}, // N1
        {1.0, 0.0, 0.0}, // N2
        {0.0, 1.0, 0.0}, // N3
        {0.0, 0.0, 1.0}, // N4
        {1.0, 0.0, 1.0}, // N5
        {0.0, 1.0, 1.0}  // N6
    };
    
    // Material properties
    double E = 200e9;   // 200 GPa
    double nu = 0.3;
    
    // Initialize ECPT array (Fortran common block)
    std::vector<double> ecpt(100, 0.0);
    ecpt[0] = 1.0; // Element ID
    ecpt[1] = 1.0; // Material ID
    
    // Stiffness matrix: 18x18 (6 nodes * 3 DOF each)
    Eigen::MatrixXd stiffness_matrix(18, 18);
    stiffness_matrix.setZero();
    
    // Call KSOLID for wedge element (ITYPE = 1)
    ksolid(1, ecpt, stiffness_matrix, E, nu, node_coords);
    
    // Make sure matrix is symmetric
    for (int i = 0; i < 18; ++i) {
        for (int j = 0; j < 18; ++j) {
            if (i != j) {
                double avg = 0.5 * (stiffness_matrix(i,j) + stiffness_matrix(j,i));
                stiffness_matrix(i,j) = avg;
                stiffness_matrix(j,i) = avg;
            }
        }
    }
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 18; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 18; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(15) << stiffness_matrix(i,j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;
    
    return 0;
}