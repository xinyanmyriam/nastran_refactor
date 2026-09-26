#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Helper function to compute determinant of 4x4 matrix
double det4x4(const Eigen::Matrix<double, 4, 4>& M) {
    return M(0,0)*(M(1,1)*M(2,2)*M(3,3) - M(1,1)*M(2,3)*M(3,2) - M(1,2)*M(2,1)*M(3,3) + M(1,2)*M(2,3)*M(3,1) + M(1,3)*M(2,1)*M(3,2) - M(1,3)*M(2,2)*M(3,1))
         - M(0,1)*(M(1,0)*M(2,2)*M(3,3) - M(1,0)*M(2,3)*M(3,2) - M(1,2)*M(2,0)*M(3,3) + M(1,2)*M(2,3)*M(3,0) + M(1,3)*M(2,0)*M(3,2) - M(1,3)*M(2,2)*M(3,0))
         + M(0,2)*(M(1,0)*M(2,1)*M(3,3) - M(1,0)*M(2,3)*M(3,1) - M(1,1)*M(2,0)*M(3,3) + M(1,1)*M(2,3)*M(3,0) + M(1,3)*M(2,0)*M(3,1) - M(1,3)*M(2,1)*M(3,0))
         - M(0,3)*(M(1,0)*M(2,1)*M(3,2) - M(1,0)*M(2,2)*M(3,1) - M(1,1)*M(2,0)*M(3,2) + M(1,1)*M(2,2)*M(3,0) + M(1,2)*M(2,0)*M(3,1) - M(1,2)*M(2,1)*M(3,0));
}

// Helper function to invert 4x4 matrix (using Eigen)
bool invert4x4(const Eigen::Matrix<double, 4, 4>& M, Eigen::Matrix<double, 4, 4>& invM, double& det) {
    det = det4x4(M);
    if (std::abs(det) < 1e-15) {
        return false;
    }
    
    // Use Eigen's inverse for robustness
    invM = M.inverse();
    return true;
}

// Helper function to compute stiffness matrix for tetrahedral element
Eigen::Matrix<double, 12, 12> computeKTETRA(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    // Node coordinates: N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // Build H matrix (4x4): [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix<double, 4, 4> H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);
    
    // Compute determinant and inverse of H
    double hdet;
    Eigen::Matrix<double, 4, 4> Hinv;
    if (!invert4x4(H, Hinv, hdet)) {
        throw std::runtime_error("Singular H matrix in KTETRA");
    }
    
    // For tetrahedron, the volume is |det(H)|/6
    double volume = std::abs(hdet) / 6.0;
    
    // Material constants
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Invalid Poisson's ratio in KTETRA");
    }
    
    double G11 = E * (1.0 - nu) / temp1;
    double G12 = E * nu / temp1;
    double G44 = E / (2.0 * (1.0 + nu));
    
    // Build 6x6 material matrix G (stress-strain)
    Eigen::Matrix<double, 6, 6> G = Eigen::Matrix<double, 6, 6>::Zero();
    G(0,0) = G11; G(0,1) = G12; G(0,2) = G12;
    G(1,0) = G12; G(1,1) = G11; G(1,2) = G12;
    G(2,0) = G12; G(2,1) = G12; G(2,2) = G11;
    G(3,3) = G44;
    G(4,4) = G44;
    G(5,5) = G44;
    
    // Build C matrices (4 matrices of size 6x3)
    // C_i = [b_i 0 0; 0 b_i 0; 0 0 b_i] where b_i is the i-th row of Hinv (without first column)
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    
    // Extract b_i vectors (rows 0-3 of Hinv, columns 1-3)
    for (int i = 0; i < 4; ++i) {
        Eigen::Vector3d b_i;
        b_i << Hinv(i,1), Hinv(i,2), Hinv(i,3);
        
        // Build 6x3 matrix: [b_i; 0; 0] for first 3 rows, [0; b_i; 0] for middle 3, [0; 0; b_i] for last 3
        C[i] << b_i.transpose(), Eigen::Vector3d::Zero().transpose(), Eigen::Vector3d::Zero().transpose(),
                Eigen::Vector3d::Zero().transpose(), b_i.transpose(), Eigen::Vector3d::Zero().transpose(),
                Eigen::Vector3d::Zero().transpose(), Eigen::Vector3d::Zero().transpose(), b_i.transpose();
    }
    
    // Compute stiffness matrix K = sum_{i=1}^4 (C_i^T * G * C_i) * (|det(H)|/6)
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Scale factor: volume * 6? Actually from NASTRAN: HDETER/6.0D0 for tetra
    double scale = std::abs(hdet) / 6.0;
    
    // For each node i (0-3), add contribution to global stiffness
    for (int i = 0; i < 4; ++i) {
        // Local 3x3 stiffness block: C_i^T * G * C_i
        Eigen::Matrix<double, 3, 3> k_local = C[i].transpose() * G * C[i];
        
        // Map to global 12x12 matrix (3 DOF per node)
        int start_row = i * 3;
        int start_col = i * 3;
        
        K.block(start_row, start_col, 3, 3) += k_local * scale;
        
        // Off-diagonal terms: C_i^T * G * C_j for i != j
        for (int j = 0; j < 4; ++j) {
            if (i == j) continue;
            
            Eigen::Matrix<double, 3, 3> k_off = C[i].transpose() * G * C[j];
            int start_row_j = j * 3;
            int start_col_j = j * 3;
            
            K.block(start_row, start_col_j, 3, 3) += k_off * scale;
            K.block(start_row_j, start_col, 3, 3) += k_off.transpose() * scale;
        }
    }
    
    return K;
}

// Convert matrix to JSON string
std::string matrixToJson(const Eigen::Matrix<double, 12, 12>& K) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(15) << K(i,j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    try {
        // Test case: Tetrahedron with nodes at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
        std::vector<Eigen::Vector3d> nodes = {
            Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(0.0, 1.0, 0.0),
            Eigen::Vector3d(0.0, 0.0, 1.0)
        };
        
        // Material properties
        double E = 200e9;      // 200 GPa
        double nu = 0.3;
        
        // Compute stiffness matrix
        Eigen::Matrix<double, 12, 12> K = computeKTETRA(nodes, E, nu);
        
        // Output as JSON
        std::cout << matrixToJson(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}