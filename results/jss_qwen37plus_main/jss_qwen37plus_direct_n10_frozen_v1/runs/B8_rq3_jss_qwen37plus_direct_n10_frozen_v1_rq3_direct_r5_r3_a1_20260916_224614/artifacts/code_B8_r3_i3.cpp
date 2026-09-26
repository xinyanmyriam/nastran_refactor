#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

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

// Helper function to compute the stiffness contribution of a tetrahedron
void computeTetStiffness(const std::vector<Eigen::Vector3d>& nodes,
                        double E, double nu,
                        Eigen::MatrixXd& K_tet) {
    // Build the H matrix: [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix4d H;
    for (int i = 0; i < 4; ++i) {
        H(i, 0) = 1.0;
        H(i, 1) = nodes[i](0);
        H(i, 2) = nodes[i](1);
        H(i, 3) = nodes[i](2);
    }

    // Invert H matrix
    Eigen::Matrix4d H_inv;
    double det;
    if (!invert4x4(H, H_inv, det)) {
        // Bad geometry - set zero matrix
        K_tet.setZero();
        return;
    }

    // Compute material matrix G (6x6)
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-15) {
        K_tet.setZero();
        return;
    }
    
    double g11 = E * (1.0 - nu) / temp1;
    double g12 = E * nu / temp1;
    double gg = E / (2.0 * (1.0 + nu)); // Shear modulus
    
    Eigen::Matrix<double, 6, 6> G;
    G.setZero();
    G(0,0) = g11; G(0,1) = g12; G(0,2) = g12;
    G(1,0) = g12; G(1,1) = g11; G(1,2) = g12;
    G(2,0) = g12; G(2,1) = g12; G(2,2) = g11;
    G(3,3) = gg; G(4,4) = gg; G(5,5) = gg;

    // Build C matrices (4 matrices of size 6x3)
    // For each node i, C_i is 6x3 with:
    // row 0: [dN_i/dx, 0, 0]
    // row 1: [0, dN_i/dy, 0]
    // row 2: [0, 0, dN_i/dz]
    // row 3: [0, dN_i/dz, dN_i/dy]
    // row 4: [dN_i/dz, 0, dN_i/dx]
    // row 5: [dN_i/dy, dN_i/dx, 0]
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    for (int i = 0; i < 4; ++i) {
        double dx = H_inv(i,1);
        double dy = H_inv(i,2);
        double dz = H_inv(i,3);
        
        C[i].setZero();
        C[i](0,0) = dx;
        C[i](1,1) = dy;
        C[i](2,2) = dz;
        C[i](3,1) = dz; C[i](3,2) = dy;
        C[i](4,0) = dz; C[i](4,2) = dx;
        C[i](5,0) = dy; C[i](5,1) = dx;
    }

    // Volume factor: |det(H)|/6.0
    double volume_factor = std::abs(det) / 6.0;

    // Initialize stiffness matrix for this tetrahedron (12x12 since 4 nodes * 3 DOF)
    K_tet.setZero();
    
    // Compute k_ij = C_i^T * G * C_j * volume_factor
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // C_i^T * G * C_j
            Eigen::Matrix<double, 3, 3> k_sub = C[i].transpose() * G * C[j];
            
            // Place in global 12x12 matrix
            for (int di = 0; di < 3; ++di) {
                for (int dj = 0; dj < 3; ++dj) {
                    int row = i*3 + di;
                    int col = j*3 + dj;
                    K_tet(row, col) += k_sub(di, dj) * volume_factor;
                }
            }
        }
    }
}

// Main KSOLID subroutine equivalent
void ksolid(std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords) {
    // Correct decomposition of wedge (6 nodes) into 3 tets:
    // Tet 1: nodes 0,1,2,3 (N1,N2,N3,N4)
    // Tet 2: nodes 1,2,4,5 (N2,N3,N5,N6)
    // Tet 3: nodes 0,2,3,5 (N1,N3,N4,N6)
    std::vector<std::vector<int>> wedge_tets = {
        {0, 1, 2, 3}, // N1,N2,N3,N4
        {1, 2, 4, 5}, // N2,N3,N5,N6
        {0, 2, 3, 5}  // N1,N3,N4,N6
    };
    
    // Initialize 18x18 stiffness matrix
    stiffness_matrix.setZero();
    
    // Process each tetrahedron
    for (const auto& tet_nodes : wedge_tets) {
        std::vector<Eigen::Vector3d> tet_coords;
        for (int idx : tet_nodes) {
            tet_coords.push_back(Eigen::Vector3d(
                node_coords[idx][0],
                node_coords[idx][1],
                node_coords[idx][2]
            ));
        }
        
        Eigen::MatrixXd K_tet(12, 12); // 4 nodes * 3 DOF = 12
        computeTetStiffness(tet_coords, E, nu, K_tet);
        
        // Map local 12x12 stiffness to global 18x18
        // Each tet has 4 nodes, each with 3 DOF
        for (int i = 0; i < 4; ++i) {
            int global_node_i = tet_nodes[i];
            for (int j = 0; j < 4; ++j) {
                int global_node_j = tet_nodes[j];
                for (int di = 0; di < 3; ++di) {
                    for (int dj = 0; dj < 3; ++dj) {
                        int local_row = i*3 + di;
                        int local_col = j*3 + dj;
                        int global_row = global_node_i*3 + di;
                        int global_col = global_node_j*3 + dj;
                        stiffness_matrix(global_row, global_col) += K_tet(local_row, local_col);
                    }
                }
            }
        }
    }
}

// Main function
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
    
    double E = 200e9;      // 200 GPa
    double nu = 0.3;
    
    // Initialize ECPT array (not fully used but for compatibility)
    std::vector<double> ecpt(100, 0.0);
    
    // Initialize stiffness matrix (18x18)
    Eigen::MatrixXd stiffness_matrix(18, 18);
    
    // Call KSOLID equivalent
    ksolid(ecpt, stiffness_matrix, E, nu, node_coords);
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 18; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 18; ++j) {
            if (j > 0) std::cout << ",";
            // Format in scientific notation with sufficient precision
            std::cout << std::scientific << std::setprecision(15) << stiffness_matrix(i, j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;
    
    return 0;
}