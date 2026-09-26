#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Forward declarations for helper functions (we'll implement only what's needed)
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det);
void gmmatd(const Eigen::MatrixXd& A, int rowsA, int colsA, int transA,
            const Eigen::MatrixXd& B, int rowsB, int colsB, int transB,
            Eigen::MatrixXd& C);

// Main function to compute tetrahedral stiffness matrix
Eigen::Matrix<double, 12, 12> ktetra() {
    // Test case: Tetrahedron N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // Material: E = 200e9, nu = 0.3
    const double E = 200e9;
    const double nu = 0.3;

    // Node coordinates (4 nodes, 3 DOF each)
    std::vector<Eigen::Vector3d> nodes = {
        Eigen::Vector3d(0.0, 0.0, 0.0), // N1
        Eigen::Vector3d(1.0, 0.0, 0.0), // N2
        Eigen::Vector3d(0.0, 1.0, 0.0), // N3
        Eigen::Vector3d(0.0, 0.0, 1.0)  // N4
    };

    // Build H matrix (4x4): [1 x y z] for each node
    Eigen::Matrix<double, 4, 4> H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);

    // Invert H and get determinant
    Eigen::Matrix<double, 4, 4> H_inv;
    double h_deter;
    invert_4x4(H, H_inv, h_deter);

    // Check for bad geometry
    if (std::abs(h_deter) < 1e-15 || h_deter <= 0.0) {
        throw std::runtime_error("Bad or reverse geometry: determinant <= 0");
    }

    // Volume factor: |det(H)| / 6.0
    double volume_factor = std::abs(h_deter) / 6.0;

    // Build material matrix G (6x6) for isotropic linear elasticity
    // G = D matrix in Voigt notation
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Illegal value of nu");
    }

    double c11 = E * (1.0 - nu) / temp1;
    double c12 = E * nu / temp1;
    double c44 = E / (2.0 * (1.0 + nu));

    Eigen::Matrix<double, 6, 6> G = Eigen::Matrix<double, 6, 6>::Zero();
    G(0,0) = c11; G(0,1) = c12; G(0,2) = c12;
    G(1,0) = c12; G(1,1) = c11; G(1,2) = c12;
    G(2,0) = c12; G(2,1) = c12; G(2,2) = c11;
    G(3,3) = c44;
    G(4,4) = c44;
    G(5,5) = c44;

    // Build C matrices (4 matrices of size 6x3)
    // For tetrahedral element, the strain-displacement matrix B is built from shape function derivatives
    // The C matrix here corresponds to the derivative matrix for each node
    // From Fortran: C(J+1) = H(I+4), etc. - this is extracting rows of H_inv
    // Actually, the C matrices are the shape function derivatives: dN_i/dx, dN_i/dy, dN_i/dz
    // For linear tetrahedron: N_i = a_i + b_i*x + c_i*y + d_i*z, so derivatives are constants
    // The derivatives are given by the last 3 columns of H_inv (since H = [1 x y z])

    std::vector<Eigen::Matrix<double, 6, 3>> C_matrices(4);
    
    // Each C matrix is 6x3: first 3 rows for epsilon_xx, epsilon_yy, epsilon_zz,
    // last 3 rows for gamma_yz, gamma_xz, gamma_xy
    for (int i = 0; i < 4; ++i) {
        // Get the i-th row of H_inv (shape function coefficients)
        // The derivatives dN_i/dx, dN_i/dy, dN_i/dz are H_inv(i,1), H_inv(i,2), H_inv(i,3)
        double dNdx = H_inv(i, 1);
        double dNdy = H_inv(i, 2);
        double dNdz = H_inv(i, 3);

        // Build C_i matrix (6x3)
        // Row 0: [dNdx, 0, 0] for epsilon_xx
        // Row 1: [0, dNdy, 0] for epsilon_yy  
        // Row 2: [0, 0, dNdz] for epsilon_zz
        // Row 3: [0, dNdz, dNdy] for gamma_yz
        // Row 4: [dNdz, 0, dNdx] for gamma_xz
        // Row 5: [dNdy, dNdx, 0] for gamma_xy
        C_matrices[i] << dNdx, 0.0, 0.0,
                         0.0, dNdy, 0.0,
                         0.0, 0.0, dNdz,
                         0.0, dNdz, dNdy,
                         dNdz, 0.0, dNdx,
                         dNdy, dNdx, 0.0;
    }

    // Initialize global stiffness matrix (12x12)
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();

    // Compute element stiffness: K = sum_i (C_i^T * G * C_i) * volume_factor
    // But note: in standard FEM, K = integral(B^T * D * B) dV
    // For linear tetrahedron, B is constant, so K = B^T * D * B * V
    // Where B is 6x12 matrix built from C_i matrices

    // Build full B matrix (6x12)
    Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        // Node i contributes to columns [3*i, 3*i+2]
        B.block<6,3>(0, 3*i) = C_matrices[i];
    }

    // Compute K = B^T * G * B * volume_factor
    Eigen::Matrix<double, 12, 6> Bt = B.transpose();
    Eigen::Matrix<double, 12, 6> BtG = Bt * G;
    Eigen::Matrix<double, 12, 12> K_local = BtG * B;
    K = K_local * volume_factor;

    return K;
}

// Helper function to invert 4x4 matrix and compute determinant
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det) {
    // Use Eigen's LU decomposition for robust inversion
    Eigen::FullPivLU<Eigen::Matrix<double, 4, 4>> lu(A);
    det = lu.determinant();
    A_inv = lu.inverse();
}

// Simplified version of GMMATD for our specific use case
// We only need: C = A * B (no transpose options needed for this problem)
void gmmatd(const Eigen::MatrixXd& A, int rowsA, int colsA, int transA,
            const Eigen::MatrixXd& B, int rowsB, int colsB, int transB,
            Eigen::MatrixXd& C) {
    Eigen::MatrixXd A_use = A;
    Eigen::MatrixXd B_use = B;
    
    if (transA) A_use = A.transpose();
    if (transB) B_use = B.transpose();
    
    C = A_use * B_use;
}

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<double, 12, 12>& K) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(12) << K(i, j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    try {
        // Compute the stiffness matrix
        Eigen::Matrix<double, 12, 12> K = ktetra();
        
        // Output as JSON
        std::cout << matrix_to_json(K) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}