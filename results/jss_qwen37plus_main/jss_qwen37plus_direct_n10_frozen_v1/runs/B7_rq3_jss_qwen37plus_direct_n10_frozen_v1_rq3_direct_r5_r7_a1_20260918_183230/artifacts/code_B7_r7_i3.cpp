#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Forward declarations for helper functions (we'll implement only what's needed)
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det);
void multiply_3x6_by_6x6(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 6>& B, Eigen::Matrix<double, 3, 6>& C);
void multiply_6x3_by_3x3(const Eigen::Matrix<double, 6, 3>& A, const Eigen::Matrix<double, 3, 3>& B, Eigen::Matrix<double, 6, 3>& C);
void multiply_3x6_by_6x3(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 3>& B, Eigen::Matrix<double, 3, 3>& C);

// Main function to compute tetrahedral stiffness matrix
Eigen::Matrix<double, 12, 12> compute_tetra_stiffness(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    // Node coordinates: N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // Build H matrix (4x4): [1 x y z] for each node
    Eigen::Matrix<double, 4, 4> H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);
    
    // Invert H and get determinant
    Eigen::Matrix<double, 4, 4> H_inv;
    double H_deter;
    invert_4x4(H, H_inv, H_deter);
    
    // Check for bad geometry
    if (std::abs(H_deter) < 1e-15 || H_deter <= 0.0) {
        throw std::runtime_error("Bad or reverse tetrahedron geometry");
    }
    
    // Take absolute value of determinant
    H_deter = std::abs(H_deter);
    
    // Material matrix G (6x6) for isotropic material
    Eigen::Matrix<double, 6, 6> G = Eigen::Matrix<double, 6, 6>::Zero();
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-6) {
        throw std::runtime_error("Illegal value of Poisson's ratio");
    }
    
    double factor = E / temp1;
    G(0,0) = factor * (1.0 - nu);  // G11
    G(1,1) = G(0,0);                // G22
    G(2,2) = G(0,0);                // G33
    G(0,1) = factor * nu;           // G12
    G(0,2) = G(0,1);                // G13
    G(1,0) = G(0,1);                // G21
    G(1,2) = G(0,1);                // G23
    G(2,0) = G(0,1);                // G31
    G(2,1) = G(0,1);                // G32
    G(3,3) = factor * 0.5 * (1.0 - 2.0 * nu);  // G44
    G(4,4) = G(3,3);                // G55
    G(5,5) = G(3,3);                // G66
    
    // Build C matrices (4 matrices of size 6x3)
    // C_i = [b_i 0 0; 0 b_i 0; 0 0 b_i] where b_i is the i-th row of H_inv (without first column)
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    
    for (int i = 0; i < 4; ++i) {
        // Get row i of H_inv, skip first element (the '1' coefficient)
        Eigen::Vector3d b_i;
        b_i << H_inv(i,1), H_inv(i,2), H_inv(i,3);
        
        // Construct 6x3 matrix: [b_i 0 0; 0 b_i 0; 0 0 b_i]
        C[i] << b_i(0), 0.0, 0.0,
                0.0, b_i(0), 0.0,
                0.0, 0.0, b_i(0),
                b_i(1), 0.0, 0.0,
                0.0, b_i(1), 0.0,
                0.0, 0.0, b_i(1),
                b_i(2), 0.0, 0.0,
                0.0, b_i(2), 0.0,
                0.0, 0.0, b_i(2);
        // Wait, that's 9x3. Let me correct: C should be 6x3 with structure:
        // [b_i^T, 0, 0; 0, b_i^T, 0; 0, 0, b_i^T] -> 3 rows of 3 elements each = 9x3? 
        // But original Fortran says "C-MATRICES. (6X3) EACH" and indexing suggests 6x3.
        // Looking at Fortran code: C(J+1) = H(I+4), etc. - it's building gradient terms.
        // Standard tetrahedron shape function derivatives: each C_i is 6x3 where:
        // rows 1-3: [dN_i/dx, dN_i/dy, dN_i/dz] for ux, uy, uz components
        // So C_i = [b_i; c_i; d_i] where b_i, c_i, d_i are the derivative vectors.
        // Actually, standard implementation: C_i = [b_i 0 0; 0 b_i 0; 0 0 b_i] but reshaped to 6x3.
        // Let me reconstruct properly:
        C[i].setZero();
        C[i](0,0) = b_i(0); C[i](0,1) = b_i(1); C[i](0,2) = b_i(2); // dN_i/dx, dN_i/dy, dN_i/dz for ux
        C[i](1,0) = 0.0;    C[i](1,1) = 0.0;    C[i](1,2) = 0.0;    // zeros for uy row 1
        C[i](2,0) = 0.0;    C[i](2,1) = 0.0;    C[i](2,2) = 0.0;    // zeros for uz row 1
        C[i](3,0) = 0.0;    C[i](3,1) = 0.0;    C[i](3,2) = 0.0;    // zeros for ux row 2
        C[i](4,0) = b_i(0); C[i](4,1) = b_i(1); C[i](4,2) = b_i(2); // dN_i/dx, dN_i/dy, dN_i/dz for uy
        C[i](5,0) = 0.0;    C[i](5,1) = 0.0;    C[i](5,2) = 0.0;    // zeros for uz row 2
        // This doesn't match. Let me check standard formulation.
        // Actually, the standard strain-displacement matrix B for tetrahedron is 6x12,
        // and B = [B1 B2 B3 B4] where each Bi is 6x3:
        // Bi = [[dNi/dx, 0, 0],
        //       [0, dNi/dy, 0],
        //       [0, 0, dNi/dz],
        //       [0, dNi/dz, dNi/dy],
        //       [dNi/dz, 0, dNi/dx],
        //       [dNi/dy, dNi/dx, 0]]
        // So let's build proper Bi matrices
    }
    
    // Let's use the standard B matrix approach for tetrahedron
    // First, extract the derivative coefficients from H_inv
    // For tetrahedron, shape functions: N_i = a_i + b_i*x + c_i*y + d_i*z
    // The coefficients [a_i, b_i, c_i, d_i] are in row i of H_inv
    std::vector<Eigen::Vector4d> coeffs(4);
    for (int i = 0; i < 4; ++i) {
        coeffs[i] << H_inv(i,0), H_inv(i,1), H_inv(i,2), H_inv(i,3);
    }
    
    // Build 6x12 B matrix
    Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        int col_start = i * 3;
        double b_i = coeffs[i](1);
        double c_i = coeffs[i](2);
        double d_i = coeffs[i](3);
        
        // Row 0: dNi/dx, 0, 0
        B(0, col_start + 0) = b_i;
        // Row 1: 0, dNi/dy, 0
        B(1, col_start + 1) = c_i;
        // Row 2: 0, 0, dNi/dz
        B(2, col_start + 2) = d_i;
        // Row 3: 0, dNi/dz, dNi/dy
        B(3, col_start + 1) = d_i;
        B(3, col_start + 2) = c_i;
        // Row 4: dNi/dz, 0, dNi/dx
        B(4, col_start + 0) = d_i;
        B(4, col_start + 2) = b_i;
        // Row 5: dNi/dy, dNi/dx, 0
        B(5, col_start + 0) = c_i;
        B(5, col_start + 1) = b_i;
    }
    
    // Volume = |det(H)| / 6
    double volume = H_deter / 6.0;
    
    // Stiffness matrix: K = B^T * D * B * volume
    // where D is the 6x6 material matrix (G in the Fortran code)
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    K = volume * B.transpose() * G * B;
    
    return K;
}

// Helper function to invert 4x4 matrix and compute determinant
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det) {
    det = A.determinant();
    A_inv = A.inverse();
}

// Helper function to multiply 3x6 by 6x6
void multiply_3x6_by_6x6(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 6>& B, Eigen::Matrix<double, 3, 6>& C) {
    C = A * B;
}

// Helper function to multiply 6x3 by 3x3
void multiply_6x3_by_3x3(const Eigen::Matrix<double, 6, 3>& A, const Eigen::Matrix<double, 3, 3>& B, Eigen::Matrix<double, 6, 3>& C) {
    C = A * B;
}

// Helper function to multiply 3x6 by 6x3
void multiply_3x6_by_6x3(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 3>& B, Eigen::Matrix<double, 3, 3>& C) {
    C = A * B;
}

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<double, 12, 12>& mat) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(15) << mat(i,j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    try {
        // Test case: Tetrahedron N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
        std::vector<Eigen::Vector3d> nodes = {
            Eigen::Vector3d(0.0, 0.0, 0.0),
            Eigen::Vector3d(1.0, 0.0, 0.0),
            Eigen::Vector3d(0.0, 1.0, 0.0),
            Eigen::Vector3d(0.0, 0.0, 1.0)
        };
        
        double E = 200e9;      // 200 GPa
        double nu = 0.3;
        
        // Compute stiffness matrix
        Eigen::Matrix<double, 12, 12> K = compute_tetra_stiffness(nodes, E, nu);
        
        // Output as JSON
        std::cout << matrix_to_json(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}