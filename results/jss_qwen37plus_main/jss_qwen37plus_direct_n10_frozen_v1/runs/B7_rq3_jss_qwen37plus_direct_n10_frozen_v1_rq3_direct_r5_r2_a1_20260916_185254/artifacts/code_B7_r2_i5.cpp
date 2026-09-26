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
void matrix_multiply(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 6>& B, Eigen::Matrix<double, 3, 6>& C);
void matrix_multiply(const Eigen::Matrix<double, 6, 3>& A, const Eigen::Matrix<double, 3, 3>& B, Eigen::Matrix<double, 6, 3>& C);
void matrix_multiply(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 3>& B, Eigen::Matrix<double, 3, 3>& C);

// Main function to compute tetrahedral stiffness matrix
Eigen::Matrix<double, 12, 12> KTETRA(
    const std::vector<std::vector<double>>& nodes, // 4 nodes, each with [x,y,z]
    double E,                                       // Young's modulus
    double nu)                                      // Poisson's ratio
{
    // Node coordinates: N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // We'll use the exact test case values
    
    // Build H matrix (4x4): [1 x y z] for each node
    Eigen::Matrix<double, 4, 4> H;
    H << 1.0, nodes[0][0], nodes[0][1], nodes[0][2],
         1.0, nodes[1][0], nodes[1][1], nodes[1][2],
         1.0, nodes[2][0], nodes[2][1], nodes[2][2],
         1.0, nodes[3][0], nodes[3][1], nodes[3][2];
    
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
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Illegal value of Poisson's ratio");
    }
    
    double factor = E / temp1;
    G(0,0) = factor * (1.0 - nu);  // G11
    G(1,1) = G(0,0);               // G22
    G(2,2) = G(0,0);               // G33
    G(0,1) = factor * nu;          // G12
    G(0,2) = G(0,1);               // G13
    G(1,0) = G(0,1);               // G21
    G(1,2) = G(0,1);               // G23
    G(2,0) = G(0,1);               // G31
    G(2,1) = G(0,1);               // G32
    G(3,3) = factor * (1.0 - 2.0 * nu) / 2.0;  // G44 = G55 = G66
    G(4,4) = G(3,3);
    G(5,5) = G(3,3);
    
    // Build C matrices (4 matrices of size 6x3 each)
    // C_i = [dN_i/dx, dN_i/dy, dN_i/dz] for i=1..4, each is 6x3
    // For tetrahedron, shape functions: N_i = a_i + b_i*x + c_i*y + d_i*z
    // From H matrix inversion, the coefficients are in rows of H_inv
    // dN_i/dx = b_i, dN_i/dy = c_i, dN_i/dz = d_i
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    
    // Each C_i is:
    // [b_i  0    0  ]
    // [0    b_i  0  ]
    // [0    0    b_i]
    // [c_i  0    0  ]
    // [0    c_i  0  ]
    // [0    0    c_i]
    // [d_i  0    0  ]
    // [0    d_i  0  ]
    // [0    0    d_i]
    // But actually for standard tetrahedron, we need the gradient of shape functions
    // The shape function gradients are constant and given by the last 3 columns of H_inv
    // Since H = [1 x y z], then H_inv * [1 x y z]^T = [N1 N2 N3 N4]^T
    // So the gradients are the last 3 rows of H_inv (excluding first column which is for constants)
    
    // Extract gradients: for each node i, grad(N_i) = [H_inv(i,1), H_inv(i,2), H_inv(i,3)]
    std::vector<Eigen::Vector3d> grads(4);
    for (int i = 0; i < 4; ++i) {
        grads[i] << H_inv(i, 1), H_inv(i, 2), H_inv(i, 3);
    }
    
    // Build C matrices: each C_i is 6x3 with block structure
    // [grad(N_i)^T, 0, 0]
    // [0, grad(N_i)^T, 0]
    // [0, 0, grad(N_i)^T]
    for (int i = 0; i < 4; ++i) {
        C[i] = Eigen::Matrix<double, 6, 3>::Zero();
        // First DOF (ux): grad(N_i) in first row
        C[i].row(0) << grads[i](0), grads[i](1), grads[i](2);
        // Second DOF (uy): grad(N_i) in second row  
        C[i].row(1) << grads[i](0), grads[i](1), grads[i](2);
        // Third DOF (uz): grad(N_i) in third row
        C[i].row(2) << grads[i](0), grads[i](1), grads[i](2);
        // But wait - standard B-matrix for tetrahedron is:
        // [dN1/dx 0 0 dN2/dx 0 0 dN3/dx 0 0 dN4/dx 0 0]
        // [0 dN1/dy 0 0 dN2/dy 0 0 dN3/dy 0 0 dN4/dy 0]
        // [0 0 dN1/dz 0 0 dN2/dz 0 0 dN3/dz 0 0 dN4/dz]
        // So actually we need a single 3x12 B matrix, not separate C_i
        
        // Let's reconstruct properly: the strain-displacement matrix B is 6x12
        // For tetrahedron: B = [B1 B2 B3 B4] where each Bi is 6x3
        // Bi = [[dNi/dx, 0, 0],
        //       [0, dNi/dy, 0],
        //       [0, 0, dNi/dz],
        //       [0, dNi/dz, dNi/dy],
        //       [dNi/dz, 0, dNi/dx],
        //       [dNi/dy, dNi/dx, 0]]
    }
    
    // Actually, let's use the standard approach for linear tetrahedron
    // B matrix construction:
    // B = [[b1 0 0 b2 0 0 b3 0 0 b4 0 0],
    //       [0 c1 0 0 c2 0 0 c3 0 0 c4 0],
    //       [0 0 d1 0 0 d2 0 0 d3 0 0 d4],
    //       [0 d1 c1 0 d2 c2 0 d3 c3 0 d4 c4],
    //       [d1 0 c1 d2 0 c2 d3 0 c3 d4 0 c4],
    //       [c1 b1 0 c2 b2 0 c3 b3 0 c4 b4 0]]
    // where for node i: b_i = dN_i/dx, c_i = dN_i/dy, d_i = dN_i/dz
    
    // From H_inv, the shape function coefficients are:
    // N_i = H_inv(i,0) + H_inv(i,1)*x + H_inv(i,2)*y + H_inv(i,3)*z
    // So b_i = H_inv(i,1), c_i = H_inv(i,2), d_i = H_inv(i,3)
    std::vector<double> b(4), c(4), d(4);
    for (int i = 0; i < 4; ++i) {
        b[i] = H_inv(i, 1);
        c[i] = H_inv(i, 2);
        d[i] = H_inv(i, 3);
    }
    
    // Build full 6x12 B matrix
    Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        int col_start = i * 3;
        B(0, col_start + 0) = b[i]; // dN_i/dx for epsilon_xx
        B(1, col_start + 1) = c[i]; // dN_i/dy for epsilon_yy
        B(2, col_start + 2) = d[i]; // dN_i/dz for epsilon_zz
        B(3, col_start + 1) = d[i]; // dN_i/dz for gamma_yz
        B(3, col_start + 2) = c[i]; // dN_i/dy for gamma_yz
        B(4, col_start + 0) = d[i]; // dN_i/dz for gamma_xz
        B(4, col_start + 2) = b[i]; // dN_i/dx for gamma_xz
        B(5, col_start + 0) = c[i]; // dN_i/dy for gamma_xy
        B(5, col_start + 1) = b[i]; // dN_i/dx for gamma_xy
    }
    
    // Volume factor: V = |det(H)| / 6.0
    double volume = std::abs(H_deter) / 6.0;
    
    // Stiffness matrix: K = B^T * G * B * volume
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    K = volume * B.transpose() * G * B;
    
    return K;
}

// Helper function to invert 4x4 matrix and compute determinant
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det) {
    // Use Eigen's built-in inverse and determinant
    det = A.determinant();
    A_inv = A.inverse();
}

// Matrix multiplication helpers
void matrix_multiply(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 6>& B, Eigen::Matrix<double, 3, 6>& C) {
    C = A * B;
}

void matrix_multiply(const Eigen::Matrix<double, 6, 3>& A, const Eigen::Matrix<double, 3, 3>& B, Eigen::Matrix<double, 6, 3>& C) {
    C = A * B;
}

void matrix_multiply(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 3>& B, Eigen::Matrix<double, 3, 3>& C) {
    C = A * B;
}

int main() {
    try {
        // Test case: Tetrahedron N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
        std::vector<std::vector<double>> nodes = {
            {0.0, 0.0, 0.0},
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };
        
        double E = 200e9;      // 200 GPa
        double nu = 0.3;
        
        // Compute stiffness matrix
        Eigen::Matrix<double, 12, 12> K = KTETRA(nodes, E, nu);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        
        for (int i = 0; i < 12; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << "[";
            for (int j = 0; j < 12; ++j) {
                if (j > 0) std::cout << ",";
                // Format in scientific notation with sufficient precision
                std::cout << std::scientific << std::setprecision(12) << K(i,j);
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