#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Forward declarations for helper functions (we'll implement only what's needed)
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det);
void matrix_multiply(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 6>& B, Eigen::Matrix<double, 3, 6>& C);
void matrix_multiply(const Eigen::Matrix<double, 6, 3>& A, const Eigen::Matrix<double, 3, 3>& B, Eigen::Matrix<double, 6, 3>& C);
void matrix_multiply(const Eigen::Matrix<double, 3, 6>& A, const Eigen::Matrix<double, 6, 3>& B, Eigen::Matrix<double, 3, 3>& C);

// Main function to compute tetrahedral stiffness matrix
Eigen::Matrix<double, 12, 12> KTETRA() {
    // Test case: Tetrahedron nodes
    // N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    std::vector<Eigen::Vector3d> nodes = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 1.0)
    };

    // Material properties
    const double E = 200e9;      // Young's modulus
    const double nu = 0.3;       // Poisson's ratio

    // Build H matrix (4x4) as in Fortran: [1 x y z] for each node
    Eigen::Matrix<double, 4, 4> H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);

    // Invert H and get determinant
    Eigen::Matrix<double, 4, 4> H_inv;
    double H_deter;
    invert_4x4(H, H_inv, H_deter);

    // Check for bad geometry (determinant <= 0)
    if (H_deter <= 0.0) {
        throw std::runtime_error("Bad or reverse tetrahedron geometry: determinant <= 0");
    }

    // Absolute value of determinant for volume calculation
    H_deter = std::abs(H_deter);

    // Compute material matrix G (6x6 stress-strain matrix for isotropic material)
    // G = [C11 C12 C12 0   0   0  ]
    //     [C12 C11 C12 0   0   0  ]
    //     [C12 C12 C11 0   0   0  ]
    //     [0   0   0   C44 0   0  ]
    //     [0   0   0   0   C44 0  ]
    //     [0   0   0   0   0   C44]
    // where C11 = E*(1-nu)/((1+nu)*(1-2*nu)), C12 = E*nu/((1+nu)*(1-2*nu)), C44 = E/(2*(1+nu))
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Illegal value of Poisson's ratio: denominator near zero");
    }

    double C11 = E * (1.0 - nu) / temp1;
    double C12 = E * nu / temp1;
    double C44 = E / (2.0 * (1.0 + nu));

    Eigen::Matrix<double, 6, 6> G;
    G.setZero();
    G(0,0) = C11; G(0,1) = C12; G(0,2) = C12; G(0,3) = 0.0;  G(0,4) = 0.0;  G(0,5) = 0.0;
    G(1,0) = C12; G(1,1) = C11; G(1,2) = C12; G(1,3) = 0.0;  G(1,4) = 0.0;  G(1,5) = 0.0;
    G(2,0) = C12; G(2,1) = C12; G(2,2) = C11; G(2,3) = 0.0;  G(2,4) = 0.0;  G(2,5) = 0.0;
    G(3,0) = 0.0;  G(3,1) = 0.0;  G(3,2) = 0.0;  G(3,3) = C44; G(3,4) = 0.0;  G(3,5) = 0.0;
    G(4,0) = 0.0;  G(4,1) = 0.0;  G(4,2) = 0.0;  G(4,3) = 0.0;  G(4,4) = C44; G(4,5) = 0.0;
    G(5,0) = 0.0;  G(5,1) = 0.0;  G(5,2) = 0.0;  G(5,3) = 0.0;  G(5,4) = 0.0;  G(5,5) = C44;

    // Build C matrices (4 matrices of size 6x3)
    // From Fortran: for i=1..4, C_i is 6x3 with specific pattern using H columns
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    
    for (int i = 0; i < 4; ++i) {
        // H(i+1, :) is row i+1 of H (1-indexed in Fortran, so row i in 0-indexed)
        // But Fortran code uses H(I+4), H(I+8), H(I+12) which are column elements
        // Looking at Fortran: C(J+1) = H(I+4), etc., where J = 18*I-18, so for I=1: J=0
        // H is stored column-wise in Fortran, but we have it as row-major in Eigen
        // So H(I+4) means element at position (I+3, 0) in 0-indexed? Let's reinterpret.
        
        // Actually, from Fortran code:
        // H(1) = 1.0D0          -> H(0,0)
        // H(2) = ECPT(8)       -> H(0,1) = x1
        // H(3) = ECPT(9)       -> H(0,2) = y1  
        // H(4) = ECPT(10)      -> H(0,3) = z1
        // H(5) = 1.0D0          -> H(1,0)
        // H(6) = ECPT(12)      -> H(1,1) = x2
        // H(7) = ECPT(13)      -> H(1,2) = y2
        // H(8) = ECPT(14)      -> H(1,3) = z2
        // So H is stored row-wise: row0: [1,x1,y1,z1], row1: [1,x2,y2,z2], etc.
        
        // Therefore H(I+4) for I=1 is H(5) = element at index 4 = row1,col0 = 1.0
        // H(I+8) for I=1 is H(9) = element at index 8 = row2,col0 = 1.0
        // H(I+12) for I=1 is H(13) = element at index 12 = row3,col0 = 1.0
        
        // So the pattern is taking column 0 of H for all rows, then column 1, etc.
        // But looking more carefully at the indexing:
        // For I=1: J=0, so C(1)=H(5)=H(1,0), C(5)=H(9)=H(2,0), C(9)=H(13)=H(3,0)
        // C(11)=H(13)=H(3,0), C(12)=H(9)=H(2,0), C(13)=H(13)=H(3,0), C(15)=H(5)=H(1,0), C(16)=H(9)=H(2,0), C(17)=H(5)=H(1,0)
        
        // This suggests the C matrix for node i is built from the i-th column of H (excluding first column)
        // Actually, let's reconstruct from the pattern:
        // The C matrix for node i should be related to the shape function derivatives.
        // Standard tetrahedral element: shape functions N_i = a_i + b_i*x + c_i*y + d_i*z
        // where [a_i b_i c_i d_i] = H^{-1} * e_i, and e_i is unit vector.
        // The strain-displacement matrix B has blocks [dN_i/dx, 0, 0; 0, dN_i/dy, 0; 0, 0, dN_i/dz; dN_i/dy, dN_i/dx, 0; dN_i/dz, 0, dN_i/dx; 0, dN_i/dz, dN_i/dy]
        
        // From standard FEM theory for tetrahedron:
        // Let H = [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
        // Then H_inv = [a1 a2 a3 a4; b1 b2 b3 b4; c1 c2 c3 c4; d1 d2 d3 d4]
        // So dN_i/dx = b_i, dN_i/dy = c_i, dN_i/dz = d_i
        
        // So for node i, the gradient is [b_i, c_i, d_i]^T
        // Then the B matrix block for node i is:
        // [b_i  0    0  ]
        // [0    c_i  0  ]
        // [0    0    d_i]
        // [c_i  b_i  0  ]
        // [d_i  0    b_i]
        // [0    d_i  c_i]
        
        // So C_i should be this 6x3 matrix.
        double b_i = H_inv(1, i);  // row 1, col i (since H_inv * H = I, so row 1 gives coefficients for x)
        double c_i = H_inv(2, i);  // row 2 gives coefficients for y
        double d_i = H_inv(3, i);  // row 3 gives coefficients for z
        
        C[i] << b_i, 0.0, 0.0,
                0.0, c_i, 0.0,
                0.0, 0.0, d_i,
                c_i, b_i, 0.0,
                d_i, 0.0, b_i,
                0.0, d_i, c_i;
    }

    // Scale factor: divide determinant by 6.0 (from Fortran: HDETER = HDETER/6.0D0)
    double scale_factor = H_deter / 6.0;

    // Initialize global stiffness matrix (12x12)
    Eigen::Matrix<double, 12, 12> K;
    K.setZero();

    // For each node i (1-4), compute contribution K_i = C_i^T * G * C_i
    for (int i = 0; i < 4; ++i) {
        // Compute G * C_i (6x3)
        Eigen::Matrix<double, 6, 3> GC;
        GC = G * C[i];
        
        // Compute C_i^T * (G * C_i) (3x3)
        Eigen::Matrix<double, 3, 3> K_local;
        K_local = C[i].transpose() * GC;
        
        // Scale by determinant factor
        K_local *= scale_factor;
        
        // Place 3x3 local stiffness into global 12x12 matrix at appropriate location
        // Node i has DOFs at positions [3*i, 3*i+1, 3*i+2] for ux, uy, uz
        int start_row = i * 3;
        int start_col = i * 3;
        
        K.block(start_row, start_col, 3, 3) += K_local;
        
        // Also add coupling terms with other nodes (since K_global = sum C_i^T * G * C_j)
        // For tetrahedron, the full stiffness is sum_{i,j} C_i^T * G * C_j
        for (int j = 0; j < 4; ++j) {
            if (i == j) continue;
            
            // Compute C_i^T * G * C_j (3x3)
            Eigen::Matrix<double, 6, 3> GCj = G * C[j];
            Eigen::Matrix<double, 3, 3> K_coupling = C[i].transpose() * GCj;
            K_coupling *= scale_factor;
            
            int start_row_j = j * 3;
            int start_col_j = j * 3;
            
            K.block(start_row, start_col_j, 3, 3) += K_coupling;
        }
    }

    return K;
}

// Helper function to invert 4x4 matrix and compute determinant
void invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& A_inv, double& det) {
    // Use Eigen's LU decomposition for robust inversion
    Eigen::FullPivLU<Eigen::Matrix<double, 4, 4>> lu(A);
    det = lu.determinant();
    A_inv = lu.inverse();
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

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<double, 12, 12>& M) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        oss << "[";
        for (int j = 0; j < 12; ++j) {
            oss << std::scientific << std::setprecision(15) << M(i, j);
            if (j < 11) oss << ",";
        }
        oss << "]";
        if (i < 11) oss << ",";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    try {
        // Compute the stiffness matrix
        Eigen::Matrix<double, 12, 12> K = KTETRA();
        
        // Output as JSON
        std::cout << matrix_to_json(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}