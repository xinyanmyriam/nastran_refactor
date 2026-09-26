#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Helper function to compute determinant of 4x4 matrix
double det4x4(const Eigen::Matrix4d& M) {
    return M.determinant();
}

// Helper function to invert 4x4 matrix and get determinant
bool invert4x4(const Eigen::Matrix4d& M, Eigen::Matrix4d& invM, double& det) {
    det = det4x4(M);
    if (std::abs(det) < 1e-15) {
        return false;
    }
    invM = M.inverse();
    return true;
}

// Main KTETRA function that computes the 12x12 stiffness matrix
Eigen::Matrix<double, 12, 12> KTETRA() {
    // Test case: Tetrahedron N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // Material properties
    const double E = 200e9;  // Pa
    const double nu = 0.3;

    // Node coordinates (4 nodes, 3 DOF each)
    Eigen::Matrix<double, 4, 3> coords;
    coords << 0.0, 0.0, 0.0,
              1.0, 0.0, 0.0,
              0.0, 1.0, 0.0,
              0.0, 0.0, 1.0;

    // Build H matrix (4x4): [1 x y z] for each node
    Eigen::Matrix4d H;
    H << 1.0, coords(0,0), coords(0,1), coords(0,2),
         1.0, coords(1,0), coords(1,1), coords(1,2),
         1.0, coords(2,0), coords(2,1), coords(2,2),
         1.0, coords(3,0), coords(3,1), coords(3,2);

    // Invert H and get determinant
    Eigen::Matrix4d H_inv;
    double H_det;
    bool invert_ok = invert4x4(H, H_inv, H_det);
    if (!invert_ok || H_det <= 0.0) {
        throw std::runtime_error("Bad or reverse tetrahedron geometry");
    }

    // Volume = |det(H)| / 6.0
    double volume = std::abs(H_det) / 6.0;

    // Compute material matrix G (6x6 stress-strain matrix for isotropic material)
    Eigen::Matrix<double, 6, 6> G = Eigen::Matrix<double, 6, 6>::Zero();
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-6) {
        throw std::runtime_error("Illegal value of Poisson's ratio");
    }
    double C1 = E * (1.0 - nu) / temp1;
    double C2 = E * nu / temp1;
    double GG = E / (2.0 * (1.0 + nu)); // Shear modulus

    G(0,0) = C1; G(1,1) = C1; G(2,2) = C1;
    G(0,1) = C2; G(0,2) = C2;
    G(1,0) = C2; G(1,2) = C2;
    G(2,0) = C2; G(2,1) = C2;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    // Build C matrices (4 matrices of size 6x3)
    // For each node i (0-indexed), C_i is 6x3
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    
    // The Fortran code builds C as:
    // For node i (1-indexed), C(i) has rows:
    // row1: H(i+4), 0, 0
    // row2: 0, H(i+8), 0  
    // row3: 0, 0, H(i+12)
    // row4: 0, H(i+12), H(i+8)
    // row5: H(i+12), 0, H(i+4)
    // row6: H(i+8), H(i+4), 0
    // But looking at the Fortran loop (I=1,4) and indexing J=18*I-18, 
    // and the assignments: C(J+1)=H(I+4), C(J+5)=H(I+8), C(J+9)=H(I+12), etc.
    // Actually, the Fortran uses column-major storage and the pattern is:
    // Each C_i is stored as 6 rows x 3 columns, with values from H_inv rows.
    // From the Fortran: H_inv contains the inverse of H, and the C matrices
    // are built from rows of H_inv (which are the shape function derivatives).
    
    // Correct approach: The gradient of shape functions is B = H_inv^T * [dN/dxi]
    // For linear tetrahedron, the shape function derivatives in natural coords are constant.
    // The matrix H_inv gives us the transformation from global to natural coordinates.
    // So the rows of H_inv (excluding first column) are the gradients.
    
    // Extract the 3x3 submatrix of H_inv (rows 0-3, cols 1-3) which contains dN_i/dx, dN_i/dy, dN_i/dz
    // Actually, for linear tetrahedron, the strain-displacement matrix B is built as:
    // B = [ dN1/dx  0      0      dN2/dx  0      0      dN3/dx  0      0      dN4/dx  0      0    ;
    //       0      dN1/dy 0      0      dN2/dy 0      0      dN3/dy 0      0      dN4/dy 0    ;
    //       0      0      dN1/dz 0      0      dN2/dz 0      0      dN3/dz 0      0      dN4/dz;
    //       dN1/dy dN1/dx 0      dN2/dy dN2/dx 0      dN3/dy dN3/dx 0      dN4/dy dN4/dx 0    ;
    //       dN1/dz 0      dN1/dx dN2/dz 0      dN2/dx dN3/dz 0      dN3/dx dN4/dz 0      dN4/dx;
    //       0      dN1/dz dN1/dy 0      dN2/dz dN2/dy 0      dN3/dz dN3/dy 0      dN4/dz dN4/dy]
    
    // The derivatives dN_i/dx, dN_i/dy, dN_i/dz are the rows of H_inv(1:4,1:3) transposed
    // because H = [1 x y z], so H_inv * [1 x y z]^T = [N1 N2 N3 N4]^T
    // Thus, dN_i/dx, dN_i/dy, dN_i/dz are in row i of H_inv, columns 1,2,3
    
    Eigen::Matrix<double, 4, 3> dNdx;
    for (int i = 0; i < 4; ++i) {
        dNdx(i, 0) = H_inv(i, 1); // dN_i/dx
        dNdx(i, 1) = H_inv(i, 2); // dN_i/dy
        dNdx(i, 2) = H_inv(i, 3); // dN_i/dz
    }
    
    // Now build the 6x12 B matrix (strain-displacement)
    Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        int col_base = i * 3;
        // Row 0: dN_i/dx, 0, 0
        B(0, col_base + 0) = dNdx(i, 0);
        // Row 1: 0, dN_i/dy, 0
        B(1, col_base + 1) = dNdx(i, 1);
        // Row 2: 0, 0, dN_i/dz
        B(2, col_base + 2) = dNdx(i, 2);
        // Row 3: dN_i/dy, dN_i/dx, 0
        B(3, col_base + 0) = dNdx(i, 1);
        B(3, col_base + 1) = dNdx(i, 0);
        // Row 4: dN_i/dz, 0, dN_i/dx
        B(4, col_base + 0) = dNdx(i, 2);
        B(4, col_base + 2) = dNdx(i, 0);
        // Row 5: 0, dN_i/dz, dN_i/dy
        B(5, col_base + 1) = dNdx(i, 2);
        B(5, col_base + 2) = dNdx(i, 1);
    }
    
    // Stiffness matrix: K = B^T * G * B * volume
    Eigen::Matrix<double, 12, 12> K = volume * B.transpose() * G * B;
    
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
            oss << std::scientific << std::setprecision(6) << K(i, j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    try {
        // Compute the stiffness matrix
        Eigen::Matrix<double, 12, 12> K = KTETRA();
        
        // Output as JSON
        std::cout << matrixToJson(K) << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}