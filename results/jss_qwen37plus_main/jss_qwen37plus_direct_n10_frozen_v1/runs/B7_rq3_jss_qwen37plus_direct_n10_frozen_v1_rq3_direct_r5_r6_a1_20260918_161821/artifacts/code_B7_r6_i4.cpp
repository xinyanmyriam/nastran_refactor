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
    
    // Extract the 3x4 submatrix of H_inv (rows 1-3, cols 0-3) which contains dN_i/dx, dN_i/dy, dN_i/dz
    // The correct relationship is: dN_i/dx = H_inv(1,i), dN_i/dy = H_inv(2,i), dN_i/dz = H_inv(3,i)
    // because N = H_inv * [1,x,y,z]^T implies dN/dx = H_inv * [0,1,0,0]^T = column 1 of H_inv
    Eigen::Matrix<double, 4, 3> dNdx;
    for (int i = 0; i < 4; ++i) {
        dNdx(i, 0) = H_inv(1, i); // dN_i/dx
        dNdx(i, 1) = H_inv(2, i); // dN_i/dy
        dNdx(i, 2) = H_inv(3, i); // dN_i/dz
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