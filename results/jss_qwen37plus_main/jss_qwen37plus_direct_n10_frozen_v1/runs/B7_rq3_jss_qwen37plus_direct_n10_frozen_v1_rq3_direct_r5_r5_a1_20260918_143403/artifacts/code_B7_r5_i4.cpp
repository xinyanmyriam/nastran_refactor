#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Helper function to compute determinant of 4x4 matrix
double det4x4(const Eigen::Matrix4d& M) {
    return M.determinant();
}

// Helper function to invert 4x4 matrix and get determinant
bool invert4x4(const Eigen::Matrix4d& M, Eigen::Matrix4d& invM, double& det) {
    det = det4x4(M);
    if (std::abs(det) < 1e-12) {
        return false;
    }
    invM = M.inverse();
    return true;
}

// Compute stiffness matrix for tetrahedral element
Eigen::Matrix<double, 12, 12> computeKTETRA(
    const std::vector<Eigen::Vector3d>& nodes,
    double E,
    double nu) {
    
    // Nodes: N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // Build H matrix (4x4): [1 x y z] for each node
    Eigen::Matrix4d H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);

    double HDETER;
    Eigen::Matrix4d H_inv;
    if (!invert4x4(H, H_inv, HDETER)) {
        throw std::runtime_error("Singular tetrahedron geometry");
    }

    // Take absolute value of determinant
    HDETER = std::abs(HDETER);

    // Material matrix G (6x6) for isotropic material
    Eigen::Matrix<double, 6, 6> G = Eigen::Matrix<double, 6, 6>::Zero();
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    
    double C1 = E * (1.0 - nu) / temp1;
    double C2 = E * nu / temp1;
    double C3 = E / (2.0 * (1.0 + nu)); // GG
    
    G(0,0) = C1; G(1,1) = C1; G(2,2) = C1;
    G(0,1) = C2; G(0,2) = C2;
    G(1,0) = C2; G(1,2) = C2;
    G(2,0) = C2; G(2,1) = C2;
    G(3,3) = C3;
    G(4,4) = C3;
    G(5,5) = C3;

    // Build C matrices (4 matrices of size 6x3)
    // For each node i (0-3), C_i is 6x3
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    
    // The Fortran code builds C as:
    // For i-th node (i=0..3), C(i) has rows:
    // row0: H(i+1,1), H(i+1,2), H(i+1,3) -> but we need derivatives
    // Actually, from the Fortran: C(J+1) = H(I+4), etc.
    // Looking at the pattern: it's building gradient matrices
    
    // Standard tetrahedral shape function gradients
    // Shape functions: N1 = a1 + b1*x + c1*y + d1*z, etc.
    // From H matrix inversion, the coefficients are in H_inv
    // The gradient of Ni is [bi, ci, di] which is row i of H_inv (columns 1-3)
    
    for (int i = 0; i < 4; ++i) {
        C[i] = Eigen::Matrix<double, 6, 3>::Zero();
        // First three rows: gradient of Ni for ux, uy, uz components
        C[i].block<3,3>(0,0) = Eigen::Matrix3d::Identity() * H_inv(i,1);
        C[i].block<3,3>(3,0) = Eigen::Matrix3d::Identity() * H_inv(i,2);
        C[i].block<3,3>(6,0) = Eigen::Matrix3d::Identity() * H_inv(i,3);
        
        // But the Fortran pattern is different - let's reconstruct from the code:
        // C(J+1) = H(I+4) -> H(4), H(5), H(6), H(7) are x coords: nodes[0].x, nodes[1].x, nodes[2].x, nodes[3].x
        // Actually, looking more carefully at Fortran:
        // H(1)=1, H(2)=x1, H(3)=y1, H(4)=z1
        // H(5)=1, H(6)=x2, H(7)=y2, H(8)=z2
        // So H(I+4) for I=1..4 means H(5)..H(8) = x coords
        // H(I+8) = y coords, H(I+12) = z coords
        
        // So for node i (0-indexed), the gradient matrix should be:
        // [dxN_i, 0, 0]
        // [0, dxN_i, 0]
        // [0, 0, dxN_i]
        // [dyN_i, 0, 0]
        // [0, dyN_i, 0]
        // [0, 0, dyN_i]
        // [dzN_i, 0, 0]
        // [0, dzN_i, 0]
        // [0, 0, dzN_i]
        // But that doesn't match 6x3...
        
        // Standard approach: for tetrahedron, B matrix is 6x12, but we need 6x3 per node
        // Each node contributes a 6x3 block to the B matrix
        // B_i = [[dN_i/dx, 0, 0],
        //         [0, dN_i/dy, 0],
        //         [0, 0, dN_i/dz],
        //         [dN_i/dy, dN_i/dx, 0],
        //         [0, dN_i/dz, dN_i/dy],
        //         [dN_i/dz, 0, dN_i/dx]]
        
        // From H_inv, the shape function coefficients are in rows of H_inv
        // Ni = H_inv(i,0) + H_inv(i,1)*x + H_inv(i,2)*y + H_inv(i,3)*z
        // So dN_i/dx = H_inv(i,1), dN_i/dy = H_inv(i,2), dN_i/dz = H_inv(i,3)
        
        double dNdx = H_inv(i,1);
        double dNdy = H_inv(i,2);
        double dNdz = H_inv(i,3);
        
        // Build 6x3 B_i matrix
        C[i] << dNdx, 0.0, 0.0,
                0.0, dNdy, 0.0,
                0.0, 0.0, dNdz,
                dNdy, dNdx, 0.0,
                0.0, dNdz, dNdy,
                dNdz, 0.0, dNdx;
    }
    
    // Scale by determinant/6.0
    double scale = HDETER / 6.0;
    
    // Assemble global stiffness matrix K (12x12)
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    
    // For each pair of nodes i,j (0-3), add contribution C_i^T * G * C_j
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // C_i is 6x3, G is 6x6, C_j is 6x3
            // So C_i^T * G * C_j is 3x3
            Eigen::Matrix3d contrib = C[i].transpose() * G * C[j];
            
            // Place 3x3 contribution in appropriate position in 12x12 matrix
            // Node i: DOFs 3*i, 3*i+1, 3*i+2 (ux,uy,uz)
            // Node j: DOFs 3*j, 3*j+1, 3*j+2
            for (int ii = 0; ii < 3; ++ii) {
                for (int jj = 0; jj < 3; ++jj) {
                    K(3*i + ii, 3*j + jj) += scale * contrib(ii, jj);
                }
            }
        }
    }
    
    return K;
}

// Convert matrix to JSON string
std::string matrixToJson(const Eigen::Matrix<double, 12, 12>& mat) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(12) << mat(i, j);
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
        Eigen::Matrix<double, 12, 12> K = computeKTETRA(nodes, E, nu);
        
        // Output as JSON
        std::cout << matrixToJson(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}