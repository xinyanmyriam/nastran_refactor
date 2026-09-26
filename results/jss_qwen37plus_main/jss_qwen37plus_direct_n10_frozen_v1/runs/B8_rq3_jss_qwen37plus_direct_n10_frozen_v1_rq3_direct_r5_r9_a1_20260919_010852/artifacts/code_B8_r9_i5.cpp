#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Use double precision throughout
using Real = double;
using Matrix18x18 = Eigen::Matrix<Real, 18, 18>;
using Vector3 = Eigen::Vector3d;
using Matrix3x3 = Eigen::Matrix3d;
using Matrix6x6 = Eigen::Matrix<Real, 6, 6>;
using Matrix3x6 = Eigen::Matrix<Real, 3, 6>;
using Matrix6x3 = Eigen::Matrix<Real, 6, 3>;

// Helper function to compute cross product
Vector3 cross(const Vector3& a, const Vector3& b) {
    return Vector3(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
}

// Helper function to compute dot product
Real dot(const Vector3& a, const Vector3& b) {
    return a.dot(b);
}

// Invert a 4x4 matrix and compute determinant
bool invert4x4(const Eigen::Matrix<Real, 4, 4>& A, Eigen::Matrix<Real, 4, 4>& A_inv, Real& det) {
    // Compute determinant using LU decomposition for stability
    Eigen::FullPivLU<Eigen::Matrix<Real, 4, 4>> lu(A);
    if (lu.rank() < 4) {
        return false;
    }
    det = lu.determinant();
    A_inv = lu.inverse();
    return true;
}

// Compute stiffness matrix for a tetrahedron element
Matrix18x18 computeTetraStiffness(
    const std::vector<Vector3>& nodes,
    Real E, Real nu,
    Real volume_factor = 1.0
) {
    // Nodes: [v1, v2, v3, v4] where each is (x,y,z)
    // H matrix: [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix<Real, 4, 4> H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);

    Real det_H;
    Eigen::Matrix<Real, 4, 4> H_inv;
    if (!invert4x4(H, H_inv, det_H)) {
        throw std::runtime_error("Singular tetrahedron geometry");
    }

    // Ensure positive volume
    if (det_H <= 0.0) {
        throw std::runtime_error("Invalid tetrahedron orientation");
    }

    // Material matrix G (6x6) for isotropic material
    // G = [C11 C12 C12 0   0   0;
    //      C12 C11 C12 0   0   0;
    //      C12 C12 C11 0   0   0;
    //      0   0   0   C44 0   0;
    //      0   0   0   0   C44 0;
    //      0   0   0   0   0   C44]
    // where C11 = E*(1-nu)/((1+nu)*(1-2*nu)), C12 = E*nu/((1+nu)*(1-2*nu)), C44 = G = E/(2*(1+nu))
    Real temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    Real C11 = E * (1.0 - nu) / temp1;
    Real C12 = E * nu / temp1;
    Real C44 = E / (2.0 * (1.0 + nu));

    Matrix6x6 G;
    G << C11, C12, C12, 0.0, 0.0, 0.0,
         C12, C11, C12, 0.0, 0.0, 0.0,
         C12, C12, C11, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, C44, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, C44, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, C44;

    // C matrices (4 matrices of size 6x3)
    // For each node i, C_i = [dN_i/dx, dN_i/dy, dN_i/dz] where N_i are shape functions
    // dN_i/dx, dN_i/dy, dN_i/dz are rows of H_inv transpose without first column
    std::vector<Matrix6x3> C(4);
    for (int i = 0; i < 4; ++i) {
        // Extract row i from H_inv (which gives coefficients for shape function N_i)
        // Shape function: N_i = H_inv_row_i * [1, x, y, z]^T
        // So dN_i/dx = H_inv(i,1), dN_i/dy = H_inv(i,2), dN_i/dz = H_inv(i,3)
        Real dNdx = H_inv(i, 1);
        Real dNdy = H_inv(i, 2);
        Real dNdz = H_inv(i, 3);
        
        // C_i is 6x3: for DOFs u,v,w at node i, the derivatives are:
        // [dN_i/dx, 0, 0;
        //  0, dN_i/dy, 0;
        //  0, 0, dN_i/dz;
        //  dN_i/dy, dN_i/dx, 0;
        //  0, dN_i/dz, dN_i/dy;
        //  dN_i/dz, 0, dN_i/dx]
        // But standard linear tetrahedron uses simpler form: only diagonal terms for displacement gradients
        // Standard: B matrix has rows: [dN_i/dx, 0, 0, dN_j/dx, 0, 0, ...] but we need per-node contribution
        
        // Actually, for stiffness assembly, we need the strain-displacement matrix B
        // For tetrahedron, B is 6x12, but we'll build contributions per node pair
        // Instead, construct C_i as 3x3 for displacement gradients at node i
        // But the Fortran code builds C as 6x3 per node with specific pattern
        // From Fortran: C(J+1) = H(I+4), C(J+5) = H(I+8), C(J+9) = H(I+12), etc.
        // This corresponds to taking rows from H_inv
        // Let's reconstruct the pattern used in Fortran:
        // For node i (0-indexed), use H_inv(i,1), H_inv(i,2), H_inv(i,3) for x,y,z derivatives
        // The Fortran pattern creates a 6x3 matrix where:
        // row0: [H_inv(i,1), 0, 0]
        // row1: [0, H_inv(i,2), 0]
        // row2: [0, 0, H_inv(i,3)]
        // row3: [H_inv(i,2), H_inv(i,1), 0]
        // row4: [0, H_inv(i,3), H_inv(i,2)]
        // row5: [H_inv(i,3), 0, H_inv(i,1)]
        
        Matrix6x3 Ci;
        Ci << dNdx, 0.0, 0.0,
              0.0, dNdy, 0.0,
              0.0, 0.0, dNdz,
              dNdy, dNdx, 0.0,
              0.0, dNdz, dNdy,
              dNdz, 0.0, dNdx;
        C[i] = Ci;
    }

    // Volume factor: divide by 6.0, and additional factors per case
    Real vol_factor_total = std::abs(det_H) / 6.0 * volume_factor;

    // Assemble global stiffness matrix (18x18)
    Matrix18x18 K = Matrix18x18::Zero();

    // For each pair of nodes (i,j), compute contribution K_ij = vol_factor * C_i^T * G * C_j
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // C_i^T * G * C_j gives 3x3 block for DOF pairing
            // But C_i is 6x3, so C_i^T is 3x6, G is 6x6, C_j is 6x3
            // Result is 3x3 matrix for the coupling between node i and node j
            Matrix3x3 K_block = C[i].transpose() * G * C[j];
            
            // Map to global DOFs: node i has DOFs [3*i, 3*i+1, 3*i+2], same for j
            for (int di = 0; di < 3; ++di) {
                for (int dj = 0; dj < 3; ++dj) {
                    int global_i = 3 * i + di;
                    int global_j = 3 * j + dj;
                    K(global_i, global_j) += vol_factor_total * K_block(di, dj);
                }
            }
        }
    }

    return K;
}

// Main KSOLID implementation for wedge element (6 nodes)
Matrix18x18 KSOLID_WEDGE(
    const std::vector<Vector3>& nodes,
    Real E, Real nu
) {
    // Wedge nodes: N1,N2,N3,N4,N5,N6
    // Mapping to tetrahedrons (12 tets as per Fortran)
    // M(1): 1,2,3,4
    // M(2): 1,2,3,5
    // M(3): 1,2,3,6
    // M(4): 1,4,5,6
    // M(5): 2,4,5,6
    // M(6): 3,4,5,6
    // M(7): 2,1,4,6
    // M(8): 2,3,4,6
    // M(9): 1,3,4,5
    // M(10): 2,3,4,5
    // M(11): 3,1,5,6
    // M(12): 2,1,5,6
    
    // Define the 12 tetrahedrons (0-indexed node indices)
    std::vector<std::vector<int>> tets = {
        {0,1,2,3}, // 1,2,3,4
        {0,1,2,4}, // 1,2,3,5
        {0,1,2,5}, // 1,2,3,6
        {0,3,4,5}, // 1,4,5,6
        {1,3,4,5}, // 2,4,5,6
        {2,3,4,5}, // 3,4,5,6
        {1,0,3,5}, // 2,1,4,6
        {1,2,3,5}, // 2,3,4,6
        {0,2,3,4}, // 1,3,4,5
        {1,2,3,4}, // 2,3,4,5
        {2,0,4,5}, // 3,1,5,6
        {1,0,4,5}  // 2,1,5,6
    };

    // Volume factors for each tet
    // From Fortran: for wedges, first 6 configs multiplied by 2, all divided by 36
    // So factors: [2/36, 2/36, 2/36, 2/36, 2/36, 2/36, 1/36, 1/36, 1/36, 1/36, 1/36, 1/36]
    std::vector<Real> factors(12, 1.0/36.0);
    for (int i = 0; i < 6; ++i) {
        factors[i] *= 2.0;
    }

    Matrix18x18 K_total = Matrix18x18::Zero();

    // Check geometry: ensure base triangle 1,2,3 and top triangle 4,5,6 are properly oriented
    Vector3 v12 = nodes[1] - nodes[0];
    Vector3 v13 = nodes[2] - nodes[0];
    Vector3 normal_base = cross(v12, v13);
    
    Vector3 v45 = nodes[4] - nodes[3];
    Vector3 v46 = nodes[5] - nodes[3];
    Vector3 normal_top = cross(v45, v46);
    
    // Check if normals point in same general direction
    if (dot(normal_base, normal_top) < 0) {
        // Flip top triangle orientation by swapping nodes 5 and 6
        std::swap(nodes[4], nodes[5]);
    }

    // Process each tetrahedron
    for (int i = 0; i < 12; ++i) {
        std::vector<Vector3> tet_nodes;
        for (int j : tets[i]) {
            tet_nodes.push_back(nodes[j]);
        }
        
        try {
            Matrix18x18 K_tet = computeTetraStiffness(tet_nodes, E, nu, factors[i]);
            K_total += K_tet;
        } catch (const std::exception& e) {
            // Skip problematic tets (shouldn't happen for valid wedge)
            continue;
        }
    }

    return K_total;
}

// Convert matrix to JSON string
std::string matrixToJson(const Matrix18x18& K) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 18; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 18; ++j) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(12) << K(i, j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    // Test case: Wedge (6 nodes): N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    std::vector<Vector3> nodes = {
        Vector3(0.0, 0.0, 0.0), // N1
        Vector3(1.0, 0.0, 0.0), // N2
        Vector3(0.0, 1.0, 0.0), // N3
        Vector3(0.0, 0.0, 1.0), // N4
        Vector3(1.0, 0.0, 1.0), // N5
        Vector3(0.0, 1.0, 1.0)  // N6
    };

    Real E = 200e9;     // 200 GPa
    Real nu = 0.3;

    try {
        Matrix18x18 K = KSOLID_WEDGE(nodes, E, nu);
        
        // Output as JSON
        std::string json = matrixToJson(K);
        std::cout << json << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}