#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Since we can't use external JSON library in strict environment, we'll implement minimal JSON serialization
// But note: the problem says "JSON must be valid and parseable", and we need to output exactly the format.
// We'll create a simple JSON string builder.

// Minimal JSON string builder for array of arrays
std::string to_json_string(const Eigen::MatrixXd& mat) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < mat.rows(); ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < mat.cols(); ++j) {
            if (j > 0) oss << ",";
            // Use scientific notation with 15 digits precision
            oss << std::scientific << std::setprecision(15) << mat(i,j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

// Cross product: r = a × b
Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.cross(b);
}

// Dot product
double dot(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.dot(b);
}

// Invert 4x4 matrix and compute determinant
bool invert_4x4(const Eigen::Matrix4d& H, Eigen::Matrix4d& H_inv, double& det) {
    det = H.determinant();
    if (std::abs(det) < 1e-15) {
        return false;
    }
    H_inv = H.inverse();
    return true;
}

// Compute stiffness matrix for a tetrahedron element
Eigen::MatrixXd tetrahedron_stiffness(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    // H matrix: [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix4d H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);
    
    Eigen::Matrix4d H_inv;
    double det_H;
    if (!invert_4x4(H, H_inv, det_H)) {
        throw std::runtime_error("Singular tetrahedron geometry");
    }
    
    // Material matrix G (6x6) in Voigt notation
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    
    double G11 = E * (1.0 - nu) / temp1;
    double G12 = E * nu / temp1;
    double GG = E / (2.0 * (1.0 + nu)); // Shear modulus
    
    // G = [ C11 C12 C12  0    0    0  ]
    //     [ C12 C11 C12  0    0    0  ]
    //     [ C12 C12 C11  0    0    0  ]
    //     [ 0   0   0   C44  0    0  ]
    //     [ 0   0   0    0  C55   0  ]
    //     [ 0   0   0    0   0   C66 ]
    // where C44 = C55 = C66 = GG
    Eigen::Matrix<double, 6, 6> G;
    G.setZero();
    G(0,0) = G11; G(1,1) = G11; G(2,2) = G11;
    G(0,1) = G12; G(0,2) = G12;
    G(1,0) = G12; G(1,2) = G12;
    G(2,0) = G12; G(2,1) = G12;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;
    
    // B matrices (4 matrices of size 6x3) — strain-displacement matrices in Voigt notation
    // For node i, gradient = [b_i, c_i, d_i] = row i of H_inv.block(0,1,4,3)
    std::vector<Eigen::Matrix<double, 6, 3>> B_matrices;
    for (int i = 0; i < 4; ++i) {
        double b_i = H_inv(i, 1); // dN_i/dx
        double c_i = H_inv(i, 2); // dN_i/dy
        double d_i = H_inv(i, 3); // dN_i/dz

        Eigen::Matrix<double, 6, 3> B_i;
        B_i.setZero();
        // ε_xx = ∂u/∂x → b_i * u
        B_i(0,0) = b_i;
        // ε_yy = ∂v/∂y → c_i * v
        B_i(1,1) = c_i;
        // ε_zz = ∂w/∂z → d_i * w
        B_i(2,2) = d_i;
        // γ_yz = ∂v/∂z + ∂w/∂y → c_i * w + d_i * v
        B_i(3,1) = d_i;
        B_i(3,2) = c_i;
        // γ_xz = ∂u/∂z + ∂w/∂x → b_i * w + d_i * u
        B_i(4,0) = d_i;
        B_i(4,2) = b_i;
        // γ_xy = ∂u/∂y + ∂v/∂x → b_i * v + c_i * u
        B_i(5,0) = c_i;
        B_i(5,1) = b_i;
        B_matrices.push_back(B_i);
    }
    
    // Build global 12x12 stiffness matrix for tetrahedron (4 nodes * 3 DOF)
    Eigen::MatrixXd K_tetra = Eigen::MatrixXd::Zero(12, 12);
    
    // Volume factor: |det(H)|/6.0
    double volume_factor = std::abs(det_H) / 6.0;
    
    // For each pair of nodes i,j: K += B_i^T * G * B_j * volume_factor
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Eigen::Matrix<double, 3, 3> contrib = B_matrices[i].transpose() * G * B_matrices[j];
            
            // Place in global matrix: rows [3*i, 3*i+2], cols [3*j, 3*j+2]
            for (int di = 0; di < 3; ++di) {
                for (int dj = 0; dj < 3; ++dj) {
                    K_tetra(3*i + di, 3*j + dj) += contrib(di, dj) * volume_factor;
                }
            }
        }
    }
    
    return K_tetra;
}

// Assemble wedge stiffness matrix using 12 tetrahedrons decomposition
Eigen::MatrixXd wedge_stiffness(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    // Wedge nodes: 0,1,2,3,4,5 correspond to N1,N2,N3,N4,N5,N6
    // N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    
    // Define the 12 tetrahedrons for wedge decomposition
    // Node indexing: 0=N1, 1=N2, 2=N3, 3=N4, 4=N5, 5=N6
    std::vector<std::vector<int>> tetrahedrons = {
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
    
    // Global stiffness matrix: 6 nodes * 3 DOF = 18x18
    Eigen::MatrixXd K_global = Eigen::MatrixXd::Zero(18, 18);
    
    // For each tetrahedron
    for (const auto& tet : tetrahedrons) {
        // Get node coordinates for this tetrahedron
        std::vector<Eigen::Vector3d> tet_nodes;
        for (int idx : tet) {
            tet_nodes.push_back(nodes[idx]);
        }
        
        // Compute tetrahedron stiffness
        Eigen::MatrixXd K_tet = tetrahedron_stiffness(tet_nodes, E, nu);
        
        // Map local DOFs (12) to global DOFs (18)
        // Each node has 3 DOFs: ux,uy,uz
        std::vector<int> global_dofs;
        for (int node_idx : tet) {
            global_dofs.push_back(3*node_idx);     // ux
            global_dofs.push_back(3*node_idx+1);   // uy
            global_dofs.push_back(3*node_idx+2);   // uz
        }
        
        // Add contribution to global matrix
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                int gi = global_dofs[i];
                int gj = global_dofs[j];
                K_global(gi, gj) += K_tet(i, j);
            }
        }
    }
    
    return K_global;
}

int main() {
    try {
        // Test case: Wedge (6 nodes): N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
        std::vector<Eigen::Vector3d> nodes = {
            Eigen::Vector3d(0.0, 0.0, 0.0), // N1
            Eigen::Vector3d(1.0, 0.0, 0.0), // N2
            Eigen::Vector3d(0.0, 1.0, 0.0), // N3
            Eigen::Vector3d(0.0, 0.0, 1.0), // N4
            Eigen::Vector3d(1.0, 0.0, 1.0), // N5
            Eigen::Vector3d(0.0, 1.0, 1.0)  // N6
        };
        
        // Material properties
        double E = 200e9;      // 200 GPa
        double nu = 0.3;
        
        // Compute wedge stiffness matrix
        Eigen::MatrixXd K = wedge_stiffness(nodes, E, nu);
        
        // Output as JSON
        std::cout << to_json_string(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}