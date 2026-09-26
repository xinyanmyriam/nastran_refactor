#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Since we cannot use external JSON library in strict MSVC compilation without linking,
// we'll implement a minimal JSON serializer for the required format.
// We'll use a simple string-based approach to generate valid JSON.

std::string double_to_string(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos) {
        // Remove trailing zeros after decimal
        size_t dot = s.find('.');
        if (dot != std::string::npos && dot < epos) {
            while (s.back() == '0' && s.size() > dot + 1) {
                s.pop_back();
            }
            if (s.back() == '.') {
                s.pop_back();
            }
        }
        // Normalize exponent sign: ensure + is present or remove extra +
        std::string exp = s.substr(epos);
        if (exp.length() >= 3 && exp[1] == '+' && exp[2] == '0') {
            // Replace "+0" with "0" but keep at least one digit
            if (exp.length() >= 4 && exp[3] != '\0') {
                s.replace(epos+1, 2, "");
            }
        }
    }
    return s;
}

std::string matrix_to_json(const Eigen::MatrixXd& mat) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    for (int i = 0; i < mat.rows(); ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < mat.cols(); ++j) {
            if (j > 0) oss << ",";
            oss << double_to_string(mat(i,j));
        }
        oss << "]";
    }
    oss << "]}";
    return oss.str();
}

// Minimal cross product implementation
Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return Eigen::Vector3d(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
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

// Tetrahedron stiffness contribution
Eigen::Matrix<double, 12, 12> tetrahedron_stiffness(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    // Nodes are 4 points: p1, p2, p3, p4
    // Build H matrix: [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix4d H;
    H.row(0) << 1.0, nodes[0](0), nodes[0](1), nodes[0](2);
    H.row(1) << 1.0, nodes[1](0), nodes[1](1), nodes[1](2);
    H.row(2) << 1.0, nodes[2](0), nodes[2](1), nodes[2](2);
    H.row(3) << 1.0, nodes[3](0), nodes[3](1), nodes[3](2);
    
    double det_H;
    Eigen::Matrix4d H_inv;
    if (!invert_4x4(H, H_inv, det_H)) {
        throw std::runtime_error("Singular tetrahedron geometry");
    }
    
    // Material matrix G (6x6) for isotropic material
    // G = [C11 C12 C12 0   0   0  ;
    //      C12 C11 C12 0   0   0  ;
    //      C12 C12 C11 0   0   0  ;
    //      0   0   0   C44 0   0  ;
    //      0   0   0   0   C44 0  ;
    //      0   0   0   0   0   C44]
    // where C11 = E*(1-nu)/((1+nu)*(1-2*nu)), C12 = E*nu/((1+nu)*(1-2*nu)), C44 = G = E/(2*(1+nu))
    double denom = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(denom) < 1e-15) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    double C11 = E * (1.0 - nu) / denom;
    double C12 = E * nu / denom;
    double C44 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 6, 6> G;
    G.setZero();
    G(0,0) = C11; G(0,1) = C12; G(0,2) = C12;
    G(1,0) = C12; G(1,1) = C11; G(1,2) = C12;
    G(2,0) = C12; G(2,1) = C12; G(2,2) = C11;
    G(3,3) = C44;
    G(4,4) = C44;
    G(5,5) = C44;
    
    // C matrices: 4 matrices of size 6x3
    // For each node i (0-indexed), C_i = [dN_i/dx, dN_i/dy, dN_i/dz] as 6x3
    // where N_i are shape functions: N_i = H_inv_row_i * [1,x,y,z]^T
    // So dN_i/dx = H_inv(i,1), dN_i/dy = H_inv(i,2), dN_i/dz = H_inv(i,3)
    std::vector<Eigen::Matrix<double, 6, 3>> C(4);
    for (int i = 0; i < 4; ++i) {
        C[i].setZero();
        double dNdx = H_inv(i,1);
        double dNdy = H_inv(i,2);
        double dNdz = H_inv(i,3);
        
        // Strain-displacement matrix for node i (6x3):
        // Row 0 (εxx): [dNdx, 0, 0]
        // Row 1 (εyy): [0, dNdy, 0]
        // Row 2 (εzz): [0, 0, dNdz]
        // Row 3 (γxy): [dNdy, dNdx, 0]
        // Row 4 (γyz): [0, dNdz, dNdy]
        // Row 5 (γzx): [dNdz, 0, dNdx]
        C[i](0,0) = dNdx;
        C[i](1,1) = dNdy;
        C[i](2,2) = dNdz;
        C[i](3,0) = dNdy;
        C[i](3,1) = dNdx;
        C[i](4,1) = dNdz;
        C[i](4,2) = dNdy;
        C[i](5,0) = dNdz;
        C[i](5,2) = dNdx;
    }
    
    // Volume factor: |det_H|/6.0
    double vol_factor = std::abs(det_H) / 6.0;
    
    // Assemble 12x12 stiffness matrix K = sum_{i,j} C_i^T * G * C_j * vol_factor
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // C_i^T * G * C_j is 3x3 block for each DOF pair
            // C_i is 6x3, G is 6x6, C_j is 6x3 -> C_i^T * G * C_j is 3x3
            Eigen::Matrix<double, 3, 3> block = C[i].transpose() * G * C[j];
            
            // Place block in K at positions (3*i, 3*j)
            for (int di = 0; di < 3; ++di) {
                for (int dj = 0; dj < 3; ++dj) {
                    K(3*i + di, 3*j + dj) += block(di, dj) * vol_factor;
                }
            }
        }
    }
    
    return K;
}

// Wedge element: 6 nodes, decomposed into 3 tetrahedrons as per NASTRAN mapping
// Mapping from Fortran: M(1,:) = {1,2,3,4}, M(2,:) = {1,2,3,5}, M(3,:) = {1,2,3,6}
Eigen::Matrix<double, 18, 18> wedge_stiffness(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    if (nodes.size() != 6) {
        throw std::runtime_error("Wedge must have exactly 6 nodes");
    }
    
    // Node indexing: 0-based: N0,N1,N2,N3,N4,N5 correspond to N1,N2,N3,N4,N5,N6 in spec
    // Spec: N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    
    // Create 3 tetrahedrons:
    // Tet1: nodes[0], nodes[1], nodes[2], nodes[3]  // 1,2,3,4
    // Tet2: nodes[0], nodes[1], nodes[2], nodes[4]  // 1,2,3,5
    // Tet3: nodes[0], nodes[1], nodes[2], nodes[5]  // 1,2,3,6
    
    std::vector<std::vector<Eigen::Vector3d>> tets;
    tets.push_back({nodes[0], nodes[1], nodes[2], nodes[3]});
    tets.push_back({nodes[0], nodes[1], nodes[2], nodes[4]});
    tets.push_back({nodes[0], nodes[1], nodes[2], nodes[5]});
    
    // Initialize 18x18 stiffness matrix
    Eigen::Matrix<double, 18, 18> K_total = Eigen::Matrix<double, 18, 18>::Zero();
    
    // Process each tetrahedron
    for (size_t t = 0; t < tets.size(); ++t) {
        Eigen::Matrix<double, 12, 12> K_tet = tetrahedron_stiffness(tets[t], E, nu);
        
        // Map local DOFs (12) to global DOFs (18)
        // Tet1: nodes 0,1,2,3 -> DOFs 0-2, 3-5, 6-8, 9-11
        // Tet2: nodes 0,1,2,4 -> DOFs 0-2, 3-5, 6-8, 12-14
        // Tet3: nodes 0,1,2,5 -> DOFs 0-2, 3-5, 6-8, 15-17
        
        std::vector<int> local_to_global;
        if (t == 0) {
            // nodes 0,1,2,3 -> DOFs 0-2, 3-5, 6-8, 9-11
            local_to_global = {0,1,2, 3,4,5, 6,7,8, 9,10,11};
        } else if (t == 1) {
            // nodes 0,1,2,4 -> DOFs 0-2, 3-5, 6-8, 12-14
            local_to_global = {0,1,2, 3,4,5, 6,7,8, 12,13,14};
        } else if (t == 2) {
            // nodes 0,1,2,5 -> DOFs 0-2, 3-5, 6-8, 15-17
            local_to_global = {0,1,2, 3,4,5, 6,7,8, 15,16,17};
        }
        
        // Add contribution
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                int gi = local_to_global[i];
                int gj = local_to_global[j];
                K_total(gi, gj) += K_tet(i, j);
            }
        }
    }
    
    return K_total;
}

int main() {
    // Test case: Wedge (6 nodes): N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    std::vector<Eigen::Vector3d> nodes;
    nodes.push_back(Eigen::Vector3d(0.0, 0.0, 0.0)); // N1
    nodes.push_back(Eigen::Vector3d(1.0, 0.0, 0.0)); // N2
    nodes.push_back(Eigen::Vector3d(0.0, 1.0, 0.0)); // N3
    nodes.push_back(Eigen::Vector3d(0.0, 0.0, 1.0)); // N4
    nodes.push_back(Eigen::Vector3d(1.0, 0.0, 1.0)); // N5
    nodes.push_back(Eigen::Vector3d(0.0, 1.0, 1.0)); // N6
    
    double E = 200e9;     // Pa
    double nu = 0.3;
    
    try {
        Eigen::Matrix<double, 18, 18> K = wedge_stiffness(nodes, E, nu);
        
        // Output as JSON
        std::cout << matrix_to_json(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}