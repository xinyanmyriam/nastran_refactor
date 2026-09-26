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
using Matrix6x6 = Eigen::Matrix<Real, 6, 6>;
using Matrix4x4 = Eigen::Matrix<Real, 4, 4>;
using Vector3 = Eigen::Vector3d;
using Matrix3x3 = Eigen::Matrix3d;
using Matrix6x3 = Eigen::Matrix<Real, 6, 3>;
using Matrix3x6 = Eigen::Matrix<Real, 3, 6>;

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

// Invert a 4x4 matrix and return determinant
bool invert4x4(const Matrix4x4& A, Matrix4x4& A_inv, Real& det) {
    // Compute determinant using cofactor expansion (simplified for 4x4)
    // We'll use Eigen's built-in inverse for robustness
    det = A.determinant();
    if (std::abs(det) < 1e-15) {
        return false;
    }
    A_inv = A.inverse();
    return true;
}

// Compute stiffness matrix for a tetrahedron element
Matrix6x6 computeTetraStiffness(
    const std::vector<Vector3>& nodes,
    Real E,
    Real nu,
    Real volume_factor = 1.0
) {
    // Nodes: [v1, v2, v3, v4] where v1,v2,v3 form base triangle, v4 is apex
    // Build H matrix: [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Matrix4x4 H;
    H << 1.0, nodes[0](0), nodes[0](1), nodes[0](2),
         1.0, nodes[1](0), nodes[1](1), nodes[1](2),
         1.0, nodes[2](0), nodes[2](1), nodes[2](2),
         1.0, nodes[3](0), nodes[3](1), nodes[3](2);

    Real det_H;
    Matrix4x4 H_inv;
    if (!invert4x4(H, H_inv, det_H)) {
        throw std::runtime_error("Singular tetrahedron geometry");
    }

    // Material matrix G (6x6) for isotropic material
    // G = [C11 C12 C12 0   0   0;
    //      C12 C11 C12 0   0   0;
    //      C12 C12 C11 0   0   0;
    //      0   0   0   C44 0   0;
    //      0   0   0   0   C44 0;
    //      0   0   0   0   0   C44]
    // where C11 = E*(1-nu)/((1+nu)*(1-2*nu)), C12 = E*nu/((1+nu)*(1-2*nu)), C44 = E/(2*(1+nu))
    Real denom = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(denom) < 1e-15) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    Real C11 = E * (1.0 - nu) / denom;
    Real C12 = E * nu / denom;
    Real C44 = E / (2.0 * (1.0 + nu));

    Matrix6x6 G;
    G << C11, C12, C12, 0.0, 0.0, 0.0,
         C12, C11, C12, 0.0, 0.0, 0.0,
         C12, C12, C11, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, C44, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, C44, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, C44;

    // Build C matrices (6x3 for each node)
    // C_i = [dN_i/dx, dN_i/dy, dN_i/dz]^T where N_i are shape functions
    // For linear tetrahedron: N_i = a_i + b_i*x + c_i*y + d_i*z
    // The gradients are constant: [b_i, c_i, d_i]^T
    // From H_inv: the last 3 columns contain the coefficients for x,y,z
    // So grad(N_i) = H_inv.row(i).tail(3)
    Matrix6x3 C;
    C.setZero();
    for (int i = 0; i < 4; ++i) {
        // grad(N_i) = [b_i, c_i, d_i] = H_inv.row(i).segment(1,3)
        C.block<3,3>(i*3, 0) = Matrix3x3::Identity() * H_inv(i, 1);
        C.block<3,3>(i*3, 0) += Matrix3x3::Identity() * H_inv(i, 2);
        C.block<3,3>(i*3, 0) += Matrix3x3::Identity() * H_inv(i, 3);
        
        // Actually, simpler: grad(N_i) = H_inv.row(i).segment(1,3)
        // So C for node i is [grad(N_i); grad(N_i); grad(N_i)] but arranged as 6x3
        // Standard: C = [B1; B2; B3; B4] where Bi = [dNi/dx, dNi/dy, dNi/dz] (3x3 block)
        Vector3 grad_Ni = H_inv.row(i).segment(1, 3);
        C.block<3,3>(i*3, 0) = grad_Ni * Matrix3x3::Identity();
    }

    // Actually, correct construction: for linear tetrahedron, B matrix is:
    // B = [dN1/dx 0 0 dN2/dx 0 0 dN3/dx 0 0 dN4/dx 0 0;
    //      0 dN1/dy 0 0 dN2/dy 0 0 dN3/dy 0 0 dN4/dy 0;
    //      0 0 dN1/dz 0 0 dN2/dz 0 0 dN3/dz 0 0 dN4/dz;
    //      dN1/dy dN1/dx 0 dN2/dy dN2/dx 0 dN3/dy dN3/dx 0 dN4/dy dN4/dx 0;
    //      0 dN1/dz dN1/dy 0 dN2/dz dN2/dy 0 dN3/dz dN3/dy 0 dN4/dz dN4/dy;
    //      dN1/dz 0 dN1/dx dN2/dz 0 dN2/dx dN3/dz 0 dN3/dx dN4/dz 0 dN4/dx]
    // But for our purpose, we'll use the standard approach
    
    // Instead, let's construct the B matrix properly
    // B is 6x12, but we need to build K = B^T * D * B * dV
    // For linear tetrahedron, B is constant, so K = V * B^T * D * B
    
    // Volume of tetrahedron = |det([v2-v1, v3-v1, v4-v1])| / 6
    Real volume = std::abs((nodes[1]-nodes[0]).dot(
        (nodes[2]-nodes[0]).cross(nodes[3]-nodes[0]))) / 6.0;
    
    // Scale factor: for wedge, each tetra gets volume_factor * volume
    Real scaled_volume = volume * volume_factor;
    
    // Build B matrix (6x12)
    Eigen::Matrix<Real, 6, 12> B;
    B.setZero();
    
    // Gradients of shape functions
    Vector3 b1 = H_inv.row(0).segment(1, 3);
    Vector3 b2 = H_inv.row(1).segment(1, 3);
    Vector3 b3 = H_inv.row(2).segment(1, 3);
    Vector3 b4 = H_inv.row(3).segment(1, 3);
    
    // B matrix assembly
    // Row 0: [b1.x, 0, 0, b2.x, 0, 0, b3.x, 0, 0, b4.x, 0, 0]
    B.row(0) << b1(0), 0, 0, b2(0), 0, 0, b3(0), 0, 0, b4(0), 0, 0;
    // Row 1: [0, b1.y, 0, 0, b2.y, 0, 0, b3.y, 0, 0, b4.y, 0]
    B.row(1) << 0, b1(1), 0, 0, b2(1), 0, 0, b3(1), 0, 0, b4(1), 0;
    // Row 2: [0, 0, b1.z, 0, 0, b2.z, 0, 0, b3.z, 0, 0, b4.z]
    B.row(2) << 0, 0, b1(2), 0, 0, b2(2), 0, 0, b3(2), 0, 0, b4(2);
    // Row 3: [b1.y, b1.x, 0, b2.y, b2.x, 0, b3.y, b3.x, 0, b4.y, b4.x, 0]
    B.row(3) << b1(1), b1(0), 0, b2(1), b2(0), 0, b3(1), b3(0), 0, b4(1), b4(0), 0;
    // Row 4: [0, b1.z, b1.y, 0, b2.z, b2.y, 0, b3.z, b3.y, 0, b4.z, b4.y]
    B.row(4) << 0, b1(2), b1(1), 0, b2(2), b2(1), 0, b3(2), b3(1), 0, b4(2), b4(1);
    // Row 5: [b1.z, 0, b1.x, b2.z, 0, b2.x, b3.z, 0, b3.x, b4.z, 0, b4.x]
    B.row(5) << b1(2), 0, b1(0), b2(2), 0, b2(0), b3(2), 0, b3(0), b4(2), 0, b4(0);
    
    // Compute element stiffness: K = B^T * G * B * scaled_volume
    Eigen::Matrix<Real, 12, 12> K_elem = B.transpose() * G * B * scaled_volume;
    
    // Extract 6x6 blocks for the global 18x18 matrix
    // But for now, return the full 12x12 for assembly
    // However, the problem asks for 18x18 (6 nodes * 3 DOF)
    // So we need to map each tetra's 4 nodes to the global 6-node wedge
    
    // Return a 6x6 stiffness contribution for this tetra's connectivity
    // Actually, we'll assemble into global later
    return Matrix6x6::Zero(); // placeholder
}

// Main KSOLID implementation for wedge element
Matrix18x18 computeWedgeStiffness(
    const std::vector<Vector3>& nodes, // 6 nodes: 0,1,2,3,4,5
    Real E,
    Real nu
) {
    // Wedge node ordering: 
    // Base triangle: nodes 0,1,2 (counter-clockwise when viewed from node 3)
    // Top triangle: nodes 3,4,5 (counter-clockwise when viewed from node 0)
    // So: N1=0, N2=1, N3=2, N4=3, N5=4, N6=5
    
    // Wedge decomposition into 3 tetrahedrons (as per Fortran M matrix)
    // M(1): [1,2,3,4] -> nodes[0], nodes[1], nodes[2], nodes[3]
    // M(2): [1,2,3,5] -> nodes[0], nodes[1], nodes[2], nodes[4]
    // M(3): [1,2,3,6] -> nodes[0], nodes[1], nodes[2], nodes[5]
    // M(4): [1,4,5,6] -> nodes[0], nodes[3], nodes[4], nodes[5]
    // M(5): [2,4,5,6] -> nodes[1], nodes[3], nodes[4], nodes[5]
    // M(6): [3,4,5,6] -> nodes[2], nodes[3], nodes[4], nodes[5]
    // M(7): [2,1,4,6] -> nodes[1], nodes[0], nodes[3], nodes[5]
    // M(8): [2,3,4,6] -> nodes[1], nodes[2], nodes[3], nodes[5]
    // M(9): [1,3,4,5] -> nodes[0], nodes[2], nodes[3], nodes[4]
    // M(10): [2,3,4,5] -> nodes[1], nodes[2], nodes[3], nodes[4]
    // M(11): [3,1,5,6] -> nodes[2], nodes[0], nodes[4], nodes[5]
    // M(12): [2,1,5,6] -> nodes[1], nodes[0], nodes[4], nodes[5]
    
    // According to NASTRAN documentation, wedge uses 3 tetrahedrons for basic case
    // But the Fortran shows 12 configurations, however for standard wedge it's 3
    // Let's use the first 3 as primary decomposition
    
    // First, validate geometry
    Vector3 v1 = nodes[0];
    Vector3 v2 = nodes[1];
    Vector3 v3 = nodes[2];
    Vector3 v4 = nodes[3];
    Vector3 v5 = nodes[4];
    Vector3 v6 = nodes[5];
    
    // Base triangle normal
    Vector3 base_normal = (v2 - v1).cross(v3 - v1);
    // Top triangle normal  
    Vector3 top_normal = (v5 - v4).cross(v6 - v4);
    
    // Check if base and top are parallel (should be for prism)
    Real base_norm = base_normal.norm();
    Real top_norm = top_normal.norm();
    if (base_norm < 1e-15 || top_norm < 1e-15) {
        throw std::runtime_error("Degenerate wedge geometry");
    }
    
    // Volume check: should be positive
    Real vol_base = std::abs(base_normal.dot(v4 - v1)) / 6.0;
    Real vol_top = std::abs(top_normal.dot(v1 - v4)) / 6.0;
    
    // Standard decomposition: 3 tetrahedrons sharing base triangle (0,1,2)
    std::vector<std::vector<int>> tetra_connectivity = {
        {0, 1, 2, 3}, // 1,2,3,4
        {0, 1, 2, 4}, // 1,2,3,5  
        {0, 1, 2, 5}, // 1,2,3,6
        {0, 3, 4, 5}, // 1,4,5,6
        {1, 3, 4, 5}, // 2,4,5,6
        {2, 3, 4, 5}, // 3,4,5,6
        {1, 0, 3, 5}, // 2,1,4,6
        {1, 2, 3, 5}, // 2,3,4,6
        {0, 2, 3, 4}, // 1,3,4,5
        {1, 2, 3, 4}, // 2,3,4,5
        {2, 0, 4, 5}, // 3,1,5,6
        {1, 0, 4, 5}  // 2,1,5,6
    };
    
    // Pre-compute all tetrahedron volumes and determinants
    std::vector<Real> tetra_volumes;
    for (const auto& conn : tetra_connectivity) {
        Vector3 p1 = nodes[conn[0]];
        Vector3 p2 = nodes[conn[1]];
        Vector3 p3 = nodes[conn[2]];
        Vector3 p4 = nodes[conn[3]];
        Real vol = std::abs((p2-p1).dot((p3-p1).cross(p4-p1))) / 6.0;
        tetra_volumes.push_back(vol);
    }
    
    // Total volume should be positive
    Real total_vol = 0.0;
    for (Real v : tetra_volumes) total_vol += v;
    if (total_vol < 1e-15) {
        throw std::runtime_error("Zero volume wedge");
    }
    
    // Initialize global stiffness matrix
    Matrix18x18 K_global = Matrix18x18::Zero();
    
    // Material properties
    Real denom = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(denom) < 1e-15) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    Real C11 = E * (1.0 - nu) / denom;
    Real C12 = E * nu / denom;
    Real C44 = E / (2.0 * (1.0 + nu));
    
    Matrix6x6 G;
    G << C11, C12, C12, 0.0, 0.0, 0.0,
         C12, C11, C12, 0.0, 0.0, 0.0,
         C12, C12, C11, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, C44, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, C44, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, C44;
    
    // Process each tetrahedron
    for (size_t t = 0; t < tetra_connectivity.size(); ++t) {
        const auto& conn = tetra_connectivity[t];
        std::vector<Vector3> tet_nodes = {
            nodes[conn[0]], nodes[conn[1]], nodes[conn[2]], nodes[conn[3]]
        };
        
        // Build H matrix
        Matrix4x4 H;
        H << 1.0, tet_nodes[0](0), tet_nodes[0](1), tet_nodes[0](2),
             1.0, tet_nodes[1](0), tet_nodes[1](1), tet_nodes[1](2),
             1.0, tet_nodes[2](0), tet_nodes[2](1), tet_nodes[2](2),
             1.0, tet_nodes[3](0), tet_nodes[3](1), tet_nodes[3](2);
        
        Real det_H;
        Matrix4x4 H_inv;
        if (!invert4x4(H, H_inv, det_H)) {
            continue; // skip degenerate tet
        }
        
        // Volume of this tetrahedron
        Real vol = tetra_volumes[t];
        
        // Volume scaling factors as per Fortran logic:
        // For wedges: first 6 configurations multiplied by 2, all divided by 36
        Real volume_factor;
        if (t < 6) {
            volume_factor = vol * 2.0 / 36.0;
        } else {
            volume_factor = vol / 36.0;
        }
        
        // Build B matrix (6x12 for 4 nodes * 3 DOF)
        Eigen::Matrix<Real, 6, 12> B;
        B.setZero();
        
        // Gradients of shape functions
        Vector3 b1 = H_inv.row(0).segment(1, 3);
        Vector3 b2 = H_inv.row(1).segment(1, 3);
        Vector3 b3 = H_inv.row(2).segment(1, 3);
        Vector3 b4 = H_inv.row(3).segment(1, 3);
        
        // B matrix assembly (standard linear tetrahedron)
        B.row(0) << b1(0), 0, 0, b2(0), 0, 0, b3(0), 0, 0, b4(0), 0, 0;
        B.row(1) << 0, b1(1), 0, 0, b2(1), 0, 0, b3(1), 0, 0, b4(1), 0;
        B.row(2) << 0, 0, b1(2), 0, 0, b2(2), 0, 0, b3(2), 0, 0, b4(2);
        B.row(3) << b1(1), b1(0), 0, b2(1), b2(0), 0, b3(1), b3(0), 0, b4(1), b4(0), 0;
        B.row(4) << 0, b1(2), b1(1), 0, b2(2), b2(1), 0, b3(2), b3(1), 0, b4(2), b4(1);
        B.row(5) << b1(2), 0, b1(0), b2(2), 0, b2(0), b3(2), 0, b3(0), b4(2), 0, b4(0);
        
        // Element stiffness: K = B^T * G * B * volume_factor
        Eigen::Matrix<Real, 12, 12> K_tet = B.transpose() * G * B * volume_factor;
        
        // Map local DOFs to global DOFs
        // Local DOFs: node0->(0,1,2), node1->(3,4,5), node2->(6,7,8), node3->(9,10,11)
        // Global DOFs for wedge: node0->(0,1,2), node1->(3,4,5), node2->(6,7,8), 
        //                        node3->(9,10,11), node4->(12,13,14), node5->(15,16,17)
        std::vector<int> global_dof_map = {
            conn[0]*3, conn[0]*3+1, conn[0]*3+2,
            conn[1]*3, conn[1]*3+1, conn[1]*3+2,
            conn[2]*3, conn[2]*3+1, conn[2]*3+2,
            conn[3]*3, conn[3]*3+1, conn[3]*3+2
        };
        
        // Add contribution to global matrix
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                int gi = global_dof_map[i];
                int gj = global_dof_map[j];
                if (gi < 18 && gj < 18) {
                    K_global(gi, gj) += K_tet(i, j);
                }
            }
        }
    }
    
    return K_global;
}

// Main function with test case
int main() {
    try {
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
        
        // Compute stiffness matrix
        Matrix18x18 K = computeWedgeStiffness(nodes, E, nu);
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 18; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << "[";
            for (int j = 0; j < 18; ++j) {
                if (j > 0) std::cout << ",";
                // Format in scientific notation with sufficient precision
                std::cout << std::scientific << std::setprecision(15) << K(i, j);
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