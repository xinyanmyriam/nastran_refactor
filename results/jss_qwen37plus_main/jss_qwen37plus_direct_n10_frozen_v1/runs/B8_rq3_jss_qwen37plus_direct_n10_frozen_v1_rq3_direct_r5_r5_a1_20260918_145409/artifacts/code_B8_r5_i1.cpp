#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <algorithm>

// Helper function to compute cross product of two 3D vectors
Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return Eigen::Vector3d(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
}

// Helper function to compute dot product of two 3D vectors
double dot(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.dot(b);
}

// Helper function to compute determinant of 4x4 matrix (for tetrahedron)
double det4x4(const Eigen::Matrix4d& m) {
    return m.determinant();
}

// Helper function to invert 4x4 matrix and get determinant
bool invert4x4(const Eigen::Matrix4d& m, Eigen::Matrix4d& inv, double& det) {
    det = det4x4(m);
    if (std::abs(det) < 1e-15) {
        return false;
    }
    inv = m.inverse();
    return true;
}

// Tetrahedron stiffness computation
Eigen::Matrix<double, 12, 12> computeTetStiffness(
    const std::vector<Eigen::Vector3d>& coords,
    double E, double nu) {
    
    // coords[0], coords[1], coords[2], coords[3] are the 4 nodes
    
    // Build H matrix: [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix4d H;
    H << 1.0, coords[0](0), coords[0](1), coords[0](2),
         1.0, coords[1](0), coords[1](1), coords[1](2),
         1.0, coords[2](0), coords[2](1), coords[2](2),
         1.0, coords[3](0), coords[3](1), coords[3](2);
    
    double detH;
    Eigen::Matrix4d Hinv;
    if (!invert4x4(H, Hinv, detH)) {
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
    double denom = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(denom) < 1e-15) {
        throw std::runtime_error("Invalid Poisson's ratio");
    }
    double C11 = E * (1.0 - nu) / denom;
    double C12 = E * nu / denom;
    double C44 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 6, 6> G;
    G << C11, C12, C12, 0.0, 0.0, 0.0,
         C12, C11, C12, 0.0, 0.0, 0.0,
         C12, C12, C11, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, C44, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, C44, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, C44;
    
    // Build C matrices (4 matrices of size 6x3 each)
    // For each node i, Ci = [bi 0 0; 0 bi 0; 0 0 bi] where bi is row i of Hinv (without first column)
    // Actually, standard formulation: Ci = [bi_x bi_y bi_z 0 0 0; 0 0 0 bi_x bi_y bi_z] but we need 6x3
    // Standard: for node i, Ci is 6x3: [bi_x 0 0; bi_y 0 0; bi_z 0 0; 0 bi_x 0; 0 bi_y 0; 0 bi_z 0] 
    // But NASTRAN uses different convention. From the Fortran, it appears they use:
    // C(J+1)=H(I+4), C(J+5)=H(I+8), C(J+9)=H(I+12), etc.
    // This suggests they're using Hinv rows, but let's use standard B-matrix approach.
    
    // Standard B-matrix for tetrahedron: B = [b1 0 0; 0 b1 0; 0 0 b1; b2 0 0; ...] but assembled differently
    // Actually, for displacement u = [u1 v1 w1 u2 v2 w2 u3 v3 w3 u4 v4 w4]^T
    // strain = B * u, where B is 6x12
    // B = [b1 0 0 b2 0 0 b3 0 0 b4 0 0;
    //      0 b1 0 0 b2 0 0 b3 0 0 b4 0;
    //      0 0 b1 0 0 b2 0 0 b3 0 0 b4;
    //      b1 b1 0 b2 b2 0 b3 b3 0 b4 b4 0;
    //      0 b1 b1 0 b2 b2 0 b3 b3 0 b4 b4;
    //      b1 0 b1 b2 0 b2 b3 0 b3 b4 0 b4]
    // where bi = [bi_x, bi_y, bi_z] and b_i = Hinv.row(i).tail(3)
    
    // Get the b coefficients from Hinv (rows 0-3, columns 1-3)
    Eigen::Matrix<double, 4, 3> b;
    for (int i = 0; i < 4; ++i) {
        b.row(i) = Hinv.row(i).segment(1, 3);
    }
    
    // Build 6x12 B matrix
    Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        int base = i * 3;
        // eps_xx = du/dx = b1_x*u1 + b2_x*u2 + b3_x*u3 + b4_x*u4
        B(0, base + 0) = b(i, 0); // du/dx contribution from ui
        B(1, base + 1) = b(i, 1); // dv/dy contribution from vi
        B(2, base + 2) = b(i, 2); // dw/dz contribution from wi
        // eps_xy = du/dy + dv/dx
        B(3, base + 0) = b(i, 1); // du/dy
        B(3, base + 1) = b(i, 0); // dv/dx
        // eps_yz = dv/dz + dw/dy
        B(4, base + 1) = b(i, 2); // dv/dz
        B(4, base + 2) = b(i, 1); // dw/dy
        // eps_xz = du/dz + dw/dx
        B(5, base + 0) = b(i, 2); // du/dz
        B(5, base + 2) = b(i, 0); // dw/dx
    }
    
    // Stiffness matrix K = B^T * G * B * detH/6.0  <-- FIXED: use detH instead of abs(detH)
    double volumeFactor = detH / 6.0;
    Eigen::Matrix<double, 12, 12> K = B.transpose() * G * B * volumeFactor;
    
    return K;
}

// Assemble global stiffness matrix for wedge element (6 nodes, 18 DOF)
Eigen::Matrix<double, 18, 18> computeWedgeStiffness(
    const std::vector<Eigen::Vector3d>& nodes,
    double E, double nu) {
    
    // Wedge is decomposed into 3 tetrahedrons according to the mapping in Fortran
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
    
    // But the Fortran says for ITYPE=1 (wedge), NTET=12 and ITET=1, so all 12 tets are used
    // However, the test case only requires the standard wedge decomposition.
    // Looking at the comment: "IOPT = 1 IMPLIES WEDGE - 3 TETRAHEDRONS"
    // So we'll use the first 3 tetrahedrons as the primary decomposition.
    // But the Fortran code actually uses all 12 for wedge? Let's check the logic.
    // In KSOLID, for ITYPE=1: ITET=1, NTET=12, so it loops from 1 to 12.
    // However, the problem states "Wedge/pentahedral solid element, 6 nodes", and the standard approach
    // is to split into 3 tets. But the Fortran has 12 mappings, so we'll use all 12.
    
    // Define the 12 tetrahedron node mappings (1-indexed in Fortran, convert to 0-indexed)
    std::vector<std::vector<int>> tetMappings = {
        {0,1,2,3}, // M1: 1,2,3,4
        {0,1,2,4}, // M2: 1,2,3,5
        {0,1,2,5}, // M3: 1,2,3,6
        {0,3,4,5}, // M4: 1,4,5,6
        {1,3,4,5}, // M5: 2,4,5,6
        {2,3,4,5}, // M6: 3,4,5,6
        {1,0,3,5}, // M7: 2,1,4,6
        {1,2,3,5}, // M8: 2,3,4,6
        {0,2,3,4}, // M9: 1,3,4,5
        {1,2,3,4}, // M10: 2,3,4,5
        {2,0,4,5}, // M11: 3,1,5,6
        {1,0,4,5}  // M12: 2,1,5,6
    };
    
    // Initialize global stiffness matrix (18x18)
    Eigen::Matrix<double, 18, 18> K_global = Eigen::Matrix<double, 18, 18>::Zero();
    
    // Process each tetrahedron
    for (int t = 0; t < 12; ++t) {
        // Get the 4 node coordinates for this tet
        std::vector<Eigen::Vector3d> tetCoords(4);
        for (int i = 0; i < 4; ++i) {
            tetCoords[i] = nodes[tetMappings[t][i]];
        }
        
        // Compute tetrahedron stiffness (12x12)
        Eigen::Matrix<double, 12, 12> K_tet = computeTetStiffness(tetCoords, E, nu);
        
        // Map local DOFs to global DOFs
        // Each node has 3 DOFs: u,v,w
        // Local DOF ordering: u1,v1,w1,u2,v2,w2,u3,v3,w3,u4,v4,w4
        // Global DOF ordering: u1,v1,w1,u2,v2,w2,u3,v3,w3,u4,v4,w4,u5,v5,w5,u6,v6,w6
        std::vector<int> localToGlobal(12);
        for (int i = 0; i < 4; ++i) {
            int nodeIdx = tetMappings[t][i];
            localToGlobal[3*i + 0] = 3 * nodeIdx + 0; // u
            localToGlobal[3*i + 1] = 3 * nodeIdx + 1; // v
            localToGlobal[3*i + 2] = 3 * nodeIdx + 2; // w
        }
        
        // Add contribution to global matrix
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                int gi = localToGlobal[i];
                int gj = localToGlobal[j];
                K_global(gi, gj) += K_tet(i, j);
            }
        }
    }
    
    return K_global;
}

// Main function
int main() {
    // Test case: Wedge (6 nodes)
    // N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    std::vector<Eigen::Vector3d> nodes = {
        Eigen::Vector3d(0.0, 0.0, 0.0), // N1
        Eigen::Vector3d(1.0, 0.0, 0.0), // N2
        Eigen::Vector3d(0.0, 1.0, 0.0), // N3
        Eigen::Vector3d(0.0, 0.0, 1.0), // N4
        Eigen::Vector3d(1.0, 0.0, 1.0), // N5
        Eigen::Vector3d(0.0, 1.0, 1.0)  // N6
    };
    
    double E = 200e9;      // 200 GPa
    double nu = 0.3;
    
    try {
        // Compute wedge stiffness matrix
        Eigen::Matrix<double, 18, 18> K = computeWedgeStiffness(nodes, E, nu);
        
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