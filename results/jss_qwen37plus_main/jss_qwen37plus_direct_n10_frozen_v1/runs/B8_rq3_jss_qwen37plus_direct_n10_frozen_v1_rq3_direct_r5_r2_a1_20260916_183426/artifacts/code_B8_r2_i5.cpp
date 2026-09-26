#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Forward declarations
void ktetra(int iopt, int jtype, std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords);

// Helper function to compute cross product of two 3D vectors
Eigen::Vector3d sAXB(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return Eigen::Vector3d(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
}

// Helper function to compute dot product of two 3D vectors
double sADOTB(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.dot(b);
}

// Helper function to invert a 4x4 matrix and compute determinant
bool invert4x4(const Eigen::Matrix4d& H, Eigen::Matrix4d& H_inv, double& det) {
    det = H.determinant();
    if (std::abs(det) < 1e-15) {
        return false;
    }
    H_inv = H.inverse();
    return true;
}

// Helper function for plane test (simplified: check if four points are coplanar)
// Returns true if points are coplanar (within tolerance)
bool kPLTST(const std::vector<double>& p1, const std::vector<double>& p2,
            const std::vector<double>& p3, const std::vector<double>& p4) {
    Eigen::Vector3d v1(p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]);
    Eigen::Vector3d v2(p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2]);
    Eigen::Vector3d v3(p4[0]-p1[0], p4[1]-p1[1], p4[2]-p1[2]);
    
    double scalar_triple = std::abs(v1.dot(sAXB(v2, v3)));
    return scalar_triple < 1e-12;
}

// Main KSOLID subroutine equivalent
void ksolid(int itype, std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords) {
    // For wedge (itype == 1), we use the 12-tetrahedron decomposition
    // The mapping matrix M for wedge: 12 tetrahedrons, each with 4 nodes
    std::vector<std::vector<int>> M = {
        {0, 1, 2, 3}, // 1,2,3,4 (0-indexed)
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

    // Geometry validation
    // Base triangle: nodes 0,1,2 (0,0,0), (1,0,0), (0,1,0)
    // Top triangle: nodes 3,4,5 (0,0,1), (1,0,1), (0,1,1)
    Eigen::Vector3d r12 = Eigen::Vector3d(node_coords[1][0]-node_coords[0][0],
                                          node_coords[1][1]-node_coords[0][1],
                                          node_coords[1][2]-node_coords[0][2]);
    Eigen::Vector3d r13 = Eigen::Vector3d(node_coords[2][0]-node_coords[0][0],
                                          node_coords[2][1]-node_coords[0][1],
                                          node_coords[2][2]-node_coords[0][2]);
    Eigen::Vector3d rxr = sAXB(r12, r13);

    Eigen::Vector3d r12_top = Eigen::Vector3d(node_coords[4][0]-node_coords[3][0],
                                              node_coords[4][1]-node_coords[3][1],
                                              node_coords[4][2]-node_coords[3][2]);
    Eigen::Vector3d r13_top = Eigen::Vector3d(node_coords[5][0]-node_coords[3][0],
                                              node_coords[5][1]-node_coords[3][1],
                                              node_coords[5][2]-node_coords[3][2]);
    Eigen::Vector3d r = sAXB(r12_top, r13_top);

    // Check orientation - should be same sign for valid wedge
    if (sADOTB(r, rxr) <= 0.0) {
        // Bad geometry
        return;
    }

    // Plane tests
    if (!kPLTST(node_coords[0], node_coords[1], node_coords[4], node_coords[3])) return;
    if (!kPLTST(node_coords[0], node_coords[3], node_coords[5], node_coords[2])) return;
    if (!kPLTST(node_coords[1], node_coords[2], node_coords[5], node_coords[4])) return;

    // Process each tetrahedron
    for (int i = 0; i < 12; ++i) {
        // Build ECPT for this tetrahedron
        // ECPT layout for tetrahedron: [el_id, mat_id, g1, g2, g3, g4, csid1, x1, y1, z1, csid2, x2, y2, z2, csid3, x3, y3, z3, csid4, x4, y4, z4, temp]
        std::vector<double> tet_ecpt(100, 0.0);
        
        // Element ID and material ID (same for all)
        tet_ecpt[0] = ecpt[0]; // element ID
        tet_ecpt[1] = ecpt[1]; // material ID
        
        // Grid point IDs (1-indexed in NASTRAN, but we'll use 0-indexed internally)
        tet_ecpt[2] = static_cast<double>(M[i][0] + 1);
        tet_ecpt[3] = static_cast<double>(M[i][1] + 1);
        tet_ecpt[4] = static_cast<double>(M[i][2] + 1);
        tet_ecpt[5] = static_cast<double>(M[i][3] + 1);
        
        // Coordinates for each grid point
        int idx = 6;
        for (int j = 0; j < 4; ++j) {
            int node_idx = M[i][j];
            tet_ecpt[idx++] = 0.0; // csid (coordinate system ID) - assume 0
            tet_ecpt[idx++] = node_coords[node_idx][0];
            tet_ecpt[idx++] = node_coords[node_idx][1];
            tet_ecpt[idx++] = node_coords[node_idx][2];
        }
        
        // Temperature (assume 0)
        tet_ecpt[23] = 0.0;
        
        // IOPT: for wedge, first 6 configurations get *2 factor, so IOPT = i+11 for i=0..5, then i+10 for i=6..11
        int iopt = (i < 6) ? (i + 11) : (i + 10);
        
        // JTYPE = 1 for wedge
        ktetra(iopt, 1, tet_ecpt, stiffness_matrix, E, nu, node_coords);
    }
}

// KTETRA subroutine equivalent
void ktetra(int iopt, int jtype, std::vector<double>& ecpt, Eigen::MatrixXd& stiffness_matrix,
            double E, double nu, const std::vector<std::vector<double>>& node_coords) {
    // Extract coordinates for the 4-node tetrahedron
    std::vector<Eigen::Vector3d> coords(4);
    for (int i = 0; i < 4; ++i) {
        int node_idx = static_cast<int>(ecpt[2+i]) - 1; // convert to 0-indexed
        coords[i] = Eigen::Vector3d(node_coords[node_idx][0], 
                                   node_coords[node_idx][1], 
                                   node_coords[node_idx][2]);
    }
    
    // Build H matrix: 4x4 [1 x y z] for each node
    Eigen::Matrix4d H;
    H << 1.0, coords[0](0), coords[0](1), coords[0](2),
         1.0, coords[1](0), coords[1](1), coords[1](2),
         1.0, coords[2](0), coords[2](1), coords[2](2),
         1.0, coords[3](0), coords[3](1), coords[3](2);
    
    double hdet;
    Eigen::Matrix4d H_inv;
    if (!invert4x4(H, H_inv, hdet)) {
        return; // bad geometry
    }
    
    hdet = std::abs(hdet);
    
    // Material matrix G (6x6 stress-strain matrix for isotropic material)
    Eigen::Matrix6d G;
    G.setZero();
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(temp1) < 1e-12) {
        return; // invalid nu
    }
    
    double c1 = E * (1.0 - nu) / temp1;
    double c2 = E * nu / temp1;
    double gg = E / (2.0 * (1.0 + nu)); // shear modulus
    
    G(0,0) = c1; G(1,1) = c1; G(2,2) = c1;
    G(0,1) = c2; G(0,2) = c2; G(1,0) = c2; G(1,2) = c2; G(2,0) = c2; G(2,1) = c2;
    G(3,3) = gg; G(4,4) = gg; G(5,5) = gg;
    
    // C matrices: 4 matrices of size 6x3
    // Each C_i = [dN_i/dx, dN_i/dy, dN_i/dz] where N_i are shape functions
    // For tetrahedron: N_i = H_inv_row_i * [1 x y z]^T, so dN_i/dx = H_inv(i,1), etc.
    std::vector<Eigen::Matrix<double,6,3>> C(4);
    for (int i = 0; i < 4; ++i) {
        C[i].setZero();
        // dN_i/dx = H_inv(i,1), dN_i/dy = H_inv(i,2), dN_i/dz = H_inv(i,3)
        double dx = H_inv(i,1);
        double dy = H_inv(i,2);
        double dz = H_inv(i,3);
        
        // For solid element, the strain-displacement matrix has:
        // [dx 0 0; 0 dy 0; 0 0 dz; dy dx 0; 0 dz dy; dz 0 dx] for each node
        // But standard B-matrix for tetrahedron is:
        // Row 0: [dx, 0, 0]
        // Row 1: [0, dy, 0]
        // Row 2: [0, 0, dz]
        // Row 3: [dy, dx, 0]
        // Row 4: [0, dz, dy]
        // Row 5: [dz, 0, dx]
        C[i](0,0) = dx;
        C[i](1,1) = dy;
        C[i](2,2) = dz;
        C[i](3,0) = dy; C[i](3,1) = dx;
        C[i](4,1) = dz; C[i](4,2) = dy;
        C[i](5,0) = dz; C[i](5,2) = dx;
    }
    
    // Scale factor for volume integration
    double scale_factor = hdet / 6.0;
    
    // For wedge: first 6 tetrahedrons get factor 2, and all get divided by 6
    // In NASTRAN, the scaling is more complex, but for our test case we use:
    // For wedge tetrahedrons: scale_factor = |det(H)| / 6.0, then apply additional factors
    if (iopt >= 11 && iopt <= 16) {
        scale_factor *= 2.0;
    }
    
    // Compute element stiffness matrix: K = sum_{i,j} C_i^T * G * C_j * V
    // where V is the volume contribution (scale_factor)
    Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(12, 12); // 4 nodes * 3 DOF = 12
    
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            // Compute C_i^T * G * C_j (3x3 matrix for DOFs i and j)
            Eigen::Matrix3d kij = C[i].transpose() * G * C[j];
            
            // Place in global stiffness matrix
            for (int di = 0; di < 3; ++di) {
                for (int dj = 0; dj < 3; ++dj) {
                    int row = i * 3 + di;
                    int col = j * 3 + dj;
                    ke(row, col) += kij(di, dj) * scale_factor;
                }
            }
        }
    }
    
    // Map local DOFs to global DOFs (6 nodes, 3 DOF each = 18 DOFs)
    // The wedge has 6 nodes, but each tetrahedron uses only 4 of them
    // We need to map the 4 local nodes to their global positions
    std::vector<int> local_to_global(4);
    for (int i = 0; i < 4; ++i) {
        int node_idx = static_cast<int>(ecpt[2+i]) - 1; // 0-indexed node number
        local_to_global[i] = node_idx;
    }
    
    // Add contribution to global stiffness matrix
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            for (int di = 0; di < 3; ++di) {
                for (int dj = 0; dj < 3; ++dj) {
                    int global_row = local_to_global[i] * 3 + di;
                    int global_col = local_to_global[j] * 3 + dj;
                    stiffness_matrix(global_row, global_col) += ke(i*3+di, j*3+dj);
                }
            }
        }
    }
}

// Convert double to string with scientific notation
std::string double_to_scientific(double value) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << value;
    std::string result = oss.str();
    // Remove trailing zeros after decimal point
    size_t pos = result.find_last_not_of('0');
    if (pos != std::string::npos && result[pos] == '.') {
        result.erase(pos + 1);
    } else if (pos != std::string::npos) {
        result.erase(pos + 1);
    }
    return result;
}

int main() {
    // Test case: Wedge (6 nodes): N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1), N5=(1,0,1), N6=(0,1,1)
    std::vector<std::vector<double>> node_coords = {
        {0.0, 0.0, 0.0}, // N1
        {1.0, 0.0, 0.0}, // N2
        {0.0, 1.0, 0.0}, // N3
        {0.0, 0.0, 1.0}, // N4
        {1.0, 0.0, 1.0}, // N5
        {0.0, 1.0, 1.0}  // N6
    };
    
    // Material properties
    double E = 200e9;      // 200 GPa
    double nu = 0.3;
    
    // Initialize ECPT array (simplified - just need element ID and material ID)
    std::vector<double> ecpt(100, 0.0);
    ecpt[0] = 1.0; // element ID
    ecpt[1] = 1.0; // material ID
    
    // Initialize 18x18 stiffness matrix
    Eigen::MatrixXd stiffness_matrix = Eigen::MatrixXd::Zero(18, 18);
    
    // Call KSOLID for wedge (itype = 1)
    ksolid(1, ecpt, stiffness_matrix, E, nu, node_coords);
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 18; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << double_to_scientific(stiffness_matrix(i, j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}