#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Constants
constexpr double DEGRA = M_PI / 180.0; // degrees to radians

// Helper function to compute cross product of two 3D vectors
Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return Eigen::Vector3d(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
}

// Helper function to compute dot product
double dot(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return a.dot(b);
}

// Helper function to normalize vector
Eigen::Vector3d normalize(const Eigen::Vector3d& v) {
    double norm = v.norm();
    if (norm == 0.0) {
        return Eigen::Vector3d::Zero();
    }
    return v / norm;
}

// Material properties structure
struct MaterialProps {
    double E;
    double nu;
    double t; // thickness
    double I; // moment of inertia = t^3/12
};

// Transformation matrix for coordinate system rotation
Eigen::Matrix3d getTransformationMatrix(double theta) {
    double cos_theta = std::cos(theta);
    double sin_theta = std::sin(theta);
    Eigen::Matrix3d T;
    T << cos_theta, -sin_theta, 0.0,
         sin_theta,  cos_theta, 0.0,
               0.0,        0.0, 1.0;
    return T;
}

// Compute the basic bending triangle stiffness contribution
// Returns a 9x9 matrix as 9 Eigen::Matrix3d blocks (row-major: KU[0] = (1,1), KU[1] = (1,2), ..., KU[8] = (3,3))
std::vector<Eigen::Matrix3d> computeTriangleStiffness(
    const Eigen::Vector3d& A, 
    const Eigen::Vector3d& B, 
    const Eigen::Vector3d& C,
    const MaterialProps& mat) {
    
    // Compute vectors
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;
    
    // Compute I-vector (AB normalized)
    double XSUBB = AB.norm();
    if (XSUBB < 1e-12) {
        throw std::runtime_error("Degenerate triangle: zero length AB");
    }
    Eigen::Vector3d I_vec = AB / XSUBB;
    
    // Compute K-vector (cross product of AB and AC, normalized)
    Eigen::Vector3d K_vec = cross(AB, AC);
    double YSUBC = K_vec.norm();
    if (YSUBC < 1e-12) {
        throw std::runtime_error("Degenerate triangle: coplanar points");
    }
    K_vec = K_vec / YSUBC;
    
    // Compute J-vector = K × I
    Eigen::Vector3d J_vec = cross(K_vec, I_vec);
    J_vec = normalize(J_vec);
    
    // Compute coordinates in local system
    double XSUBC = dot(I_vec, AC);
    double YSUBC_local = dot(K_vec, AC);
    
    // Area of triangle
    double AREA = XSUBB * YSUBC_local / 2.0;
    
    // Centroid coordinates
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC_local / 3.0;
    
    // Precompute terms
    double XCSQ = XSUBC * XSUBC;
    double YCSQ = YSUBC_local * YSUBC_local;
    double XBSQ = XSUBB * XSUBB;
    double XCYC = XSUBC * YSUBC_local;
    double PX2 = (XBSQ + XSUBB * XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC_local * (XSUBB + 2.0 * XSUBC) / 12.0;
    double XBAR3 = 3.0 * XBAR;
    double YBAR3 = 3.0 * YBAR;
    double YBAR2 = 2.0 * YBAR;
    
    // Material matrix D = E*I/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2]
    double D_factor = mat.E * mat.I / (1.0 - mat.nu * mat.nu);
    double D11 = D_factor;
    double D12 = D_factor * mat.nu;
    double D22 = D_factor;
    double D33 = D_factor * (1.0 - mat.nu) / 2.0;
    
    // Build K^X matrix (6x6) stored as 36 elements
    std::vector<double> KX(36, 0.0);
    
    KX[0] = D11;                    // (1,1)
    KX[1] = D12;                    // (1,2)
    KX[2] = D12;                    // (1,3)
    KX[3] = D11 * XBAR3;            // (1,4)
    KX[4] = D12 * XBAR + YBAR2 * D12; // (1,5) - note: D12 appears twice in original
    KX[5] = D12 * YBAR3;            // (1,6)
    
    KX[6] = D12;                    // (2,1)
    KX[7] = D33;                    // (2,2)
    KX[8] = D33;                    // (2,3)
    KX[9] = D12 * XBAR3;            // (2,4)
    KX[10] = D33 * XBAR + YBAR2 * D33; // (2,5)
    KX[11] = D33 * YBAR3;           // (2,6)
    
    KX[12] = D12;                   // (3,1)
    KX[13] = D33;                   // (3,2)
    KX[14] = D22;                   // (3,3)
    KX[15] = D12 * XBAR3;           // (3,4)
    KX[16] = D22 * XBAR + YBAR2 * D33; // (3,5)
    KX[17] = D22 * YBAR3;           // (3,6)
    
    KX[18] = D11 * XBAR3;           // (4,1)
    KX[19] = D12 * XBAR3;           // (4,2)
    KX[20] = D12 * XBAR3;           // (4,3)
    KX[21] = D11 * 9.0 * PX2;       // (4,4)
    KX[22] = D12 * 3.0 * PX2 + 6.0 * PXY2 * D12; // (4,5)
    KX[23] = D12 * 9.0 * PXY2;      // (4,6)
    
    KX[24] = D12 * XBAR + YBAR2 * D12; // (5,1)
    KX[25] = D33 * XBAR + YBAR2 * D33; // (5,2)
    KX[26] = D22 * XBAR + YBAR2 * D33; // (5,3)
    KX[27] = D12 * 3.0 * PX2 + 6.0 * PXY2 * D12; // (5,4)
    KX[28] = D22 * PX2 + 4.0 * PXY2 * D33 + 4.0 * PY2 * D33; // (5,5)
    KX[29] = D22 * 3.0 * PXY2 + 6.0 * PY2 * D33; // (5,6)
    
    KX[30] = D12 * YBAR3;           // (6,1)
    KX[31] = D33 * YBAR3;           // (6,2)
    KX[32] = D22 * YBAR3;           // (6,3)
    KX[33] = D12 * 9.0 * PXY2;      // (6,4)
    KX[34] = D22 * 3.0 * PXY2 + 6.0 * PY2 * D33; // (6,5)
    KX[35] = D22 * 9.0 * PY2;       // (6,6)
    
    // Scale by 4*AREA
    double scale = 4.0 * AREA;
    for (auto& val : KX) {
        val *= scale;
    }
    
    // Convert to 3x3 blocks for K^U (9x9 matrix): 9 blocks of 3x3
    // Layout: KU[i] corresponds to block (i/3 + 1, i%3 + 1) i.e., row-major 3x3 block grid
    std::vector<Eigen::Matrix3d> KU(9);
    
    // Block (1,1): rows 0-2, cols 0-2
    KU[0] << KX[0], KX[1], KX[2],
              KX[6], KX[7], KX[8],
              KX[12], KX[13], KX[14];
    
    // Block (1,2): rows 0-2, cols 3-5
    KU[1] << KX[3], KX[4], KX[5],
              KX[9], KX[10], KX[11],
              KX[15], KX[16], KX[17];
    
    // Block (1,3): rows 0-2, cols 6-8 → but KX only has 6 columns → zeros
    KU[2].setZero();
    
    // Block (2,1): rows 3-5, cols 0-2
    KU[3] << KX[18], KX[19], KX[20],
              KX[24], KX[25], KX[26],
              KX[30], KX[31], KX[32];
    
    // Block (2,2): rows 3-5, cols 3-5
    KU[4] << KX[21], KX[22], KX[23],
              KX[27], KX[28], KX[29],
              KX[33], KX[34], KX[35];
    
    // Block (2,3): rows 3-5, cols 6-8 → zeros
    KU[5].setZero();
    
    // Block (3,1): rows 6-8, cols 0-2 → zeros (KX has no rows 6-8)
    KU[6].setZero();
    
    // Block (3,2): rows 6-8, cols 3-5 → zeros
    KU[7].setZero();
    
    // Block (3,3): rows 6-8, cols 6-8 → zeros
    KU[8].setZero();
    
    return KU;
}

// Main KQDPLT routine
Eigen::Matrix<double, 12, 12> KQDPLT(
    const std::vector<Eigen::Vector3d>& nodes,
    const MaterialProps& mat) {
    
    // Validate input: 4 nodes
    if (nodes.size() != 4) {
        throw std::runtime_error("KQDPLT requires exactly 4 nodes");
    }
    
    // M-matrix mapping for sub-triangles: 12 elements (1-based indices)
    // Original Fortran: M = [2,4,1, 3,1,2, 4,2,3, 1,3,4]
    std::vector<int> M = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    
    // Initialize full 12x12 stiffness matrix
    Eigen::Matrix<double, 12, 12> K_total = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Process each of the 4 sub-triangles (j=0 to 3)
    for (int j = 0; j < 4; ++j) {
        int km = 3 * j; // base index in M for this triangle
        int SUBSCA = M[km];     // 1-based node index
        int SUBSCB = M[km+1];   // 1-based node index
        int SUBSCC = M[km+2];   // 1-based node index
        
        // Convert to 0-based indices for nodes vector
        int idxA = SUBSCA - 1;
        int idxB = SUBSCB - 1;
        int idxC = SUBSCC - 1;
        
        // Validate indices
        if (idxA < 0 || idxA >= 4 || idxB < 0 || idxB >= 4 || idxC < 0 || idxC >= 4) {
            throw std::runtime_error("Invalid node index in M array");
        }
        
        Eigen::Vector3d A = nodes[idxA];
        Eigen::Vector3d B = nodes[idxB];
        Eigen::Vector3d C = nodes[idxC];
        
        try {
            auto KU = computeTriangleStiffness(A, B, C, mat);
            
            // KU has 9 blocks: KU[0] = (1,1), KU[1] = (1,2), KU[2] = (1,3),
            //                 KU[3] = (2,1), KU[4] = (2,2), KU[5] = (2,3),
            //                 KU[6] = (3,1), KU[7] = (3,2), KU[8] = (3,3)
            // Map local node indices (0,1,2) to global node indices (idxA, idxB, idxC)
            std::vector<int> local_to_global = {idxA, idxB, idxC};
            
            for (int i = 0; i < 3; ++i) {
                for (int j_block = 0; j_block < 3; ++j_block) {
                    int global_i = local_to_global[i];
                    int global_j = local_to_global[j_block];
                    int block_idx = i * 3 + j_block; // KU[block_idx] is (i+1, j_block+1)
                    
                    // Each node has 3 DOFs: w, theta_x, theta_y -> start at global_i*3, global_j*3
                    int dof_i_start = global_i * 3;
                    int dof_j_start = global_j * 3;
                    
                    for (int r = 0; r < 3; ++r) {
                        for (int c = 0; c < 3; ++c) {
                            K_total(dof_i_start + r, dof_j_start + c) += KU[block_idx](r, c);
                        }
                    }
                }
            }
        } catch (...) {
            // Skip degenerate triangle
        }
    }
    
    // Enforce symmetry
    K_total = (K_total + K_total.transpose()) * 0.5;
    
    return K_total;
}

int main() {
    try {
        // Test case: Quad A=(0,0,0), B=(1,0,0), C=(1,1,0), D=(0,1,0)
        std::vector<Eigen::Vector3d> nodes = {
            Eigen::Vector3d(0.0, 0.0, 0.0), // A (node 1)
            Eigen::Vector3d(1.0, 0.0, 0.0), // B (node 2)
            Eigen::Vector3d(1.0, 1.0, 0.0), // C (node 3)
            Eigen::Vector3d(0.0, 1.0, 0.0)  // D (node 4)
        };
        
        // Material properties
        MaterialProps mat;
        mat.E = 200e9;      // Pa
        mat.nu = 0.3;
        mat.t = 0.01;       // m
        mat.I = mat.t * mat.t * mat.t / 12.0; // m^4
        
        // Compute stiffness matrix
        Eigen::Matrix<double, 12, 12> K = KQDPLT(nodes, mat);
        
        // Output as plain text in JSON-like array format (scientific notation)
        std::cout << "{\n  \"stiffness_matrix\": [\n";
        for (int i = 0; i < 12; ++i) {
            std::cout << "    [";
            for (int j = 0; j < 12; ++j) {
                std::cout << std::setprecision(12) << std::scientific << K(i, j);
                if (j < 11) std::cout << ", ";
            }
            std::cout << "]";
            if (i < 11) std::cout << ",";
            std::cout << "\n";
        }
        std::cout << "  ]\n}\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}