#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Define M_PI if not already defined (for portability)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations for required subroutines (we'll implement only what's needed)
void ktrmem(int flag);
void sma1b(const Eigen::Matrix<double, 3, 3>& kij, int grid_id, int npvt, int ifkgg, double temp);

// Global data structures (mimicking Fortran COMMON blocks)
struct CommonData {
    // /CONDAS/
    std::vector<double> consts; // 5 elements
    
    // /SMA1HT/
    bool heat;
    
    // /SMA1ET/ - ECPT array (100 elements)
    std::vector<double> ecpt;
    
    // /SMA1IO/
    std::vector<double> dum1;
    int ifkgg;
    std::vector<double> dum2;
    int if4gg;
    std::vector<double> dum3;
    
    // /SMA1CL/
    int iopt4;
    int k4ggsw;
    int npvt;
    std::vector<double> dumcl;
    std::vector<int> link;
    int idetck;
    bool dodet;
    int nogO;
    
    // /SMA1DP/
    std::vector<double> kij;      // 6x6 stiffness matrix storage
    std::vector<double> dum7;
    std::vector<double> ksum;     // 6x6 sum storage
    double temp;
    double cosang;
    double sinang;
    double vecl;
    std::vector<double> ivec;
    std::vector<double> jvec;
    std::vector<double> kvec;
    std::vector<double> pvec;
    std::vector<double> vsubk;
    std::vector<double> v;
    std::vector<double> si;
    int npivot;
    int mpoint;
    int mi;
    int nsubsc;
    std::vector<int> ngrid;
    double u1;
    double u2;
    std::vector<double> coord;
    std::vector<double> dummy88; // renamed from dum88 to avoid conflict
    
    // /MATIN/
    int matid;
    int inflag;
    double eltemp;
    double stress;
    double sinth;
    double costh;
    
    // /MATOUT/
    std::vector<double> dum99;
    double gsube;
    std::vector<double> dum88_matout; // renamed from dum88 to avoid conflict
    
    // K3X3 array (27 elements) - stored in same memory as kij for first 9 elements
    // We'll use a separate vector for clarity
    std::vector<double> k3x3;
    
    // Constructor to initialize all members
    CommonData() : 
        consts(5, 0.0),
        heat(false),
        ecpt(100, 0.0),
        dum1(10, 0.0),
        ifkgg(0),
        dum2(1, 0.0),
        if4gg(0),
        dum3(23, 0.0),
        iopt4(0),
        k4ggsw(0),
        npvt(0),
        dumcl(7, 0.0),
        link(10, 0),
        idetck(0),
        dodet(false),
        nogO(0),
        kij(36, 0.0),
        dum7(156, 0.0),
        ksum(36, 0.0),
        temp(0.0),
        cosang(0.0),
        sinang(0.0),
        vecl(0.0),
        ivec(3, 0.0),
        jvec(3, 0.0),
        kvec(3, 0.0),
        pvec(3, 0.0),
        vsubk(3, 0.0),
        v(3, 0.0),
        si(3, 0.0),
        npivot(0),
        mpoint(0),
        mi(0),
        nsubsc(0),
        ngrid(4, 0),
        u1(0.0),
        u2(0.0),
        coord(16, 0.0),
        dummy88(248, 0.0),
        matid(0),
        inflag(0),
        eltemp(0.0),
        stress(0.0),
        sinth(0.0),
        costh(0.0),
        dum99(11, 0.0),
        gsube(0.0),
        dum88_matout(6, 0.0),
        k3x3(27, 0.0) {}
};

CommonData common;

// Helper function to get DEGRA (degrees to radians conversion factor)
double get_degra() {
    return M_PI / 180.0;
}

// Helper function to compute triangular membrane stiffness (simplified version)
// This implements the core physics of KTRMEM for isotropic material
void ktrmem(int flag) {
    // For CQUAD4, we need to compute the 3x3 stiffness matrix for a triangle
    // Using standard linear triangular membrane element formulation
    
    // Extract triangle coordinates from ECPT
    // ECPT(10-12): point A (x,y,z)
    // ECPT(14-16): point B (x,y,z)  
    // ECPT(18-20): point C (x,y,z)
    double x1 = common.ecpt[9];  // index 9 for ECPT(10)
    double y1 = common.ecpt[10]; // ECPT(11)
    double z1 = common.ecpt[11]; // ECPT(12)
    
    double x2 = common.ecpt[13]; // ECPT(14)
    double y2 = common.ecpt[14]; // ECPT(15)
    double z2 = common.ecpt[15]; // ECPT(16)
    
    double x3 = common.ecpt[17]; // ECPT(18)
    double y3 = common.ecpt[18]; // ECPT(19)
    double z3 = common.ecpt[19]; // ECPT(20)
    
    // Material properties from CONSTS
    double e = common.consts[0]; // Young's modulus
    double nu = common.consts[1]; // Poisson's ratio
    double t = common.ecpt[6];   // thickness (ECPT(7))
    
    // Compute area of triangle
    double area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    
    if (area == 0.0) {
        // Error handling
        common.nogO = 1;
        return;
    }
    
    // Compute strain-displacement matrix B (3x6) for plane stress membrane
    // Nodes: A, B, C → DOFs: (Ax,Ay, Bx,By, Cx,Cy)
    // b_i = y_j - y_k, c_i = x_k - x_j
    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;
    
    // B = [b1 0 b2 0 b3 0; 0 c1 0 c2 0 c3; c1 b1 c2 b2 c3 b3] / (2*area)
    double denom = 2.0 * area;
    Eigen::Matrix<double, 3, 6> B;
    B << b1/denom, 0.0,      b2/denom, 0.0,      b3/denom, 0.0,
         0.0,      c1/denom, 0.0,      c2/denom, 0.0,      c3/denom,
         c1/denom, b1/denom, c2/denom, b2/denom, c3/denom, b3/denom;
    
    // Material matrix D for plane stress
    double d11 = e / (1.0 - nu*nu);
    double d12 = nu * d11;
    double d22 = d11;
    double d33 = e / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3> D;
    D << d11, d12, 0.0,
         d12, d22, 0.0,
         0.0, 0.0, d33;
    
    // Stiffness matrix K = t * B^T * D * B * area
    Eigen::Matrix<double, 6, 6> K_tri = t * area * B.transpose() * D * B;
    
    // Now extract three 3x3 blocks for nodes A, B, C
    // Each node has 2 in-plane DOFs (x,y) and 1 out-of-plane (z) — but z is zero for membrane.
    // NASTRAN K3X3 stores for each node a 3x3 block: [Kxx, Kyx, Kzx, Kxy, Kyy, Kzy, Kxz, Kyz, Kzz]
    // So for node A (DOFs 0,1 → x,y), we take submatrix K_tri(0:1,0:1) and pad with zeros for z.
    // We'll build three 3x3 blocks: blockA, blockB, blockC
    std::vector<Eigen::Matrix<double, 3, 3>> blocks(3);
    
    // Block for node A (DOFs 0,1 → x,y)
    blocks[0].setZero();
    blocks[0](0,0) = K_tri(0,0); // Ax-Ax
    blocks[0](0,1) = K_tri(0,1); // Ax-Ay
    blocks[0](1,0) = K_tri(1,0); // Ay-Ax
    blocks[0](1,1) = K_tri(1,1); // Ay-Ay
    // z-rows/cols remain zero
    
    // Block for node B (DOFs 2,3 → x,y)
    blocks[1].setZero();
    blocks[1](0,0) = K_tri(2,2); // Bx-Bx
    blocks[1](0,1) = K_tri(2,3); // Bx-By
    blocks[1](1,0) = K_tri(3,2); // By-Bx
    blocks[1](1,1) = K_tri(3,3); // By-By
    
    // Block for node C (DOFs 4,5 → x,y)
    blocks[2].setZero();
    blocks[2](0,0) = K_tri(4,4); // Cx-Cx
    blocks[2](0,1) = K_tri(4,5); // Cx-Cy
    blocks[2](1,0) = K_tri(5,4); // Cy-Cx
    blocks[2](1,1) = K_tri(5,5); // Cy-Cy
    
    // Build triangle-local to element-local rotation matrix T (3x3)
    // Triangle local x-axis = vector AB, y-axis = in-plane normal to AB in triangle plane
    Eigen::Vector3d ab(x2-x1, y2-y1, z2-z1);
    Eigen::Vector3d ac(x3-x1, y3-y1, z3-z1);
    Eigen::Vector3d tri_normal = ab.cross(ac);
    tri_normal.normalize();
    
    // Triangle x-axis = AB normalized
    Eigen::Vector3d tri_x = ab.normalized();
    // Triangle y-axis = tri_normal × tri_x (in-plane, orthogonal to x)
    Eigen::Vector3d tri_y = tri_normal.cross(tri_x);
    tri_y.normalize();
    // Triangle z-axis = tri_normal
    
    // Element local axes are stored in common.ivec, common.jvec, common.kvec
    // They are orthonormal by construction
    Eigen::Matrix<double, 3, 3> T;
    T.col(0) = Eigen::Vector3d(common.ivec[0], common.ivec[1], common.ivec[2]);
    T.col(1) = Eigen::Vector3d(common.jvec[0], common.jvec[1], common.jvec[2]);
    T.col(2) = Eigen::Vector3d(common.kvec[0], common.kvec[1], common.kvec[2]);
    
    // Rotation from triangle-local to element-local: T_tri_to_elem = T^T * T_tri
    // But T_tri = [tri_x, tri_y, tri_z], so T_tri_to_elem = T.transpose() * [tri_x, tri_y, tri_z]
    Eigen::Matrix<double, 3, 3> T_tri;
    T_tri.col(0) = tri_x;
    T_tri.col(1) = tri_y;
    T_tri.col(2) = tri_normal;
    
    Eigen::Matrix<double, 3, 3> T_rot = T.transpose() * T_tri; // 3x3 rotation
    
    // Apply T_rot^T * block * T_rot to each block
    for (int idx = 0; idx < 3; ++idx) {
        Eigen::Matrix<double, 3, 3> rotated = T_rot.transpose() * blocks[idx] * T_rot;
        
        // Store in K3X3 in COLUMN-MAJOR order (Fortran convention)
        // K3X3(1:9) = block0, K3X3(10:18) = block1, K3X3(19:27) = block2
        int base = idx * 9;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                common.k3x3[base + i + 3*j] = rotated(i,j); // column-major: i + 3*j
            }
        }
    }
}

// Simplified SMA1B that accumulates into ksum
// In real NASTRAN, this would insert into global stiffness matrix
// Here, we accumulate into a global 8x8 matrix passed by reference
void sma1b(const Eigen::Matrix<double, 3, 3>& kij, int grid_id, int npvt, int ifkgg, double temp, 
           Eigen::Matrix<double, 8, 8>& K_global, const std::vector<int>& node_to_dof) {
    // grid_id is 1-based node ID (1,2,3,4)
    // Map to 0-based index in node_to_dof
    int node_idx = -1;
    for (int i = 0; i < 4; ++i) {
        if (static_cast<int>(common.ngrid[i]) == grid_id) {
            node_idx = i;
            break;
        }
    }
    if (node_idx == -1) return;
    
    // Get global DOF start index for this node: 2*node_idx (0,2,4,6)
    int dof_start = node_to_dof[node_idx];
    
    // Add kij(0,0), kij(0,1), kij(1,0), kij(1,1) to K_global[dof_start + i, dof_start + j]
    K_global(dof_start + 0, dof_start + 0) += kij(0,0);
    K_global(dof_start + 0, dof_start + 1) += kij(0,1);
    K_global(dof_start + 1, dof_start + 0) += kij(1,0);
    K_global(dof_start + 1, dof_start + 1) += kij(1,1);
}

// Main CQUAD4 stiffness computation function
Eigen::Matrix<double, 8, 8> compute_cquad4_stiffness(
    const std::vector<std::vector<double>>& nodes, // 4 nodes, each with [x,y]
    double E, double nu, double t) {
    
    // Initialize common data
    common.consts[0] = E;
    common.consts[1] = nu;
    common.consts[3] = get_degra(); // DEGRA
    
    // Set up ECPT array with test case data
    // ECPT(1) = element ID (we'll use 1)
    common.ecpt[0] = 1.0;
    
    // ECPT(2-5) = grid point IDs (1,2,3,4)
    common.ecpt[1] = 1.0; // node 1
    common.ecpt[2] = 2.0; // node 2
    common.ecpt[3] = 3.0; // node 3
    common.ecpt[4] = 4.0; // node 4
    
    // ECPT(6) = theta (0.0 for test case)
    common.ecpt[5] = 0.0;
    
    // ECPT(7) = material ID (1)
    common.ecpt[6] = 1.0;
    
    // ECPT(8) = thickness (t)
    common.ecpt[7] = t;
    
    // ECPT(9) = non-structural mass (0.0)
    common.ecpt[8] = 0.0;
    
    // ECPT(10-12) = coordinates of node 1
    common.ecpt[9]  = nodes[0][0]; // x1
    common.ecpt[10] = nodes[0][1]; // y1
    common.ecpt[11] = 0.0;        // z1
    
    // ECPT(13-15) = coordinates of node 2
    common.ecpt[12] = 0.0;        // CSID for node 1 (we'll ignore)
    common.ecpt[13] = nodes[1][0]; // x2
    common.ecpt[14] = nodes[1][1]; // y2
    common.ecpt[15] = 0.0;        // z2
    
    // ECPT(16-18) = coordinates of node 3
    common.ecpt[16] = 0.0;        // CSID for node 2
    common.ecpt[17] = nodes[2][0]; // x3
    common.ecpt[18] = nodes[2][1]; // y3
    common.ecpt[19] = 0.0;        // z3
    
    // ECPT(19-21) = coordinates of node 4
    common.ecpt[20] = 0.0;        // CSID for node 3
    common.ecpt[21] = nodes[3][0]; // x4
    common.ecpt[22] = nodes[3][1]; // y4
    common.ecpt[23] = 0.0;        // z4
    
    // ECPT(24-26) = more CSID info (set to 0)
    common.ecpt[24] = 0.0;
    common.ecpt[25] = 0.0;
    common.ecpt[26] = 0.0;
    
    // Set other parameters
    common.npvt = 1; // pivot point (node 1)
    common.heat = false;
    common.iopt4 = 0;
    
    // Save grid point numbers
    common.ngrid[0] = static_cast<int>(common.ecpt[1]); // node 1
    common.ngrid[1] = static_cast<int>(common.ecpt[2]); // node 2
    common.ngrid[2] = static_cast<int>(common.ecpt[3]); // node 3
    common.ngrid[3] = static_cast<int>(common.ecpt[4]); // node 4
    
    // Copy coordinates to COORD array
    for (int i = 0; i < 4; ++i) {
        int base_idx = i * 4;
        common.coord[base_idx]     = 0.0; // CSID
        common.coord[base_idx + 1] = nodes[i][0]; // x
        common.coord[base_idx + 2] = nodes[i][1]; // y
        common.coord[base_idx + 3] = 0.0; // z
    }
    
    // Set up coordinate system vectors (simplified for 2D)
    // IVEC = vector from node 1 to node 2
    common.ivec[0] = nodes[1][0] - nodes[0][0];
    common.ivec[1] = nodes[1][1] - nodes[0][1];
    common.ivec[2] = 0.0;
    
    double vecl = std::sqrt(common.ivec[0]*common.ivec[0] + 
                           common.ivec[1]*common.ivec[1] + 
                           common.ivec[2]*common.ivec[2]);
    if (vecl != 0.0) {
        common.ivec[0] /= vecl;
        common.ivec[1] /= vecl;
        common.ivec[2] /= vecl;
    }
    
    // KVEC = normal vector (cross product of IVEC and vector from node 1 to node 3)
    double v3x = nodes[2][0] - nodes[0][0];
    double v3y = nodes[2][1] - nodes[0][1];
    double v3z = 0.0;
    
    common.kvec[0] = common.ivec[1]*v3z - common.ivec[2]*v3y;
    common.kvec[1] = common.ivec[2]*v3x - common.ivec[0]*v3z;
    common.kvec[2] = common.ivec[0]*v3y - common.ivec[1]*v3x;
    
    vecl = std::sqrt(common.kvec[0]*common.kvec[0] + 
                    common.kvec[1]*common.kvec[1] + 
                    common.kvec[2]*common.kvec[2]);
    if (vecl != 0.0) {
        common.kvec[0] /= vecl;
        common.kvec[1] /= vecl;
        common.kvec[2] /= vecl;
    }
    
    // JVEC = cross product of KVEC and IVEC
    common.jvec[0] = common.kvec[1]*common.ivec[2] - common.kvec[2]*common.ivec[1];
    common.jvec[1] = common.kvec[2]*common.ivec[0] - common.kvec[0]*common.ivec[2];
    common.jvec[2] = common.kvec[0]*common.ivec[1] - common.kvec[1]*common.ivec[0];
    
    // PVEC = COSANG*IVEC + SINANG*JVEC (for theta=0, PVEC=IVEC)
    common.cosang = 1.0;
    common.sinang = 0.0;
    common.pvec[0] = common.cosang * common.ivec[0] + common.sinang * common.jvec[0];
    common.pvec[1] = common.cosang * common.ivec[1] + common.sinang * common.jvec[1];
    common.pvec[2] = common.cosang * common.ivec[2] + common.sinang * common.jvec[2];
    
    // Mapping matrix M (12 elements) — not used in fixed version
    // Find pivot point
    int npivot = 0;
    for (int i = 0; i < 4; ++i) {
        if (common.ngrid[i] == common.npvt) {
            npivot = i + 1; // 1-based indexing
            break;
        }
    }
    
    // Zero out KSUM (not used, but keep for compatibility)
    for (int i = 0; i < 36; ++i) {
        common.ksum[i] = 0.0;
    }
    
    // Build global 8x8 stiffness matrix
    Eigen::Matrix<double, 8, 8> K_global = Eigen::Matrix<double, 8, 8>::Zero();
    
    // Node-to-DOF mapping: node i → DOFs (2*i, 2*i+1) → indices 0,1,2,3,4,5,6,7
    std::vector<int> node_to_dof = {0, 2, 4, 6}; // node 0→DOFs 0,1; node 1→2,3; node 2→4,5; node 3→6,7
    
    // Process 4 triangles (standard CQUAD4 decomposition)
    // Triangles: (1,2,4), (2,3,1), (3,4,2), (4,1,3) — using 0-based node indices
    std::vector<std::vector<int>> triangles = {
        {0, 1, 3}, // triangle 1: nodes 1,2,4 → indices 0,1,3
        {1, 2, 0}, // triangle 2: nodes 2,3,1 → indices 1,2,0
        {2, 3, 1}, // triangle 3: nodes 3,4,2 → indices 2,3,1
        {3, 0, 2}  // triangle 4: nodes 4,1,3 → indices 3,0,2
    };
    
    // For each triangle, set up ECPT and call ktrmem
    for (int j = 0; j < 4; ++j) {
        // Set up triangle coordinates in ECPT
        const auto& tri = triangles[j];
        
        // ECPT(10-12): point A (first node of triangle)
        common.ecpt[9]  = nodes[tri[0]][0];
        common.ecpt[10] = nodes[tri[0]][1];
        common.ecpt[11] = 0.0;
        
        // ECPT(14-16): point B (second node of triangle)
        common.ecpt[13] = nodes[tri[1]][0];
        common.ecpt[14] = nodes[tri[1]][1];
        common.ecpt[15] = 0.0;
        
        // ECPT(18-20): point C (third node of triangle)
        common.ecpt[17] = nodes[tri[2]][0];
        common.ecpt[18] = nodes[tri[2]][1];
        common.ecpt[19] = 0.0;
        
        // Call ktrmem to compute and rotate K3X3
        ktrmem(1);
        
        // Now extract and assemble the three 3x3 blocks
        // Blocks are stored in K3X3: indices 0-8, 9-17, 18-26
        for (int blk_idx = 0; blk_idx < 3; ++blk_idx) {
            int base = blk_idx * 9;
            Eigen::Matrix<double, 3, 3> block;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    block(i,j) = common.k3x3[base + i + 3*j]; // column-major read
                }
            }
            
            // Map triangle node index to global node ID
            int tri_node_global_id = static_cast<int>(common.ecpt[1 + blk_idx]);
            // But ECPT(2-5) holds node IDs: ECPT(2)=node1, ECPT(3)=node2, ECPT(4)=node3, ECPT(5)=node4
            // So for block 0 → ECPT(2), block 1 → ECPT(3), block 2 → ECPT(4)
            int global_node_id;
            if (blk_idx == 0) global_node_id = static_cast<int>(common.ecpt[1]);
            else if (blk_idx == 1) global_node_id = static_cast<int>(common.ecpt[2]);
            else global_node_id = static_cast<int>(common.ecpt[3]);
            
            // Convert global node ID to 0-based node index (1→0, 2→1, 3→2, 4→3)
            int node_idx = -1;
            for (int i = 0; i < 4; ++i) {
                if (common.ngrid[i] == global_node_id) {
                    node_idx = i;
                    break;
                }
            }
            if (node_idx == -1) continue;
            
            // Accumulate into global matrix
            sma1b(block, global_node_id, common.npvt, common.ifkgg, common.temp, K_global, node_to_dof);
        }
    }
    
    return K_global;
}

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<double, 8, 8>& mat) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 8; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 8; ++j) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(15) << mat(i,j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    // Test case: Nodes (0,0), (2,0), (2,1.5), (0,1.5)
    std::vector<std::vector<double>> nodes = {
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.5},
        {0.0, 1.5}
    };
    
    double E = 200e9;   // Pa
    double nu = 0.3;
    double t = 0.01;    // m
    
    // Compute stiffness matrix
    Eigen::Matrix<double, 8, 8> K = compute_cquad4_stiffness(nodes, E, nu, t);
    
    // Output as JSON
    std::cout << matrix_to_json(K) << std::endl;
    
    return 0;
}