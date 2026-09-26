#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Using double precision throughout
using Real = double;

// Mock of KTRMEM subroutine - computes 3x3 stiffness for a triangle
// Based on standard linear triangular membrane element (CST)
void ktrmem(const std::vector<Real>& coords, Real E, Real nu, Real t, 
            Real sinth, Real costh, Eigen::Matrix<Real, 3, 3>& k3x3) {
    // coords: [x1,y1, x2,y2, x3,y3] - 6 values
    Real x1 = coords[0], y1 = coords[1];
    Real x2 = coords[2], y2 = coords[3];
    Real x3 = coords[4], y3 = coords[5];

    // Compute area of triangle
    Real area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    if (area == 0.0) {
        k3x3.setZero();
        return;
    }

    // Material matrix D for plane stress (isotropic)
    // D = (E/(1-nu^2)) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
    Real denom = 1.0 - nu*nu;
    Real d11 = E / denom;
    Real d12 = E * nu / denom;
    Real d33 = E * (1.0 - nu) / (2.0 * denom);

    // Strain-displacement matrix B
    // For CST: B = [[b1,0,b2,0,b3,0], [0,c1,0,c2,0,c3], [c1,b1,c2,b2,c3,b3]] / (2*area)
    Real b1 = y2 - y3;
    Real b2 = y3 - y1;
    Real b3 = y1 - y2;
    Real c1 = x3 - x2;
    Real c2 = x1 - x3;
    Real c3 = x2 - x1;

    // Rotation matrix for material orientation
    Eigen::Matrix<Real, 3, 3> T;
    T << costh*costh, sinth*sinth, 2.0*costh*sinth,
         sinth*sinth, costh*costh, -2.0*costh*sinth,
         -costh*sinth, costh*sinth, costh*costh - sinth*sinth;

    // Material matrix in rotated coordinates: D_rot = T * D * T^T
    Eigen::Matrix<Real, 3, 3> D;
    D << d11, d12, 0.0,
         d12, d11, 0.0,
         0.0, 0.0, d33;

    Eigen::Matrix<Real, 3, 3> D_rot = T * D * T.transpose();

    // B matrix (3x6)
    Eigen::Matrix<Real, 3, 6> B;
    B.setZero();
    Real inv_2A = 1.0 / (2.0 * area);
    B(0,0) = b1 * inv_2A; B(0,2) = b2 * inv_2A; B(0,4) = b3 * inv_2A;
    B(1,1) = c1 * inv_2A; B(1,3) = c2 * inv_2A; B(1,5) = c3 * inv_2A;
    B(2,0) = c1 * inv_2A; B(2,1) = b1 * inv_2A;
    B(2,2) = c2 * inv_2A; B(2,3) = b2 * inv_2A;
    B(2,4) = c3 * inv_2A; B(2,5) = b3 * inv_2A;

    // Element stiffness: k = t * area * B^T * D_rot * B
    Eigen::Matrix<Real, 6, 6> k6x6 = t * area * B.transpose() * D_rot * B;

    // Extract the 3x3 submatrix for the pivot node (first 2 DOFs of each node -> 3x3 for 3 nodes? But we need 3x3 per triangle for assembly)
    // Actually, for CQUAD4 decomposition, we need the 3x3 stiffness contribution for the three nodes of the triangle
    // But the Fortran code assembles into KSUM as 4 blocks of 3x3, then maps to 8x8
    // So we return the full 6x6 but the calling code will extract appropriate parts
    // However, the Fortran uses K3X3(27) to store 3 separate 3x3 matrices
    // Looking at the mapping, it appears they compute 3x3 contributions for each triangle's three nodes
    
    // Instead, we'll compute the standard CST 6x6 and let the assembly handle mapping
    // But the Fortran expects K3X3 to contain three 3x3 matrices (27 elements total)
    // From the assembly logic: for triangle J with nodes i,j,k, they add K3X3(1..9) to KSUM for node i,
    // K3X3(10..18) for node j, K3X3(19..27) for node k
    
    // So we need to compute the 3x3 contribution for each node pair
    // Standard approach: k3x3 for node i is k6x6.block(0,0,2,2), for node j is k6x6.block(2,2,2,2), for node k is k6x6.block(4,4,2,2)
    // But the Fortran mapping suggests different layout
    
    // Actually, looking at the Fortran assembly:
    // They have K3X3(1..27) and loop I=1,3 and K=1,9 to add K3X3(NPT1+K) to KSUM(NPT3)
    // Where NPT1 = 9*I-9, so for I=1: indices 1..9, I=2: 10..18, I=3: 19..27
    // And they map to KSUM positions based on M(MI) which gives node index
    
    // So K3X3 contains three 3x3 matrices: one for each node in the triangle
    // Each 3x3 is the contribution to that node's 2 DOFs (so actually 2x2, but stored in 3x3 layout with zeros?)
    
    // Given complexity, we'll compute the full 6x6 and extract the diagonal 2x2 blocks
    // Since each node has 2 DOFs, the 3x3 in Fortran is likely a misnomer and should be 2x2 stored in 3x3 array
    // But the data says "3 3X3 ARRAYS", so perhaps it's 3 matrices of size 3x3 for some reason
    
    // Simpler approach: implement the standard analytical formula for CST stiffness
    // and return the three 2x2 diagonal blocks as 3x3 matrices (with last row/col zero)
    
    // Let's compute the three 2x2 diagonal blocks
    Eigen::Matrix<Real, 2, 2> k11 = k6x6.block(0,0,2,2);
    Eigen::Matrix<Real, 2, 2> k22 = k6x6.block(2,2,2,2);
    Eigen::Matrix<Real, 2, 2> k33 = k6x6.block(4,4,2,2);
    
    // Store in k3x3 as three 3x3 matrices (padded with zeros)
    k3x3.setZero();
    k3x3.block(0,0,2,2) = k11;
    k3x3.block(3,3,2,2) = k22;
    k3x3.block(6,6,2,2) = k33;
}

// Mock SMA1B subroutine - accumulates into global stiffness
void sma1b(Eigen::Matrix<Real, 8, 8>& K_global, const Eigen::Matrix<Real, 3, 3>& k_local, 
           int node_idx, int pivot_node) {
    // Map local node index (1-based) to global DOF indices
    // Each node has 2 DOFs: ux=2*i-1, uy=2*i
    int dof_start = 2 * (node_idx - 1);
    
    // Add k_local(1:2,1:2) to K_global(dof_start+0:dof_start+1, dof_start+0:dof_start+1)
    // But k_local is 3x3, we use top-left 2x2
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            K_global(dof_start + i, dof_start + j) += k_local(i, j);
        }
    }
}

// Main CQUAD4 stiffness computation
Eigen::Matrix<Real, 8, 8> compute_cquad4_stiffness(
    const std::vector<std::vector<Real>>& nodes, // 4 nodes, each [x,y]
    Real E, Real nu, Real t) {
    
    // Validate input
    if (nodes.size() != 4) {
        throw std::runtime_error("CQUAD4 requires exactly 4 nodes");
    }
    
    // Extract node coordinates
    Real x1 = nodes[0][0], y1 = nodes[0][1];
    Real x2 = nodes[1][0], y2 = nodes[1][1];
    Real x3 = nodes[2][0], y3 = nodes[2][1];
    Real x4 = nodes[3][0], y4 = nodes[3][1];
    
    // Global stiffness matrix (8x8)
    Eigen::Matrix<Real, 8, 8> K_global = Eigen::Matrix<Real, 8, 8>::Zero();
    
    // Define the four triangles (as in NASTRAN: divide quad into 4 triangles from pivot)
    // Triangle mapping: M = [1,2,4, 2,3,1, 3,4,2, 4,1,3]
    // So triangles are: (1,2,4), (2,3,1), (3,4,2), (4,1,3)
    std::vector<std::vector<int>> triangles = {
        {1, 2, 4}, // triangle 1
        {2, 3, 1}, // triangle 2  
        {3, 4, 2}, // triangle 3
        {4, 1, 3}  // triangle 4
    };
    
    // Precompute coordinate system vectors for material orientation
    // IVEC = vector from node1 to node2
    Eigen::Vector3d ivec(x2 - x1, y2 - y1, 0.0);
    Real vecl = ivec.norm();
    if (vecl == 0.0) {
        throw std::runtime_error("Zero length IVEC");
    }
    ivec /= vecl;
    
    // VSUBK = cross product of IVEC and vector from node1 to node4
    Eigen::Vector3d v41(x4 - x1, y4 - y1, 0.0);
    Eigen::Vector3d vsubk = ivec.cross(v41);
    vecl = vsubk.norm();
    if (vecl == 0.0) {
        throw std::runtime_error("Zero length VSUBK");
    }
    Eigen::Vector3d kvec = vsubk / vecl;
    
    // JVEC = cross product of KVEC and IVEC
    Eigen::Vector3d jvec = kvec.cross(ivec);
    
    // Material orientation angle (theta = 0 for isotropic, so sin=0, cos=1)
    Real sinth = 0.0;
    Real costh = 1.0;
    
    // Process each triangle
    for (int tri_idx = 0; tri_idx < 4; ++tri_idx) {
        const auto& tri = triangles[tri_idx];
        
        // Get coordinates for this triangle
        std::vector<Real> coords(6);
        coords[0] = nodes[tri[0]-1][0]; coords[1] = nodes[tri[0]-1][1]; // node 1
        coords[2] = nodes[tri[1]-1][0]; coords[3] = nodes[tri[1]-1][1]; // node 2  
        coords[4] = nodes[tri[2]-1][0]; coords[5] = nodes[tri[2]-1][1]; // node 3
        
        // Compute triangle stiffness
        Eigen::Matrix<Real, 3, 3> k3x3;
        ktrmem(coords, E, nu, t, sinth, costh, k3x3);
        
        // Assemble into global matrix
        // For each node in triangle, add its 2x2 contribution
        for (int node_in_tri = 0; node_in_tri < 3; ++node_in_tri) {
            int global_node = tri[node_in_tri];
            // The k3x3 contains three 2x2 blocks at positions (0,0), (3,3), (6,6)
            // We need to add the appropriate block based on node_in_tri
            int block_offset = node_in_tri * 3;
            
            // Add 2x2 block to global stiffness
            int dof_start = 2 * (global_node - 1);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    K_global(dof_start + i, dof_start + j) += k3x3(block_offset + i, block_offset + j);
                }
            }
        }
    }
    
    // But the above doesn't match NASTRAN's approach. Let's implement the exact NASTRAN method:
    // NASTRAN divides the quad into 4 triangles sharing a common pivot node
    // and only computes 3 triangles (excluding the one opposite the pivot)
    
    // Find pivot node (we'll use node 1 as pivot for simplicity, as NASTRAN does)
    int pivot_node = 1;
    
    // Determine which triangle to skip (JNOT)
    // In NASTRAN: if pivot is node 1, JNOT = 3; if pivot is node 2, JNOT = 4; etc.
    // For pivot node 1: JNOT = 3 (triangle 3: nodes 3,4,2)
    int jnot = 3;
    
    // Recompute with proper triangle selection
    K_global.setZero();
    
    // The four triangles in order: (1,2,4), (2,3,1), (3,4,2), (4,1,3)
    // Skip triangle jnot (1-indexed, so skip index jnot-1)
    for (int tri_idx = 0; tri_idx < 4; ++tri_idx) {
        if (tri_idx == jnot - 1) continue; // skip JNOT
        
        const auto& tri = triangles[tri_idx];
        
        // Get coordinates for this triangle
        std::vector<Real> coords(6);
        coords[0] = nodes[tri[0]-1][0]; coords[1] = nodes[tri[0]-1][1];
        coords[2] = nodes[tri[1]-1][0]; coords[3] = nodes[tri[1]-1][1];
        coords[4] = nodes[tri[2]-1][0]; coords[5] = nodes[tri[2]-1][1];
        
        // Compute triangle stiffness
        Eigen::Matrix<Real, 3, 3> k3x3;
        ktrmem(coords, E, nu, t, sinth, costh, k3x3);
        
        // Assemble: for each node in triangle, add its contribution
        // NASTRAN stores three 3x3 matrices in k3x3: indices 1-9, 10-18, 19-27
        // We'll extract the three 2x2 blocks
        for (int node_in_tri = 0; node_in_tri < 3; ++node_in_tri) {
            int global_node = tri[node_in_tri];
            int block_start = node_in_tri * 3;
            
            int dof_start = 2 * (global_node - 1);
            for (int i = 0; i < 2; ++i) {
                for (int j = 0; j < 2; ++j) {
                    K_global(dof_start + i, dof_start + j) += k3x3(block_start + i, block_start + j);
                }
            }
        }
    }
    
    // Also need to add off-diagonal terms - the full 6x6 triangle stiffness contributes to coupling
    // Let's do a more accurate implementation: compute full 6x6 for each triangle and assemble properly
    
    K_global.setZero();
    
    // For each triangle, compute full 6x6 and assemble
    for (int tri_idx = 0; tri_idx < 4; ++tri_idx) {
        if (tri_idx == jnot - 1) continue;
        
        const auto& tri = triangles[tri_idx];
        
        // Get coordinates
        std::vector<Real> coords(6);
        coords[0] = nodes[tri[0]-1][0]; coords[1] = nodes[tri[0]-1][1];
        coords[2] = nodes[tri[1]-1][0]; coords[3] = nodes[tri[1]-1][1];
        coords[4] = nodes[tri[2]-1][0]; coords[5] = nodes[tri[2]-1][1];
        
        // Compute full 6x6 stiffness for triangle
        Real x1t = coords[0], y1t = coords[1];
        Real x2t = coords[2], y2t = coords[3];
        Real x3t = coords[4], y3t = coords[5];
        
        Real area = 0.5 * std::abs((x2t-x1t)*(y3t-y1t) - (x3t-x1t)*(y2t-y1t));
        if (area == 0.0) continue;
        
        Real denom = 1.0 - nu*nu;
        Real d11 = E / denom;
        Real d12 = E * nu / denom;
        Real d33 = E * (1.0 - nu) / (2.0 * denom);
        
        Real b1 = y2t - y3t;
        Real b2 = y3t - y1t;
        Real b3 = y1t - y2t;
        Real c1 = x3t - x2t;
        Real c2 = x1t - x3t;
        Real c3 = x2t - x1t;
        
        Real inv_2A = 1.0 / (2.0 * area);
        
        // B matrix (3x6)
        Eigen::Matrix<Real, 3, 6> B;
        B.setZero();
        B(0,0) = b1 * inv_2A; B(0,2) = b2 * inv_2A; B(0,4) = b3 * inv_2A;
        B(1,1) = c1 * inv_2A; B(1,3) = c2 * inv_2A; B(1,5) = c3 * inv_2A;
        B(2,0) = c1 * inv_2A; B(2,1) = b1 * inv_2A;
        B(2,2) = c2 * inv_2A; B(2,3) = b2 * inv_2A;
        B(2,4) = c3 * inv_2A; B(2,5) = b3 * inv_2A;
        
        // D matrix (3x3)
        Eigen::Matrix<Real, 3, 3> D;
        D << d11, d12, 0.0,
             d12, d11, 0.0,
             0.0, 0.0, d33;
        
        // k = t * area * B' * D * B
        Eigen::Matrix<Real, 6, 6> k6x6 = t * area * B.transpose() * D * B;
        
        // Assemble into 8x8 global matrix
        // Map triangle nodes to global DOFs
        for (int i = 0; i < 3; ++i) {
            int node_i = tri[i];
            int dof_i_start = 2 * (node_i - 1);
            for (int j = 0; j < 3; ++j) {
                int node_j = tri[j];
                int dof_j_start = 2 * (node_j - 1);
                // Add 2x2 block from k6x6
                for (int di = 0; di < 2; ++di) {
                    for (int dj = 0; dj < 2; ++dj) {
                        int k6x6_i = i * 2 + di;
                        int k6x6_j = j * 2 + dj;
                        K_global(dof_i_start + di, dof_j_start + dj) += k6x6(k6x6_i, k6x6_j);
                    }
                }
            }
        }
    }
    
    return K_global;
}

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<Real, 8, 8>& mat) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 8; ++i) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 8; ++j) {
            if (j > 0) oss << ",";
            // Use scientific notation with 15 digits precision
            oss << std::scientific << std::setprecision(15) << mat(i, j);
        }
        oss << "]";
    }
    
    oss << "]}";
    return oss.str();
}

int main() {
    try {
        // Test case: Nodes (0,0), (2,0), (2,1.5), (0,1.5)
        std::vector<std::vector<Real>> nodes = {
            {0.0, 0.0},
            {2.0, 0.0},
            {2.0, 1.5},
            {0.0, 1.5}
        };
        
        // Material properties
        Real E = 200e9;      // Pa
        Real nu = 0.3;       // Poisson's ratio
        Real t = 0.01;       // thickness in meters
        
        // Compute stiffness matrix
        Eigen::Matrix<Real, 8, 8> K = compute_cquad4_stiffness(nodes, E, nu, t);
        
        // Output as JSON
        std::cout << matrix_to_json(K) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}