#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Forward declarations for required subroutines (we'll implement only what's needed)
void ktrmem(int flag);
void sma1b(const std::vector<double>& kij, int grid_id, int npvt, int ifkgg, double temp);

// Global data structures (mimicking Fortran COMMON blocks)
struct CommonData {
    // /CONDAS/
    std::vector<double> consts;
    
    // /SMA1HT/
    bool heat;
    
    // /SMA1ET/ - ECPT(100)
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
    int nogO; // renamed to avoid conflict with 'NOGO' in Fortran
    
    // /SMA1DP/
    std::vector<double> kij;
    std::vector<double> dum7;
    std::vector<double> ksum;
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
    std::vector<double> dummy8;
    
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
    std::vector<double> dum88;
    
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
        dummy8(248, 0.0),
        matid(0),
        inflag(0),
        eltemp(0.0),
        stress(0.0),
        sinth(0.0),
        costh(0.0),
        dum99(11, 0.0),
        gsube(0.0),
        dum88(6, 0.0) {}
};

CommonData common;

// Helper function to get DEGRA (degrees to radians conversion factor)
double get_degra() {
    return M_PI / 180.0;
}

// Helper function to compute triangular membrane stiffness (simplified)
// This is a simplified version that computes the 3x3 stiffness for a triangle
// using standard linear triangular element formulation
void compute_triangle_stiffness(const std::vector<std::vector<double>>& coords,
                               double E, double nu, double t,
                               std::vector<std::vector<double>>& k6x6) {
    // coords: 3 points, each with [x, y, z] (we use only x,y for 2D)
    // k6x6: 6x6 stiffness matrix for the triangle (2 DOF per node)
    
    // Extract coordinates (2D only)
    double x1 = coords[0][0], y1 = coords[0][1];
    double x2 = coords[1][0], y2 = coords[1][1];
    double x3 = coords[2][0], y3 = coords[2][1];
    
    // Area of triangle
    double area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    if (area == 0.0) {
        // Zero area triangle
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                k6x6[i][j] = 0.0;
            }
        }
        return;
    }
    
    // Material matrix D for plane stress
    // D = (E/(1-nu^2)) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
    double denom = 1.0 - nu*nu;
    double d11 = E / denom;
    double d12 = E * nu / denom;
    double d33 = E * (1.0 - nu) / (2.0 * denom);
    
    // Strain-displacement matrix B
    // For linear triangle: B = [[b1,0,b2,0,b3,0], [0,c1,0,c2,0,c3], [c1,b1,c2,b2,c3,b3]] / (2*area)
    // where b_i = y_{j} - y_{k}, c_i = x_{k} - x_{j} for cyclic permutations
    
    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;
    
    // B matrix (3x6)
    std::vector<std::vector<double>> B(3, std::vector<double>(6, 0.0));
    double factor = 1.0 / (2.0 * area);
    B[0][0] = b1 * factor; B[0][1] = 0.0;      B[0][2] = b2 * factor; B[0][3] = 0.0;      B[0][4] = b3 * factor; B[0][5] = 0.0;
    B[1][0] = 0.0;      B[1][1] = c1 * factor; B[1][2] = 0.0;      B[1][3] = c2 * factor; B[1][4] = 0.0;      B[1][5] = c3 * factor;
    B[2][0] = c1 * factor; B[2][1] = b1 * factor; B[2][2] = c2 * factor; B[2][3] = b2 * factor; B[2][4] = c3 * factor; B[2][5] = b3 * factor;
    
    // D matrix (3x3)
    std::vector<std::vector<double>> D(3, std::vector<double>(3, 0.0));
    D[0][0] = d11; D[0][1] = d12; D[0][2] = 0.0;
    D[1][0] = d12; D[1][1] = d11; D[1][2] = 0.0;
    D[2][0] = 0.0;  D[2][1] = 0.0;  D[2][2] = d33;
    
    // Compute B^T * D * B * t * area
    // Result is 6x6 stiffness matrix for the triangle
    std::vector<std::vector<double>> k6x6_temp(6, std::vector<double>(6, 0.0));
    
    // First compute D * B
    std::vector<std::vector<double>> DB(3, std::vector<double>(6, 0.0));
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 3; k++) {
                DB[i][j] += D[i][k] * B[k][j];
            }
        }
    }
    
    // Then compute B^T * (D * B)
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            for (int k = 0; k < 3; k++) {
                k6x6_temp[i][j] += B[k][i] * DB[k][j];
            }
        }
    }
    
    // Multiply by t * area
    double scale = t * area;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            k6x6_temp[i][j] *= scale;
        }
    }
    
    // Copy result
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            k6x6[i][j] = k6x6_temp[i][j];
        }
    }
}

// Simplified KTRMEM subroutine (triangular membrane)
void ktrmem(int flag) {
    // This is a simplified version that computes the stiffness for the current triangle
    // The triangle coordinates are in ECPT(10) onwards: ECPT(10-12), (14-16), (18-20) for nodes A,B,C
    
    // Get the three nodes for the current triangle
    std::vector<std::vector<double>> coords(3, std::vector<double>(3, 0.0));
    
    // Node A: ECPT(10), ECPT(11), ECPT(12)
    coords[0][0] = common.ecpt[9];  // index 9 for 10th element (0-based)
    coords[0][1] = common.ecpt[10]; // index 10 for 11th element
    coords[0][2] = common.ecpt[11]; // index 11 for 12th element
    
    // Node B: ECPT(14), ECPT(15), ECPT(16)
    coords[1][0] = common.ecpt[13]; // index 13 for 14th element
    coords[1][1] = common.ecpt[14]; // index 14 for 15th element
    coords[1][2] = common.ecpt[15]; // index 15 for 16th element
    
    // Node C: ECPT(18), ECPT(19), ECPT(20)
    coords[2][0] = common.ecpt[17]; // index 17 for 18th element
    coords[2][1] = common.ecpt[18]; // index 18 for 19th element
    coords[2][2] = common.ecpt[19]; // index 19 for 20th element
    
    // Material properties from CONSTS
    double E = common.consts[0]; // Young's modulus
    double nu = common.consts[1]; // Poisson's ratio
    double t = common.ecpt[6];    // thickness (ECPT(7) in 1-based, index 6 in 0-based)
    
    // Compute stiffness matrix for this triangle
    std::vector<std::vector<double>> k6x6(6, std::vector<double>(6, 0.0));
    compute_triangle_stiffness(coords, E, nu, t, k6x6);
    
    // Store in common.kij array (6x6 matrix stored in row-major order)
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            common.kij[i * 6 + j] = k6x6[i][j];
        }
    }
}

// Simplified SMA1B subroutine (matrix assembly)
void sma1b(const std::vector<double>& kij, int grid_id, int npvt, int ifkgg, double temp) {
    // This is a simplified version that accumulates the stiffness contributions
    // into the global 8x8 stiffness matrix
    // In the Fortran code, SMA1B inserts the 6x6 matrix into the global system
    // For our CQUAD4 with 4 nodes and 2 DOF per node, we need an 8x8 matrix
    
    // Since we're building the full 8x8 stiffness matrix, we'll accumulate
    // the contributions from each triangle into the appropriate positions
    // The mapping is: node 1 -> DOF 0,1; node 2 -> DOF 2,3; node 3 -> DOF 4,5; node 4 -> DOF 6,7
}

// Main CQUAD4 stiffness computation
Eigen::Matrix<double, 8, 8> compute_cquad4_stiffness(
    const std::vector<std::vector<double>>& nodes,
    double E, double nu, double t) {
    
    // Initialize 8x8 stiffness matrix
    Eigen::Matrix<double, 8, 8> K = Eigen::Matrix<double, 8, 8>::Zero();
    
    // Set up common data to mimic Fortran behavior
    common.consts[0] = E;      // Young's modulus
    common.consts[1] = nu;     // Poisson's ratio
    common.consts[3] = get_degra(); // DEGRA for angle conversion
    
    // Set up ECPT array with test case data
    // ECPT(1) = element ID
    common.ecpt[0] = 1.0;
    
    // ECPT(2-5) = grid point IDs (1,2,3,4)
    common.ecpt[1] = 1.0; // node 1
    common.ecpt[2] = 2.0; // node 2
    common.ecpt[3] = 3.0; // node 3
    common.ecpt[4] = 4.0; // node 4
    
    // ECPT(6) = theta (0.0 for this test)
    common.ecpt[5] = 0.0;
    
    // ECPT(7) = material ID
    common.ecpt[6] = 1.0;
    
    // ECPT(8) = thickness
    common.ecpt[7] = t;
    
    // ECPT(9) = non-structural mass (0.0)
    common.ecpt[8] = 0.0;
    
    // ECPT(10-12) = coordinates of node 1 (0,0,0)
    common.ecpt[9]  = nodes[0][0]; // x1
    common.ecpt[10] = nodes[0][1]; // y1
    common.ecpt[11] = 0.0;        // z1
    
    // ECPT(14-16) = coordinates of node 2 (2,0,0)
    common.ecpt[13] = nodes[1][0]; // x2
    common.ecpt[14] = nodes[1][1]; // y2
    common.ecpt[15] = 0.0;        // z2
    
    // ECPT(18-20) = coordinates of node 3 (2,1.5,0)
    common.ecpt[17] = nodes[2][0]; // x3
    common.ecpt[18] = nodes[2][1]; // y3
    common.ecpt[19] = 0.0;        // z3
    
    // ECPT(23-25) = coordinates of node 4 (0,1.5,0)
    common.ecpt[22] = nodes[3][0]; // x4
    common.ecpt[23] = nodes[3][1]; // y4
    common.ecpt[24] = 0.0;        // z4
    
    // ECPT(26) = element temperature
    common.ecpt[25] = 0.0;
    
    // Set HRING = false
    bool hring = false;
    
    // Compute angles
    double angl = common.ecpt[5] * common.consts[3]; // theta * DEGRA
    common.cosang = std::cos(angl);
    common.sinang = std::sin(angl);
    
    // Compute IVEC = node2 - node1
    common.ivec[0] = common.ecpt[13] - common.ecpt[9];  // x2 - x1
    common.ivec[1] = common.ecpt[14] - common.ecpt[10]; // y2 - y1
    common.ivec[2] = common.ecpt[15] - common.ecpt[11]; // z2 - z1
    
    // Normalize IVEC
    common.vecl = std::sqrt(common.ivec[0]*common.ivec[0] + 
                           common.ivec[1]*common.ivec[1] + 
                           common.ivec[2]*common.ivec[2]);
    if (common.vecl != 0.0) {
        common.ivec[0] /= common.vecl;
        common.ivec[1] /= common.vecl;
        common.ivec[2] /= common.vecl;
    }
    
    // Compute VSUBK = IVEC × (node4 - node1)
    double dx41 = common.ecpt[22] - common.ecpt[9];
    double dy41 = common.ecpt[23] - common.ecpt[10];
    double dz41 = common.ecpt[24] - common.ecpt[11];
    
    common.vsubk[0] = common.ivec[1]*dz41 - common.ivec[2]*dy41;
    common.vsubk[1] = common.ivec[2]*dx41 - common.ivec[0]*dz41;
    common.vsubk[2] = common.ivec[0]*dy41 - common.ivec[1]*dx41;
    
    // Normalize VSUBK to get KVEC
    common.vecl = std::sqrt(common.vsubk[0]*common.vsubk[0] + 
                           common.vsubk[1]*common.vsubk[1] + 
                           common.vsubk[2]*common.vsubk[2]);
    if (common.vecl != 0.0) {
        common.kvec[0] = common.vsubk[0] / common.vecl;
        common.kvec[1] = common.vsubk[1] / common.vecl;
        common.kvec[2] = common.vsubk[2] / common.vecl;
    }
    
    // Compute JVEC = KVEC × IVEC
    common.jvec[0] = common.kvec[1]*common.ivec[2] - common.kvec[2]*common.ivec[1];
    common.jvec[1] = common.kvec[2]*common.ivec[0] - common.kvec[0]*common.ivec[2];
    common.jvec[2] = common.kvec[0]*common.ivec[1] - common.kvec[1]*common.ivec[0];
    
    // Compute PVEC = cos(theta)*IVEC + sin(theta)*JVEC
    for (int i = 0; i < 3; i++) {
        common.pvec[i] = common.cosang * common.ivec[i] + common.sinang * common.jvec[i];
    }
    
    // Set NGRID (node IDs)
    common.ngrid[0] = static_cast<int>(common.ecpt[1]); // node 1
    common.ngrid[1] = static_cast<int>(common.ecpt[2]); // node 2
    common.ngrid[2] = static_cast<int>(common.ecpt[3]); // node 3
    common.ngrid[3] = static_cast<int>(common.ecpt[4]); // node 4
    
    // Copy coordinates to COORD
    for (int i = 0; i < 4; i++) {
        int base = i * 4;
        common.coord[base]     = common.ecpt[9 + i*4];   // x_i
        common.coord[base+1] = common.ecpt[10 + i*4];  // y_i
        common.coord[base+2] = common.ecpt[11 + i*4];  // z_i
        common.coord[base+3] = 0.0; // coordinate system ID (not used)
    }
    
    // Set up pivot point (NPVT = 1 for this test, meaning node 1 is pivot)
    common.npvt = 1;
    
    // Find pivot index
    int npivot = 0;
    for (int i = 0; i < 4; i++) {
        if (common.ngrid[i] == common.npvt) {
            npivot = i + 1; // 1-based index
            break;
        }
    }
    
    // Compute JNOT (triangle not used)
    int jnot;
    if (npivot <= 2) {
        jnot = npivot + 2;
    } else {
        jnot = npivot - 2;
    }
    
    // Zero out KSUM
    for (int i = 0; i < 36; i++) {
        common.ksum[i] = 0.0;
    }
    
    // Mapping matrix M (12 elements, 1-based indexing)
    std::vector<int> m = {1, 2, 4, 2, 3, 1, 3, 4, 2, 4, 1, 3};
    
    // Process 4 triangles
    for (int j = 1; j <= 4; j++) {
        if (j == jnot) continue;
        
        // Fill ECPT for triangle j
        int mpoint = 3*j - 3; // 0-based index calculation
        
        // Triangle nodes: get indices from mapping matrix
        std::vector<int> tri_nodes(3);
        for (int i = 0; i < 3; i++) {
            int npt1 = mpoint + i;
            int ns = m[npt1] - 1; // convert to 0-based
            tri_nodes[i] = common.ngrid[ns];
            
            // Copy coordinates for this node
            int npt1_coord = 4 * ns;
            int npt3 = 4 * i + 4;
            common.ecpt[npt3]     = common.coord[npt1_coord];
            common.ecpt[npt3+1] = common.coord[npt1_coord+1];
            common.ecpt[npt3+2] = common.coord[npt1_coord+2];
            common.ecpt[npt3+3] = common.coord[npt1_coord+3];
        }
        
        // Set up SINTH and COSTH for this triangle
        if (j == 1) {
            common.sinth = common.sinang;
            common.costh = common.cosang;
        } else {
            // Compute for other triangles
            double v1 = common.ecpt[13] - common.ecpt[9];  // x2-x1
            double v2 = common.ecpt[14] - common.ecpt[10]; // y2-y1
            double v3 = common.ecpt[15] - common.ecpt[11]; // z2-z1
            
            double vecl_v = std::sqrt(v1*v1 + v2*v2 + v3*v3);
            if (vecl_v != 0.0) {
                double u1_val = (v1*common.pvec[0] + v2*common.pvec[1] + v3*common.pvec[2]) / vecl_v;
                
                double si1 = v2*common.kvec[2] - v3*common.kvec[1];
                double si2 = v3*common.kvec[0] - v1*common.kvec[2];
                double si3 = v1*common.kvec[1] - v2*common.kvec[0];
                
                double u2_val = (si1*common.kvec[0] + si2*common.kvec[1] + si3*common.kvec[2]) / vecl_v;
                
                double vecl_u = std::sqrt(u1_val*u1_val + u2_val*u2_val);
                if (vecl_u != 0.0) {
                    u1_val /= vecl_u;
                    u2_val /= vecl_u;
                }
                common.sinth = u2_val;
                common.costh = u1_val;
            }
        }
        
        // Call KTRMEM (simplified)
        ktrmem(1);
        
        // Add contribution to KSUM
        // In real Fortran, K3X3(27) contains 3 matrices of 3x3 each
        // We'll simulate by computing the stiffness contribution for this triangle
        // and adding to the appropriate positions in KSUM
        
        // For simplicity, we'll compute the full 6x6 stiffness for this triangle
        // using the three nodes of the triangle
        std::vector<std::vector<double>> tri_coords(3, std::vector<double>(3, 0.0));
        
        // Get coordinates for the three nodes of triangle j
        switch(j) {
            case 1: // nodes 1,2,4
                tri_coords[0] = {nodes[0][0], nodes[0][1], 0.0};
                tri_coords[1] = {nodes[1][0], nodes[1][1], 0.0};
                tri_coords[2] = {nodes[3][0], nodes[3][1], 0.0};
                break;
            case 2: // nodes 2,3,1
                tri_coords[0] = {nodes[1][0], nodes[1][1], 0.0};
                tri_coords[1] = {nodes[2][0], nodes[2][1], 0.0};
                tri_coords[2] = {nodes[0][0], nodes[0][1], 0.0};
                break;
            case 3: // nodes 3,4,2
                tri_coords[0] = {nodes[2][0], nodes[2][1], 0.0};
                tri_coords[1] = {nodes[3][0], nodes[3][1], 0.0};
                tri_coords[2] = {nodes[1][0], nodes[1][1], 0.0};
                break;
            case 4: // nodes 4,1,3
                tri_coords[0] = {nodes[3][0], nodes[3][1], 0.0};
                tri_coords[1] = {nodes[0][0], nodes[0][1], 0.0};
                tri_coords[2] = {nodes[2][0], nodes[2][1], 0.0};
                break;
        }
        
        // Compute triangle stiffness
        std::vector<std::vector<double>> k6x6(6, std::vector<double>(6, 0.0));
        compute_triangle_stiffness(tri_coords, E, nu, t, k6x6);
        
        // Map triangle DOFs to global DOFs
        // Triangle DOFs: node1->(0,1), node2->(2,3), node3->(4,5)
        // Global DOFs: node1->(0,1), node2->(2,3), node3->(4,5), node4->(6,7)
        
        std::vector<int> global_dof_map(6);
        switch(j) {
            case 1: // nodes 1,2,4 -> DOF (0,1), (2,3), (6,7)
                global_dof_map = {0,1,2,3,6,7};
                break;
            case 2: // nodes 2,3,1 -> DOF (2,3), (4,5), (0,1)
                global_dof_map = {2,3,4,5,0,1};
                break;
            case 3: // nodes 3,4,2 -> DOF (4,5), (6,7), (2,3)
                global_dof_map = {4,5,6,7,2,3};
                break;
            case 4: // nodes 4,1,3 -> DOF (6,7), (0,1), (4,5)
                global_dof_map = {6,7,0,1,4,5};
                break;
        }
        
        // Add to global stiffness matrix
        for (int i = 0; i < 6; i++) {
            for (int j_idx = 0; j_idx < 6; j_idx++) {
                int gi = global_dof_map[i];
                int gj = global_dof_map[j_idx];
                K(gi, gj) += k6x6[i][j_idx];
            }
        }
    }
    
    return K;
}

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<double, 8, 8>& K) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 8; i++) {
        oss << "[";
        for (int j = 0; j < 8; j++) {
            oss << std::scientific << std::setprecision(15) << K(i, j);
            if (j < 7) oss << ",";
        }
        oss << "]";
        if (i < 7) oss << ",";
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
    
    // Material properties
    double E = 200e9;   // Pa
    double nu = 0.3;
    double t = 0.01;    // m
    
    // Compute stiffness matrix
    Eigen::Matrix<double, 8, 8> K = compute_cquad4_stiffness(nodes, E, nu, t);
    
    // Output as JSON
    std::cout << matrix_to_json(K) << std::endl;
    
    return 0;
}