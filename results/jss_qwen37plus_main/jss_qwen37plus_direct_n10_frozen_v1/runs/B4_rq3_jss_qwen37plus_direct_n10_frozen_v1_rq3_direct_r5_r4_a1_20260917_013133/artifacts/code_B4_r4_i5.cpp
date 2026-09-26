#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Define M_PI for Windows compatibility
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations
void ktrmem(int flag);
void sma1b(const Eigen::Matrix<double, 6, 6>& k6x6, const std::vector<int>& tri_nodes, int npvt, int ifkgg, double temp);

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
    int nogO; // renamed to avoid conflict with 'no' keyword
    
    // /SMA1DP/
    std::vector<double> kij;      // 36 elements for 6x6 stiffness
    std::vector<double> dum7;
    std::vector<double> ksum;     // 36 elements for 6x6 sum matrix storage
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
    std::vector<int> ngrid;       // 4 nodes: 1,2,3,4
    double u1;
    double u2;
    std::vector<double> coord;    // 16 elements: x,y,z for 4 nodes (each 3 values, but z=0 for membrane)
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
        dum88(6, 0.0)
    {}
    
    // Getter for k6x6 (equivalent to kij)
    std::vector<double>& get_k6x6() { return kij; }
    const std::vector<double>& get_k6x6() const { return kij; }
    
    // Getter for necpt (equivalent to ecpt)
    std::vector<double>& get_necpt() { return ecpt; }
    const std::vector<double>& get_necpt() const { return ecpt; }
};

CommonData common;

// Helper function to get DEGRA (degrees to radians conversion factor)
double get_degra() {
    return M_PI / 180.0;
}

// Helper function to compute triangular membrane stiffness
// Triangle nodes: ngrid[0], ngrid[1], ngrid[2]
void ktrmem(int flag) {
    // Get material properties
    double E = common.consts[0]; // Young's modulus
    double nu = common.consts[1]; // Poisson's ratio
    double t = common.ecpt[6];    // thickness (ECPT(7) in 1-based indexing)
    
    // Get coordinates of triangle vertices: ngrid[0], ngrid[1], ngrid[2]
    // coord layout: [x0,y0,z0, x1,y1,z1, x2,y2,z2, x3,y3,z3] → 4 nodes × 3 = 12 values
    // But we allocated 16 → safe to use first 12
    auto get_node_coord = [&](int idx) -> std::vector<double> {
        int base = idx * 3;
        return {common.coord[base], common.coord[base+1], common.coord[base+2]};
    };
    
    auto A = get_node_coord(common.ngrid[0]-1); // 1-based node ID → 0-based index
    auto B = get_node_coord(common.ngrid[1]-1);
    auto C = get_node_coord(common.ngrid[2]-1);
    
    double x1 = A[0], y1 = A[1];
    double x2 = B[0], y2 = B[1];
    double x3 = C[0], y3 = C[1];
    
    // Compute area of triangle
    double area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    if (area < 1e-15) {
        for (int i = 0; i < 36; i++) common.kij[i] = 0.0;
        return;
    }
    
    // Material matrix D for plane stress
    double denom = 1.0 - nu*nu;
    double d11 = E / denom;
    double d12 = E * nu / denom;
    double d33 = E * (1.0 - nu) / (2.0 * denom);
    
    Eigen::Matrix<double, 3, 3> D;
    D << d11, d12, 0.0,
         d12, d11, 0.0,
         0.0, 0.0, d33;
    
    // Shape function derivatives
    double b1 = (y2 - y3) / (2.0 * area);
    double c1 = (x3 - x2) / (2.0 * area);
    double b2 = (y3 - y1) / (2.0 * area);
    double c2 = (x1 - x3) / (2.0 * area);
    double b3 = (y1 - y2) / (2.0 * area);
    double c3 = (x2 - x1) / (2.0 * area);
    
    // B matrix (3x6)
    Eigen::Matrix<double, 3, 6> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    
    // Ke = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6> ke_full = t * area * B.transpose() * D * B;
    
    // Store in common.kij (row-major order: ke_full(i,j) → kij[i*6+j])
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            common.kij[i*6 + j] = ke_full(i, j);
        }
    }
}

// SMA1B: assemble 6x6 triangle stiffness into global 8x8 matrix
// tri_nodes: list of 3 node IDs (1-based) in this triangle
void sma1b(const Eigen::Matrix<double, 6, 6>& k6x6, const std::vector<int>& tri_nodes, int npvt, int ifkgg, double temp) {
    // Build 6x6 transformation matrix T6 = block_diag(R, R, R)
    // R = [[cosang, -sinang], [sinang, cosang]]
    Eigen::Matrix<double, 2, 2> R;
    R << common.cosang, -common.sinang,
         common.sinang,  common.cosang;
    
    Eigen::Matrix<double, 6, 6> T6 = Eigen::Matrix<double, 6, 6>::Zero();
    for (int i = 0; i < 3; i++) {
        T6.block<2,2>(i*2, i*2) = R;
    }
    
    // Rotate: k6x6_global = T6^T * k6x6 * T6
    Eigen::Matrix<double, 6, 6> k6x6_global = T6.transpose() * k6x6 * T6;
    
    // Map local triangle nodes to global DOF indices (0-based DOF index)
    std::vector<int> dof_indices;
    for (int node_id : tri_nodes) {
        // Find position of node_id in ngrid (1-based)
        int node_idx = -1;
        for (int i = 0; i < 4; i++) {
            if (common.ngrid[i] == node_id) {
                node_idx = i;
                break;
            }
        }
        if (node_idx == -1) continue;
        // Each node has 2 DOF: ux=2*node_idx, uy=2*node_idx+1
        dof_indices.push_back(node_idx * 2);
        dof_indices.push_back(node_idx * 2 + 1);
    }
    
    if (dof_indices.size() != 6) return;
    
    // Accumulate into global 8x8 matrix (static)
    static Eigen::Matrix<double, 8, 8> global_K = Eigen::Matrix<double, 8, 8>::Zero();
    
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            int gi = dof_indices[i];
            int gj = dof_indices[j];
            if (gi >= 0 && gi < 8 && gj >= 0 && gj < 8) {
                global_K(gi, gj) += k6x6_global(i, j);
            }
        }
    }
    
    // Store result for output
    static int call_count = 0;
    call_count++;
    if (call_count == 2) {
        // After both triangles assembled, output
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 8; i++) {
            std::cout << "[";
            for (int j = 0; j < 8; j++) {
                std::cout << std::scientific << std::setprecision(15) << global_K(i, j);
                if (j < 7) std::cout << ",";
            }
            std::cout << "]";
            if (i < 7) std::cout << ",";
        }
        std::cout << "]}" << std::endl;
    }
}

// Main driver
int main() {
    // Set up test mesh: CQUAD4 nodes (0,0), (2,0), (2,1.5), (0,1.5)
    std::vector<std::vector<double>> node_coords = {
        {0.0, 0.0, 0.0},
        {2.0, 0.0, 0.0},
        {2.0, 1.5, 0.0},
        {0.0, 1.5, 0.0}
    };
    
    // Populate common.ngrid (1-based node IDs)
    common.ngrid = {1, 2, 3, 4};
    
    // Populate common.coord: flatten [x0,y0,z0, x1,y1,z1, ...]
    common.coord.clear();
    for (const auto& nc : node_coords) {
        common.coord.push_back(nc[0]);
        common.coord.push_back(nc[1]);
        common.coord.push_back(nc[2]);
    }
    
    // Material properties: E=200e9, nu=0.3, t=0.01
    common.consts[0] = 200e9;   // E
    common.consts[1] = 0.3;     // nu
    common.ecpt[6] = 0.01;      // thickness t (ECPT(7))
    
    // Rotation angles: assume element local x-axis aligned with global x → cos=1, sin=0
    common.cosang = 1.0;
    common.sinang = 0.0;
    
    // Triangle 1: nodes 1,2,3
    common.ngrid[0] = 1; common.ngrid[1] = 2; common.ngrid[2] = 3;
    ktrmem(0);
    Eigen::Matrix<double, 6, 6> k1;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            k1(i, j) = common.kij[i*6 + j];
        }
    }
    sma1b(k1, {1,2,3}, common.npvt, common.ifkgg, common.temp);
    
    // Triangle 2: nodes 1,3,4
    common.ngrid[0] = 1; common.ngrid[1] = 3; common.ngrid[2] = 4;
    ktrmem(0);
    Eigen::Matrix<double, 6, 6> k2;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            k2(i, j) = common.kij[i*6 + j];
        }
    }
    sma1b(k2, {1,3,4}, common.npvt, common.ifkgg, common.temp);
    
    return 0;
}