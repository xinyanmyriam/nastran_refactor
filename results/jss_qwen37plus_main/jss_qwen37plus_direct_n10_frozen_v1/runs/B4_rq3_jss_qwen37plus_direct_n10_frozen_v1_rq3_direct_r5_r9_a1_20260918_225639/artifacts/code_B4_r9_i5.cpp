#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Forward declarations
void ktrmem(int flag);
void sma1b(const std::vector<double>& kij, int grid_id, int npvt, int ifkgg, double temp);

// Global data structures (mimicking Fortran COMMON blocks)
struct CommonData {
    // /CONDAS/
    std::vector<double> consts;

    // /SMA1HT/
    bool heat;

    // /SMA1ET/
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
    int noGo; // renamed from nogO

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
    std::vector<double> dum88; // now unique name

    // /K3X3/ (equivalenced to kij) — REFERENCE → must be initialized in ctor
    std::vector<double>& k3x3;

    // /NECPT/ (equivalenced to ecpt) — REFERENCE → must be initialized in ctor
    std::vector<double>& necpt;

    // Constructor initializes ALL members, especially references
    explicit CommonData()
        : consts(5, 0.0)
        , heat(false)
        , ecpt(100, 0.0)
        , dum1(10, 0.0)
        , ifkgg(0)
        , dum2(1, 0.0)
        , if4gg(0)
        , dum3(23, 0.0)
        , iopt4(0)
        , k4ggsw(0)
        , npvt(0)
        , dumcl(7, 0.0)
        , link(10, 0)
        , idetck(0)
        , dodet(false)
        , noGo(0)
        , kij(36, 0.0)  // 6x6 for triangle
        , dum7(156, 0.0)
        , ksum(36, 0.0)
        , temp(0.0)
        , cosang(0.0)
        , sinang(0.0)
        , vecl(0.0)
        , ivec(3, 0.0)
        , jvec(3, 0.0)
        , kvec(3, 0.0)
        , pvec(3, 0.0)
        , vsubk(3, 0.0)
        , v(3, 0.0)
        , si(3, 0.0)
        , npivot(0)
        , mpoint(0)
        , mi(0)
        , nsubsc(0)
        , ngrid(4, 0)
        , u1(0.0)
        , u2(0.0)
        , coord(16, 0.0)  // 4 nodes × 2 coords = 8, but pad to 16 for safety
        , dummy88(248, 0.0)
        , matid(0)
        , inflag(0)
        , eltemp(0.0)
        , stress(0.0)
        , sinth(0.0)
        , costh(0.0)
        , dum99(11, 0.0)
        , gsube(0.0)
        , dum88(6, 0.0)
        , k3x3(kij)   // bound to kij
        , necpt(ecpt) // bound to ecpt
    {}
};

CommonData common;

// Helper function to compute cross product
std::vector<double> cross(const std::vector<double>& a, const std::vector<double>& b) {
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

// Helper function to compute dot product
double dot(const std::vector<double>& a, const std::vector<double>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Helper function to compute vector norm
double norm(const std::vector<double>& v) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Helper function to normalize vector
std::vector<double> normalize(const std::vector<double>& v) {
    double n = norm(v);
    if (n == 0.0) return {0.0, 0.0, 0.0};
    return {v[0]/n, v[1]/n, v[2]/n};
}

// Mock MESAGE subroutine (just sets error flag)
void mesage(int code, int subcode, double elid) {
    common.noGo = 1;
}

// Mock SMA1B subroutine - accumulates stiffness contributions into global K matrix
Eigen::Matrix<double, 8, 8> global_k_matrix = Eigen::Matrix<double, 8, 8>::Zero();

void sma1b(const std::vector<double>& kij, int grid_id, int npvt, int ifkgg, double temp) {
    // We'll do assembly inline in kqdmem; this is just a placeholder
}

// Compute 6x6 CST triangle stiffness in LOCAL coordinate system (x' along edge 0->1, y' in-plane)
// Nodes: p0, p1, p2 in local 2D (x',y') — all z=0
// Area = |(p1-p0) × (p2-p0)| / 2
// B = [[dN1/dx' 0 dN2/dx' 0 dN3/dx' 0],
//      [0 dN1/dy' 0 dN2/dy' 0 dN3/dy'],
//      [dN1/dy' dN1/dx' dN2/dy' dN2/dx' dN3/dy' dN3/dx']]  (3x6)
// D = [[D11,D12,0],[D12,D22,0],[0,0,D66]]
// K = t * A * B^T * D * B
void ktrmem(int flag) {
    // We assume flag==1 means compute and store in common.kij (6x6)
    if (flag != 1) return;

    // Get triangle nodes from common.coord: 4 nodes, but we use only 3
    // For now, assume triangle is nodes 0,1,2 (first triangle)
    // common.coord layout: [x0,y0, x1,y1, x2,y2, x3,y3] → 8 elements
    std::vector<std::vector<double>> tri_nodes(3, std::vector<double>(2));
    for (int i = 0; i < 3; ++i) {
        tri_nodes[i][0] = common.coord[2*i];
        tri_nodes[i][1] = common.coord[2*i + 1];
    }

    // Compute vectors
    double x0 = tri_nodes[0][0], y0 = tri_nodes[0][1];
    double x1 = tri_nodes[1][0], y1 = tri_nodes[1][1];
    double x2 = tri_nodes[2][0], y2 = tri_nodes[2][1];

    double dx10 = x1 - x0, dy10 = y1 - y0;
    double dx20 = x2 - x0, dy20 = y2 - y0;

    double area = 0.5 * std::abs(dx10*dy20 - dy10*dx20);
    if (area < 1e-15) {
        for (int i = 0; i < 36; ++i) common.kij[i] = 0.0;
        return;
    }

    // Material properties
    double E = 200e9;      // Pa
    double nu = 0.3;
    double t = 0.01;       // m

    double D11 = E / (1.0 - nu*nu);
    double D12 = nu * E / (1.0 - nu*nu);
    double D22 = D11;
    double D66 = E / (2.0 * (1.0 + nu));

    // B matrix coefficients (constant for CST)
    // N1 = a1 + b1*x + c1*y, etc.
    // b1 = (y2-y3)/2A = (y2-y1)/2A, c1 = (x1-x2)/2A, etc.
    double inv2A = 1.0 / (2.0 * area);
    double b1 = (y2 - y1) * inv2A;
    double c1 = (x1 - x2) * inv2A;
    double b2 = (y0 - y2) * inv2A;
    double c2 = (x2 - x0) * inv2A;
    double b3 = (y1 - y0) * inv2A;
    double c3 = (x0 - x1) * inv2A;

    // B = [ b1 0 b2 0 b3 0;
    //       0 c1 0 c2 0 c3;
    //       c1 b1 c2 b2 c3 b3 ]
    Eigen::Matrix<double, 3, 6> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;

    Eigen::Matrix<double, 3, 3> D;
    D << D11, D12, 0.0,
         D12, D22, 0.0,
         0.0, 0.0, D66;

    // K_local = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6> K_local = t * area * B.transpose() * D * B;

    // Copy to common.kij (row-major)
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            common.kij[i*6 + j] = K_local(i,j);
        }
    }
}

// Build 6x6 rotation matrix for triangle: T6 = blockdiag(T2,T2,T2)
// where T2 = [[c, s], [-s, c]], c=cosang, s=sinang
Eigen::Matrix<double, 6, 6> build_T6x6(double c, double s) {
    Eigen::Matrix<double, 6, 6> T6 = Eigen::Matrix<double, 6, 6>::Zero();
    T6(0,0) = c; T6(0,1) = s;
    T6(1,0) = -s; T6(1,1) = c;

    T6(2,2) = c; T6(2,3) = s;
    T6(3,2) = -s; T6(3,3) = c;

    T6(4,4) = c; T6(4,5) = s;
    T6(5,4) = -s; T6(5,5) = c;

    return T6;
}

// Main CQUAD4 stiffness computation function — now using two triangles
Eigen::Matrix<double, 8, 8> kqdmem() {
    // Reset global matrix
    Eigen::Matrix<double, 8, 8> K = Eigen::Matrix<double, 8, 8>::Zero();

    // Set up nodes: (0,0), (2,0), (2,1.5), (0,1.5)
    std::vector<std::vector<double>> nodes = {
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.5},
        {0.0, 1.5}
    };

    // Store in common.coord (8 elements: x0,y0,x1,y1,x2,y2,x3,y3)
    for (int i = 0; i < 4; ++i) {
        common.coord[2*i]     = nodes[i][0];
        common.coord[2*i+1] = nodes[i][1];
    }

    // Material properties — set in common for consistency
    // (already hardcoded in ktrmem, but for clarity)
    double E = 200e9;
    double nu = 0.3;
    double t = 0.01;

    // Triangle 1: nodes 0,1,2 → global DOFs: [0,1, 2,3, 4,5]
    // Compute angle of edge 0->1
    double dx = nodes[1][0] - nodes[0][0];
    double dy = nodes[1][1] - nodes[0][1];
    double len = std::sqrt(dx*dx + dy*dy);
    if (len > 1e-15) {
        common.cosang = dx / len;
        common.sinang = dy / len;
    } else {
        common.cosang = 1.0;
        common.sinang = 0.0;
    }

    // Call ktrmem to fill common.kij with 6x6 local K
    ktrmem(1);

    // Build T6x6
    Eigen::Matrix<double, 6, 6> T6 = build_T6x6(common.cosang, common.sinang);

    // Convert common.kij to Eigen matrix
    Eigen::Matrix<double, 6, 6> K_tri_local;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K_tri_local(i,j) = common.kij[i*6 + j];
        }
    }

    // Rotate: K_tri_global = T6^T * K_tri_local * T6
    Eigen::Matrix<double, 6, 6> K_tri_global = T6.transpose() * K_tri_local * T6;

    // Assemble into 8x8: map triangle DOFs [0,1,2,3,4,5] → global [0,1,2,3,4,5]
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            int gi = i; // node0:0,1; node1:2,3; node2:4,5
            int gj = j;
            K(gi, gj) += K_tri_global(i,j);
        }
    }

    // Triangle 2: nodes 0,2,3 → global DOFs: [0,1, 4,5, 6,7]
    // Edge 0->2
    dx = nodes[2][0] - nodes[0][0];
    dy = nodes[2][1] - nodes[0][1];
    len = std::sqrt(dx*dx + dy*dy);
    if (len > 1e-15) {
        common.cosang = dx / len;
        common.sinang = dy / len;
    } else {
        common.cosang = 1.0;
        common.sinang = 0.0;
    }

    // Update common.coord to triangle 0,2,3
    // common.coord reused: set positions for p0,p2,p3
    // p0 = nodes[0], p2 = nodes[2], p3 = nodes[3]
    common.coord[0] = nodes[0][0]; common.coord[1] = nodes[0][1];
    common.coord[2] = nodes[2][0]; common.coord[3] = nodes[2][1];
    common.coord[4] = nodes[3][0]; common.coord[5] = nodes[3][1];
    // (we don't use index 6,7 for triangle)

    ktrmem(1);
    T6 = build_T6x6(common.cosang, common.sinang);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            K_tri_local(i,j) = common.kij[i*6 + j];
        }
    }
    K_tri_global = T6.transpose() * K_tri_local * T6;

    // Map triangle DOFs [0,1,2,3,4,5] → global [0,1, 4,5, 6,7]
    std::vector<int> map = {0,1,4,5,6,7};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            int gi = map[i];
            int gj = map[j];
            K(gi, gj) += K_tri_global(i,j);
        }
    }

    return K;
}

// Convert matrix to JSON string
std::string matrix_to_json(const Eigen::Matrix<double, 8, 8>& K) {
    std::ostringstream oss;
    oss << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 8; i++) {
        if (i > 0) oss << ",";
        oss << "[";
        for (int j = 0; j < 8; j++) {
            if (j > 0) oss << ",";
            oss << std::scientific << std::setprecision(12) << K(i,j);
        }
        oss << "]";
    }
    oss << "]}";
    
    return oss.str();
}

int main() {
    // Compute the CQUAD4 stiffness matrix via triangle assembly with rotation
    Eigen::Matrix<double, 8, 8> K = kqdmem();
    
    // Output as JSON
    std::cout << matrix_to_json(K) << std::endl;
    
    return 0;
}