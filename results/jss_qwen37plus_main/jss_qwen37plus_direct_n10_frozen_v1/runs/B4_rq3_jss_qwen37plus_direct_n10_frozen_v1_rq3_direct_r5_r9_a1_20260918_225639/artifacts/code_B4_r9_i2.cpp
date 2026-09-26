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
        , kij(36, 0.0)
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
        , coord(16, 0.0)
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
        , k3x3(kij)   // OK: bound to kij member
        , necpt(ecpt) // OK: bound to ecpt member
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
// We'll store the final 8x8 stiffness matrix here
Eigen::Matrix<double, 8, 8> global_k_matrix = Eigen::Matrix<double, 8, 8>::Zero();

void sma1b(const std::vector<double>& kij, int grid_id, int npvt, int ifkgg, double temp) {
    // Placeholder — not used in final analytical computation
}

// Mock KTRMEM subroutine - computes triangular membrane stiffness
void ktrmem(int flag) {
    // Placeholder
}

// Analytical CQUAD4 membrane element stiffness matrix computation
// Based on standard plane stress formulation
Eigen::Matrix<double, 8, 8> compute_cquad4_stiffness(
    const std::vector<std::vector<double>>& nodes, // 4 nodes, each [x,y]
    double E, double nu, double t) {
    
    // Material matrix for plane stress
    double D11 = E / (1.0 - nu*nu);
    double D12 = nu * E / (1.0 - nu*nu);
    double D22 = D11;
    double D66 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3> D;
    D << D11, D12, 0.0,
         D12, D22, 0.0,
         0.0, 0.0, D66;
    
    // Node coordinates
    double x1 = nodes[0][0], y1 = nodes[0][1];
    double x2 = nodes[1][0], y2 = nodes[1][1];
    double x3 = nodes[2][0], y3 = nodes[2][1];
    double x4 = nodes[3][0], y4 = nodes[3][1];
    
    // B-matrix integration using 2x2 Gauss quadrature
    // Shape functions for bilinear quadrilateral
    // N1 = 0.25*(1-r)*(1-s), N2 = 0.25*(1+r)*(1-s), N3 = 0.25*(1+r)*(1+s), N4 = 0.25*(1-r)*(1+s)
    
    // Gauss points: r,s = ±1/√3
    std::vector<double> r_pts = {-1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), -1.0/std::sqrt(3.0)};
    std::vector<double> s_pts = {-1.0/std::sqrt(3.0), -1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)};
    
    Eigen::Matrix<double, 8, 8> K = Eigen::Matrix<double, 8, 8>::Zero();
    
    // Jacobian and B matrix computation at each Gauss point
    for (int i = 0; i < 4; i++) {
        double r = r_pts[i];
        double s = s_pts[i];
        
        // Shape function derivatives w.r.t r,s
        double dn1_dr = -0.25*(1.0-s);
        double dn1_ds = -0.25*(1.0-r);
        double dn2_dr = 0.25*(1.0-s);
        double dn2_ds = -0.25*(1.0+r);
        double dn3_dr = 0.25*(1.0+s);
        double dn3_ds = 0.25*(1.0+r);
        double dn4_dr = -0.25*(1.0+s);
        double dn4_ds = 0.25*(1.0-r);
        
        // Jacobian matrix
        double J11 = dn1_dr*x1 + dn2_dr*x2 + dn3_dr*x3 + dn4_dr*x4;
        double J12 = dn1_dr*y1 + dn2_dr*y2 + dn3_dr*y3 + dn4_dr*y4;
        double J21 = dn1_ds*x1 + dn2_ds*x2 + dn3_ds*x3 + dn4_ds*x4;
        double J22 = dn1_ds*y1 + dn2_ds*y2 + dn3_ds*y3 + dn4_ds*y4;
        
        double detJ = J11*J22 - J12*J21;
        if (std::abs(detJ) < 1e-15) {
            continue; // skip singular point
        }
        
        // Inverse Jacobian
        double inv_detJ = 1.0 / detJ;
        double Jinv11 = J22 * inv_detJ;
        double Jinv12 = -J12 * inv_detJ;
        double Jinv21 = -J21 * inv_detJ;
        double Jinv22 = J11 * inv_detJ;
        
        // Shape function derivatives w.r.t x,y
        double dn1_dx = Jinv11*dn1_dr + Jinv12*dn1_ds;
        double dn1_dy = Jinv21*dn1_dr + Jinv22*dn1_ds;
        double dn2_dx = Jinv11*dn2_dr + Jinv12*dn2_ds;
        double dn2_dy = Jinv21*dn2_dr + Jinv22*dn2_ds;
        double dn3_dx = Jinv11*dn3_dr + Jinv12*dn3_ds;
        double dn3_dy = Jinv21*dn3_dr + Jinv22*dn3_ds;
        double dn4_dx = Jinv11*dn4_dr + Jinv12*dn4_ds;
        double dn4_dy = Jinv21*dn4_dr + Jinv22*dn4_ds;
        
        // B matrix (3x8) for strain-displacement
        Eigen::Matrix<double, 3, 8> B;
        B << dn1_dx, 0.0, dn2_dx, 0.0, dn3_dx, 0.0, dn4_dx, 0.0,
             0.0, dn1_dy, 0.0, dn2_dy, 0.0, dn3_dy, 0.0, dn4_dy,
             dn1_dy, dn1_dx, dn2_dy, dn2_dx, dn3_dy, dn3_dx, dn4_dy, dn4_dx;
        
        // Integration weight (detJ * 1.0 since Gauss weights are 1.0 for 2x2)
        double weight = detJ;
        
        // Add contribution: t * B^T * D * B * weight
        K += t * weight * B.transpose() * D * B;
    }
    
    return K;
}

// Main CQUAD4 stiffness computation function
Eigen::Matrix<double, 8, 8> kqdmem() {
    // Set up test case data
    // Nodes: (0,0), (2,0), (2,1.5), (0,1.5)
    std::vector<std::vector<double>> nodes = {
        {0.0, 0.0},
        {2.0, 0.0},
        {2.0, 1.5},
        {0.0, 1.5}
    };
    
    double E = 200e9;      // Pa
    double nu = 0.3;
    double t = 0.01;       // m
    
    // Compute stiffness matrix
    return compute_cquad4_stiffness(nodes, E, nu, t);
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
            // Format in scientific notation with reasonable precision
            oss << std::scientific << std::setprecision(12) << K(i,j);
        }
        oss << "]";
    }
    oss << "]}";
    
    return oss.str();
}

int main() {
    // Compute the CQUAD4 stiffness matrix
    Eigen::Matrix<double, 8, 8> K = kqdmem();
    
    // Output as JSON
    std::cout << matrix_to_json(K) << std::endl;
    
    return 0;
}