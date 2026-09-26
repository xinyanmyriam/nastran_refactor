#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Simplified version of MAT subroutine - just returns E and G for given material ID
void MAT(int mat_id, double& E, double& G) {
    // For test case: assume material ID 1 corresponds to steel-like properties
    // In real NASTRAN, this would look up in material database
    if (mat_id == 1) {
        E = 200.0e9;  // Pa
        G = 76.923e9; // Pa
    } else {
        E = 200.0e9;
        G = 76.923e9;
    }
}

int main() {
    // Test case parameters
    const double xA = 0.0, yA = 0.0, zA = 0.0;
    const double xB = 2.0, yB = 0.0, zB = 0.0;
    const double E = 200.0e9;      // Pa
    const double A = 0.01;         // m^2
    const double G = 76.923e9;     // Pa
    const double J = 5.0e-6;       // m^4 (polar moment of inertia)
    
    // Simulate ECPT array (17 elements) as described in Fortran
    // ECPT(1): element ID = 1
    // ECPT(2): grid point A ID = 1
    // ECPT(3): grid point B ID = 2
    // ECPT(4): material ID = 1
    // ECPT(5): area A = 0.01
    // ECPT(6): polar moment J = 5e-6
    // ECPT(7): torsional stress coeff (not used) = 0.0
    // ECPT(8): non-structural mass (not used) = 0.0
    // ECPT(9): coord sys ID for A = 0 (basic)
    // ECPT(10): xA = 0.0
    // ECPT(11): yA = 0.0
    // ECPT(12): zA = 0.0
    // ECPT(13): coord sys ID for B = 0 (basic)
    // ECPT(14): xB = 2.0
    // ECPT(15): yB = 0.0
    // ECPT(16): zB = 0.0
    // ECPT(17): element temperature = 0.0
    std::vector<double> ECPT(17);
    ECPT[0] = 1.0;   // element ID
    ECPT[1] = 1.0;   // grid point A ID
    ECPT[2] = 2.0;   // grid point B ID
    ECPT[3] = 1.0;   // material ID
    ECPT[4] = A;     // area
    ECPT[5] = J;     // polar moment
    ECPT[6] = 0.0;   // torsional stress coeff
    ECPT[7] = 0.0;   // non-structural mass
    ECPT[8] = 0.0;   // coord sys ID for A (0 = basic)
    ECPT[9] = xA;    // xA
    ECPT[10] = yA;   // yA
    ECPT[11] = zA;   // zA
    ECPT[12] = 0.0;  // coord sys ID for B (0 = basic)
    ECPT[13] = xB;   // xB
    ECPT[14] = yB;   // yB
    ECPT[15] = zB;   // zB
    ECPT[16] = 0.0;  // temperature

    // Simulate common block variables — use non-conflicting names
    int NPVT = 1;           // pivot node ID (node A)
    int node_B = 2;         // non-pivot node ID (node B) — renamed from 'J'
    // bool HEAT = false;     // unused
    // int IOPT4 = 0;         // unused
    // double DAMPC = 0.0;    // unused
    // int K4GGSW = 0;       // unused

    // Compute length XL
    double X, Y, Z, XL;
    X = ECPT[9] - ECPT[13];  // xA - xB
    Y = ECPT[10] - ECPT[14]; // yA - yB
    Z = ECPT[11] - ECPT[15]; // zA - zB
    XL = std::sqrt(X*X + Y*Y + Z*Z);

    // Defensive: avoid division by zero
    if (XL == 0.0) {
        std::cerr << "Error: Zero length rod element." << std::endl;
        return 1;
    }

    // Normalize direction vector
    Eigen::Vector3d XN;
    XN << X/XL, Y/XL, Z/XL;

    // Get material properties
    double E_val, G_val;
    MAT(static_cast<int>(ECPT[3]), E_val, G_val);

    // Compute DSCL = A * E / XL and DSCR = J * G / XL
    double DSCL = ECPT[4] * E_val / XL;
    double DSCR = ECPT[5] * G_val / XL;

    // Build N matrix (3x3 outer product of XN)
    Eigen::Matrix3d N;
    N << XN(0)*XN(0), XN(0)*XN(1), XN(0)*XN(2),
         XN(1)*XN(0), XN(1)*XN(1), XN(1)*XN(2),
         XN(2)*XN(0), XN(2)*XN(1), XN(2)*XN(2);

    // Build 6x6 stiffness matrix KE as per Fortran KROD logic
    Eigen::Matrix<double, 6, 6> K_6x6 = Eigen::Matrix<double, 6, 6>::Zero();

    // Top-left 3x3 block: DSCL * N
    K_6x6.block<3,3>(0,0) = DSCL * N;
    // Top-right 3x3 block: -DSCL * N
    K_6x6.block<3,3>(0,3) = -DSCL * N;
    // Bottom-left 3x3 block: -DSCL * N
    K_6x6.block<3,3>(3,0) = -DSCL * N;
    // Bottom-right 3x3 block: DSCL * N + DSCR * N
    K_6x6.block<3,3>(3,3) = DSCL * N + DSCR * N;

    // Convert to JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << std::scientific << std::setprecision(6) << K_6x6(i,j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}