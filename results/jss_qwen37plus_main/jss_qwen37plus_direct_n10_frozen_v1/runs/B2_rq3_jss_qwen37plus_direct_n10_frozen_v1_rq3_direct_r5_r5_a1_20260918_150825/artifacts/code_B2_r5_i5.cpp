#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Helper function to format double in scientific notation with fixed precision
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Remove '+' from exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos && s[epos + 1] == '+') {
        s.erase(epos + 1, 1);
    }
    return s;
}

// Main function to compute CBAR stiffness matrix
Eigen::Matrix<double, 12, 12> compute_cbar_stiffness() {
    // Test case parameters
    const double E = 200e9;           // Pa
    const double G = 76.923e9;        // Pa
    const double A = 0.01;            // m^2
    const double Iy = 8.333e-6;       // m^4 (I1 in Fortran)
    const double Iz = 8.333e-6;       // m^4 (I2 in Fortran)
    const double J = 1.667e-5;        // m^4 (FJ in Fortran)
    
    // Node coordinates
    Eigen::Vector3d nodeA(0.0, 0.0, 0.0);
    Eigen::Vector3d nodeB(2.0, 0.0, 0.0);
    
    // Compute element length L
    double L = (nodeB - nodeA).norm();
    
    // In the Fortran code, VECI is computed as -(nodeA - nodeB) = nodeB - nodeA
    // So our local x-axis points from A to B
    Eigen::Vector3d veci = (nodeB - nodeA).normalized();
    
    // Reference vector SMALLV: use y-axis as reference for plane
    Eigen::Vector3d smallv(0.0, 1.0, 0.0);
    double fl = smallv.norm();
    Eigen::Vector3d smalv0 = smallv / fl;
    
    // Compute veck = veci × smalv0
    Eigen::Vector3d veck = veci.cross(smalv0);
    double fll = veck.norm();
    if (fll == 0.0) {
        // Fallback: if veci is parallel to y-axis, use z-axis
        smalv0 = Eigen::Vector3d(0.0, 0.0, 1.0);
        veck = veci.cross(smalv0);
        fll = veck.norm();
    }
    veck /= fll;
    
    // Compute vecj = veck × veci
    Eigen::Vector3d vecj = veck.cross(veci);
    fll = vecj.norm();
    if (fll == 0.0) {
        // Should not happen with proper cross products
        vecj = Eigen::Vector3d(0.0, 1.0, 0.0);
        vecj -= vecj.dot(veci) * veci;
        fll = vecj.norm();
        if (fll > 1e-15) {
            vecj /= fll;
        } else {
            // Last resort: construct orthogonal basis
            if (std::abs(veci(0)) > std::abs(veci(1)) && std::abs(veci(0)) > std::abs(veci(2))) {
                vecj << 0.0, -veci(2), veci(1);
            } else if (std::abs(veci(1)) > std::abs(veci(2))) {
                vecj << -veci(2), 0.0, veci(0);
            } else {
                vecj << -veci(1), veci(0), 0.0;
            }
            fll = vecj.norm();
            vecj /= fll;
        }
    }
    
    // Now we have orthonormal basis: veci (x), vecj (y), veck (z)
    // Build transformation matrix T from local to global: columns are veci, vecj, veck
    Eigen::Matrix<double, 3, 3> T;
    T.col(0) = veci;
    T.col(1) = vecj;
    T.col(2) = veck;
    
    // Material properties
    double EI1 = E * Iy;  // bending about y-axis (local y -> global vecj)
    double EI2 = E * Iz;  // bending about z-axis (local z -> global veck)
    double AEL = E * A / L;
    double GJL = G * J / L;
    
    // Shear correction factors (K1, K2) - set to 1.0 as not provided in test case
    // Use basic Euler-Bernoulli (no shear)
    double LSQ = L * L;
    double LCUBE = LSQ * L;
    double R1 = 12.0 * EI1 / LCUBE;  // 12EI1/L^3
    double R2 = 12.0 * EI2 / LCUBE;  // 12EI2/L^3
    
    // Standard Euler-Bernoulli coefficients
    double LR1 = 6.0 * EI1 / LSQ;    // 6EI1/L^2
    double LR2 = 6.0 * EI2 / LSQ;    // 6EI2/L^2
    double SK1 = 4.0 * EI1 / L;      // 4EI1/L
    double SK2 = 4.0 * EI2 / L;      // 4EI2/L
    double SK3 = 2.0 * EI1 / L;      // 2EI1/L
    double SK4 = 2.0 * EI2 / L;      // 2EI2/L
    
    // Initialize 12x12 stiffness matrix (in local coordinates)
    Eigen::Matrix<double, 12, 12> K_local = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Fill the local stiffness matrix according to standard CBAR formulation
    // DOF ordering: [TX1, TY1, TZ1, RX1, RY1, RZ1, TX2, TY2, TZ2, RX2, RY2, RZ2]
    
    // Axial (TX1-TX2)
    K_local(0,0) = AEL;
    K_local(0,6) = -AEL;
    K_local(6,0) = -AEL;
    K_local(6,6) = AEL;
    
    // Y-bending (TY1, RZ1, TY2, RZ2) - uses Iy (EI1)
    K_local(1,1) = R1;
    K_local(1,5) = LR1;
    K_local(1,7) = -R1;
    K_local(1,11) = LR1;
    
    K_local(5,1) = LR1;
    K_local(5,5) = SK1;
    K_local(5,7) = -LR1;
    K_local(5,11) = SK3;
    
    K_local(7,1) = -R1;
    K_local(7,5) = -LR1;
    K_local(7,7) = R1;
    K_local(7,11) = -LR1;
    
    K_local(11,1) = LR1;
    K_local(11,5) = SK3;
    K_local(11,7) = -LR1;
    K_local(11,11) = SK1;
    
    // Z-bending (TZ1, RY1, TZ2, RY2) - uses Iz (EI2)
    K_local(2,2) = R2;
    K_local(2,4) = -LR2;
    K_local(2,8) = -R2;
    K_local(2,10) = LR2;
    
    K_local(4,2) = -LR2;
    K_local(4,4) = SK2;
    K_local(4,8) = LR2;
    K_local(4,10) = SK4;
    
    K_local(8,2) = -R2;
    K_local(8,4) = LR2;
    K_local(8,8) = R2;
    K_local(8,10) = -LR2;
    
    K_local(10,2) = LR2;
    K_local(10,4) = SK4;
    K_local(10,8) = -LR2;
    K_local(10,10) = SK2;
    
    // Torsion (RX1, RX2)
    K_local(3,3) = GJL;
    K_local(3,9) = -GJL;
    K_local(9,3) = -GJL;
    K_local(9,9) = GJL;
    
    // Now transform to global coordinates
    // The transformation matrix for a beam element is:
    // [T  0  0  0]
    // [0  T  0  0]
    // [0  0  T  0]
    // [0  0  0  T]
    // where T is the 3x3 rotation matrix
    
    Eigen::Matrix<double, 12, 12> T_global = Eigen::Matrix<double, 12, 12>::Zero();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T_global(i, j) = T(i, j);
            T_global(i+3, j+3) = T(i, j);
            T_global(i+6, j+6) = T(i, j);
            T_global(i+9, j+9) = T(i, j);
        }
    }
    
    // K_global = T_global^T * K_local * T_global
    Eigen::Matrix<double, 12, 12> K_global = T_global.transpose() * K_local * T_global;
    
    return K_global;
}

int main() {
    // Compute the stiffness matrix
    Eigen::Matrix<double, 12, 12> K = compute_cbar_stiffness();
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << format_double(K(i, j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}