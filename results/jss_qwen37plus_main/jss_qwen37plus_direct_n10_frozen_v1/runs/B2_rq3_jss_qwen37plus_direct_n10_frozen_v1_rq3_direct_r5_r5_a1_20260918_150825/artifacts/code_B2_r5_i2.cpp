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
    // In Fortran: if K1==0 or I12!=0, use R1 = 12*EI1/L^3, else compute with shear
    // Since test case doesn't specify K1/K2, use basic Euler-Bernoulli (no shear)
    double LSQ = L * L;
    double LCUBE = LSQ * L;
    double R1 = 12.0 * EI1 / LCUBE;
    double R2 = 12.0 * EI2 / LCUBE;
    
    // SK terms - FIXED: use 1/12 instead of 0.25 (which is 1/4)
    // Standard CBAR: SK1 = 4*EI1/L = EI1/L + R1*L^2/12
    //               SK2 = 4*EI2/L = EI2/L + R2*L^2/12  
    //               SK3 = 2*EI1/L = -EI1/L + R1*L^2/12
    //               SK4 = 2*EI2/L = -EI2/L + R2*L^2/12
    double SK1 = R1 * LSQ / 12.0 + EI1 / L;
    double SK2 = R2 * LSQ / 12.0 + EI2 / L;
    double SK3 = R1 * LSQ / 12.0 - EI1 / L;
    double SK4 = R2 * LSQ / 12.0 - EI2 / L;
    
    // LR1 and LR2 terms - these are correct: 6*EI/L^2 = R1*L/2
    double LR1 = L * R1 / 2.0;
    double LR2 = L * R2 / 2.0;
    
    // Initialize 12x12 stiffness matrix (in local coordinates)
    Eigen::Matrix<double, 12, 12> K_local = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Fill the local stiffness matrix according to Fortran indexing (1-based)
    // We'll fill row by row, column by column (0-based indexing)
    
    // Row 0 (index 0): [AEL, 0, 0, 0, 0, 0, -AEL, 0, 0, 0, 0, 0]
    K_local(0,0) = AEL;
    K_local(0,6) = -AEL;
    
    // Row 1: [0, R1, 0, 0, 0, LR1, 0, -R1, 0, 0, 0, LR1]
    K_local(1,1) = R1;
    K_local(1,5) = LR1;
    K_local(1,7) = -R1;
    K_local(1,11) = LR1;
    
    // Row 2: [0, 0, R2, 0, -LR2, 0, 0, 0, -R2, 0, LR2, 0]
    K_local(2,2) = R2;
    K_local(2,4) = -LR2;
    K_local(2,8) = -R2;
    K_local(2,10) = LR2;
    
    // Row 3: [0, 0, 0, GJL, 0, 0, 0, 0, 0, -GJL, 0, 0]
    K_local(3,3) = GJL;
    K_local(3,9) = -GJL;
    
    // Row 4: [0, 0, -LR2, 0, SK2, 0, 0, 0, LR2, 0, SK4, 0]
    K_local(4,2) = -LR2;
    K_local(4,4) = SK2;
    K_local(4,8) = LR2;
    K_local(4,10) = SK4;
    
    // Row 5: [0, LR1, 0, 0, 0, SK1, 0, -LR1, 0, 0, 0, SK3]
    K_local(5,1) = LR1;
    K_local(5,5) = SK1;
    K_local(5,7) = -LR1;
    K_local(5,11) = SK3;
    
    // Row 6: [-AEL, 0, 0, 0, 0, 0, AEL, 0, 0, 0, 0, 0]
    K_local(6,0) = -AEL;
    K_local(6,6) = AEL;
    
    // Row 7: [0, -R1, 0, 0, 0, -LR1, 0, R1, 0, 0, 0, -LR1]
    K_local(7,1) = -R1;
    K_local(7,5) = -LR1;
    K_local(7,7) = R1;
    K_local(7,11) = -LR1;
    
    // Row 8: [0, 0, -R2, 0, LR2, 0, 0, 0, R2, 0, -LR2, 0]
    K_local(8,2) = -R2;
    K_local(8,4) = LR2;
    K_local(8,8) = R2;
    K_local(8,10) = -LR2;
    
    // Row 9: [0, 0, 0, -GJL, 0, 0, 0, 0, 0, GJL, 0, 0]
    K_local(9,3) = -GJL;
    K_local(9,9) = GJL;
    
    // Row 10: [0, 0, LR2, 0, SK4, 0, 0, 0, -LR2, 0, SK2, 0]
    K_local(10,2) = LR2;
    K_local(10,4) = SK4;
    K_local(10,8) = -LR2;
    K_local(10,10) = SK2;
    
    // Row 11: [0, -LR1, 0, 0, 0, SK3, 0, LR1, 0, 0, 0, SK1]
    K_local(11,1) = -LR1;
    K_local(11,5) = SK3;
    K_local(11,7) = LR1;
    K_local(11,11) = SK1;
    
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
    
    // CORRECTION: K_global = T_global * K_local * T_global^T
    // because u_local = T_global^T * u_global, so K_global = T_global * K_local * T_global^T
    Eigen::Matrix<double, 12, 12> K_global = T_global * K_local * T_global.transpose();
    
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