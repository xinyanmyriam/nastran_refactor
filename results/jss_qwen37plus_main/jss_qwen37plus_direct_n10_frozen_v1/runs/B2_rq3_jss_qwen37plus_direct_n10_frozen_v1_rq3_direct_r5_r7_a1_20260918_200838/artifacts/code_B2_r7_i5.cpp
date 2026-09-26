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
    size_t e_pos = s.find('e');
    if (e_pos != std::string::npos) {
        // Remove trailing zeros after decimal point
        size_t dot_pos = s.find('.');
        if (dot_pos != std::string::npos && dot_pos < e_pos) {
            // Find last non-zero digit before exponent
            size_t last_nonzero = e_pos;
            for (size_t i = e_pos; i > dot_pos; --i) {
                if (s[i] != '0' && s[i] != '.') {
                    last_nonzero = i;
                    break;
                }
            }
            if (last_nonzero < e_pos) {
                s.erase(last_nonzero + 1, e_pos - last_nonzero);
                // Ensure we keep at least one digit after decimal
                if (s.back() == '.') {
                    s.pop_back();
                }
            }
        }
    }
    // Remove '+' from exponent
    size_t plus_pos = s.find("e+");
    if (plus_pos != std::string::npos) {
        s.replace(plus_pos + 1, 1, "");
    }
    return s;
}

// Main function to compute CBAR stiffness matrix
Eigen::Matrix<double, 12, 12> compute_cbar_stiffness() {
    // Test case parameters
    const double nodeA[3] = {0.0, 0.0, 0.0};
    const double nodeB[3] = {2.0, 0.0, 0.0};
    
    const double E = 200e9;           // Pa
    const double G = 76.923e9;        // Pa
    const double A = 0.01;            // m^2
    const double Iy = 8.333e-6;       // m^4 (I1 in Fortran)
    const double Iz = 8.333e-6;       // m^4 (I2 in Fortran)
    const double J = 1.667e-5;        // m^4 (FJ in Fortran)
    const double I12 = 0.0;           // No warping (I12 = 0)
    
    // Reference vector (SMALLV) - use default [0,0,1] for z-axis
    const double smallv[3] = {0.0, 0.0, 1.0};
    
    // Material properties
    const double K1 = 0.0;  // Shear correction factor (not used in test case)
    const double K2 = 0.0;  // Shear correction factor (not used in test case)
    
    // Compute element length L
    double dx = nodeB[0] - nodeA[0];
    double dy = nodeB[1] - nodeA[1];
    double dz = nodeB[2] - nodeA[2];
    double L = std::sqrt(dx*dx + dy*dy + dz*dz);
    
    // Normalize reference vector
    double fl = std::sqrt(smallv[0]*smallv[0] + smallv[1]*smallv[1] + smallv[2]*smallv[2]);
    double smalv0[3] = {smallv[0]/fl, smallv[1]/fl, smallv[2]/fl};
    
    // Compute VECI = vector from A to B (note: corrected — Fortran uses A->B for x-axis)
    // Earlier code used B->A; that gave -x. We want x_local = A->B.
    double veci[3] = {dx, dy, dz};
    fl = std::sqrt(veci[0]*veci[0] + veci[1]*veci[1] + veci[2]*veci[2]);
    
    // Normalize VECI
    veci[0] /= fl;
    veci[1] /= fl;
    veci[2] /= fl;
    
    // Compute VECK = VECI × SMALV0
    double veck[3];
    veck[0] = veci[1]*smalv0[2] - veci[2]*smalv0[1];
    veck[1] = veci[2]*smalv0[0] - veci[0]*smalv0[2];
    veck[2] = veci[0]*smalv0[1] - veci[1]*smalv0[0];
    
    double fll = std::sqrt(veck[0]*veck[0] + veck[1]*veck[1] + veck[2]*veck[2]);
    if (fll < 1e-12) {
        // fallback: if veci || smalv0, pick arbitrary orthogonal vector
        if (std::abs(veci[0]) > 1e-12 || std::abs(veci[1]) > 1e-12) {
            veck[0] = -veci[1]; veck[1] = veci[0]; veck[2] = 0.0;
        } else {
            veck[0] = 1.0; veck[1] = 0.0; veck[2] = 0.0;
        }
        fll = std::sqrt(veck[0]*veck[0] + veck[1]*veck[1] + veck[2]*veck[2]);
    }
    veck[0] /= fll;
    veck[1] /= fll;
    veck[2] /= fll;
    
    // Compute VECJ = VECK × VECI
    double vecj[3];
    vecj[0] = veck[1]*veci[2] - veck[2]*veci[1];
    vecj[1] = veck[2]*veci[0] - veck[0]*veci[2];
    vecj[2] = veck[0]*veci[1] - veck[1]*veci[0];
    
    fll = std::sqrt(vecj[0]*vecj[0] + vecj[1]*vecj[1] + vecj[2]*vecj[2]);
    if (fll < 1e-12) {
        // fallback
        vecj[0] = veci[1]; vecj[1] = -veci[0]; vecj[2] = 0.0;
        fll = std::sqrt(vecj[0]*vecj[0] + vecj[1]*vecj[1] + vecj[2]*vecj[2]);
    }
    vecj[0] /= fll;
    vecj[1] /= fll;
    vecj[2] /= fll;
    
    // Build 3x3 direction cosine matrix T3: local -> global
    // Rows = local x,y,z axes expressed in global coordinates
    Eigen::Matrix<double, 3, 3> T3;
    T3 << veci[0], veci[1], veci[2],
          vecj[0], vecj[1], vecj[2],
          veck[0], veck[1], veck[2];

    // Intermediate variables
    double LSQ = L * L;
    double LCUBE = LSQ * L;
    
    // EI1 = E * Iy, EI2 = E * Iz
    double EI1 = E * Iy;
    double EI2 = E * Iz;
    
    // R1 and R2 calculations (shear deformation terms)
    double R1, R2;
    if (K1 == 0.0 || I12 != 0.0) {
        R1 = 12.0 * EI1 / LCUBE;
    } else {
        double GAK1 = G * A * K1;
        R1 = (12.0 * EI1 * GAK1) / (GAK1 * LCUBE + 12.0 * L * EI1);
    }
    
    if (K2 == 0.0 || I12 != 0.0) {
        R2 = 12.0 * EI2 / LCUBE;
    } else {
        double GAK2 = G * A * K2;
        R2 = (12.0 * EI2 * GAK2) / (GAK2 * LCUBE + 12.0 * L * EI2);
    }
    
    // SK1, SK2, SK3, SK4
    double SK1 = 0.25 * R1 * LSQ + EI1 / L;
    double SK2 = 0.25 * R2 * LSQ + EI2 / L;
    double SK3 = 0.25 * R1 * LSQ - EI1 / L;
    double SK4 = 0.25 * R2 * LSQ - EI2 / L;
    
    // Other terms
    double AEL = E * A / L;
    double LR1 = L * R1 / 2.0;
    double LR2 = L * R2 / 2.0;
    double GJL = G * J / L;
    
    // Initialize 12x12 LOCAL stiffness matrix
    Eigen::Matrix<double, 12, 12> K_local = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Fill the local stiffness matrix (standard CBAR in local frame: x along bar)
    // Row 0 (index 0): [AEL, 0, 0, 0, 0, 0, -AEL, 0, 0, 0, 0, 0]
    K_local(0,0) = AEL;
    K_local(0,6) = -AEL;
    
    // Row 1 (index 1): [0, R1, 0, 0, 0, LR1, 0, -R1, 0, 0, 0, LR1]
    K_local(1,1) = R1;
    K_local(1,5) = LR1;
    K_local(1,7) = -R1;
    K_local(1,11) = LR1;
    
    // Row 2 (index 2): [0, 0, R2, 0, -LR2, 0, 0, 0, -R2, 0, LR2, 0]
    K_local(2,2) = R2;
    K_local(2,4) = -LR2;
    K_local(2,8) = -R2;
    K_local(2,10) = LR2;
    
    // Row 3 (index 3): [0, 0, 0, GJL, 0, 0, 0, 0, 0, -GJL, 0, 0]
    K_local(3,3) = GJL;
    K_local(3,9) = -GJL;
    
    // Row 4 (index 4): [0, 0, -LR2, 0, SK2, 0, 0, 0, LR2, 0, SK4, 0]
    K_local(4,2) = -LR2;
    K_local(4,4) = SK2;
    K_local(4,8) = LR2;
    K_local(4,10) = SK4;
    
    // Row 5 (index 5): [0, LR1, 0, 0, 0, SK1, 0, -LR1, 0, 0, 0, SK3]
    K_local(5,1) = LR1;
    K_local(5,5) = SK1;
    K_local(5,7) = -LR1;
    K_local(5,11) = SK3;
    
    // Row 6 (index 6): [-AEL, 0, 0, 0, 0, 0, AEL, 0, 0, 0, 0, 0]
    K_local(6,0) = -AEL;
    K_local(6,6) = AEL;
    
    // Row 7 (index 7): [0, -R1, 0, 0, 0, -LR1, 0, R1, 0, 0, 0, -LR1]
    K_local(7,1) = -R1;
    K_local(7,5) = -LR1;
    K_local(7,7) = R1;
    K_local(7,11) = -LR1;
    
    // Row 8 (index 8): [0, 0, -R2, 0, LR2, 0, 0, 0, R2, 0, -LR2, 0]
    K_local(8,2) = -R2;
    K_local(8,4) = LR2;
    K_local(8,8) = R2;
    K_local(8,10) = -LR2;
    
    // Row 9 (index 9): [0, 0, 0, -GJL, 0, 0, 0, 0, 0, GJL, 0, 0]
    K_local(9,3) = -GJL;
    K_local(9,9) = GJL;
    
    // Row 10 (index 10): [0, 0, LR2, 0, SK4, 0, 0, 0, -LR2, 0, SK2, 0]
    K_local(10,2) = LR2;
    K_local(10,4) = SK4;
    K_local(10,8) = -LR2;
    K_local(10,10) = SK2;
    
    // Row 11 (index 11): [0, -LR1, 0, 0, 0, SK3, 0, LR1, 0, 0, 0, SK1]
    K_local(11,1) = -LR1;
    K_local(11,5) = SK3;
    K_local(11,7) = LR1;
    K_local(11,11) = SK1;
    
    // Build 12x12 transformation matrix T12
    // T12 = blockdiag(T3, T3, T3, T3) — for [u1,v1,w1,rx1,ry1,rz1, u2,v2,w2,rx2,ry2,rz2]
    Eigen::Matrix<double, 12, 12> T12 = Eigen::Matrix<double, 12, 12>::Zero();
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            T12(i, j) = T3(i, j);           // node1 translation
            T12(i+3, j+3) = T3(i, j);       // node1 rotation
            T12(i+6, j+6) = T3(i, j);       // node2 translation
            T12(i+9, j+9) = T3(i, j);       // node2 rotation
        }
    }
    
    // Transform to global frame: K_global = T12^T * K_local * T12
    Eigen::Matrix<double, 12, 12> K_global = T12.transpose() * K_local * T12;
    
    return K_global;
}

int main() {
    // Compute the stiffness matrix
    Eigen::Matrix<double, 12, 12> K = compute_cbar_stiffness();
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << format_double(K(i,j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    
    std::cout << "]}" << std::endl;
    
    return 0;
}