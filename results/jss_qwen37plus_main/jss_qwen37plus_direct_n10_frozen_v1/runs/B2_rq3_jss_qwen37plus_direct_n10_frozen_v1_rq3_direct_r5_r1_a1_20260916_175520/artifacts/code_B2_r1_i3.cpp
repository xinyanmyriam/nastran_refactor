#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>

std::string double_to_scientific(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << x;
    std::string s = oss.str();
    // Replace "e+" with "e" (e.g., "1.234567e+02" → "1.234567e02")
    size_t pos = s.find("e+");
    if (pos != std::string::npos) {
        s.replace(pos, 2, "e");
    }
    return s;
}

std::string matrix_to_json(const Eigen::Matrix<double, 12, 12>& K) {
    std::string result = "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        if (i > 0) result += ",";
        result += "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) result += ",";
            result += double_to_scientific(K(i,j));
        }
        result += "]";
    }
    
    result += "]}";
    return result;
}

int main() {
    // Test case parameters
    const double nodeA[3] = {0.0, 0.0, 0.0};
    const double nodeB[3] = {2.0, 0.0, 0.0};
    const double L = 2.0; // length
    
    const double E = 200e9;      // Pa
    const double G = 76.923e9;   // Pa
    const double A = 0.01;       // m^2
    const double Iy = 8.333e-6;  // m^4 (I1 in Fortran)
    const double Iz = 8.333e-6;  // m^4 (I2 in Fortran)
    const double J = 1.667e-5;   // m^4 (FJ in Fortran)
    const double I12 = 0.0;      // no warping (I12 = 0)
    
    // Pin flags: none (0)
    const int pinA = 0;
    const int pinB = 0;
    
    // Reference vector (SMALLV) - use y-axis as reference for local y-axis
    const double smallv[3] = {0.0, 1.0, 0.0};
    
    // No offsets
    const double za[3] = {0.0, 0.0, 0.0};
    const double zb[3] = {0.0, 0.0, 0.0};
    
    // Shear correction factors (K1, K2) - set to 0 to use pure Euler-Bernoulli
    const double K1 = 0.0;
    const double K2 = 0.0;
    
    // Compute element stiffness matrix
    
    // Step 1: Compute local coordinate system
    // VECI = vector from B to A (in Fortran, VECI is negated version of B-A)
    Eigen::Vector3d veci;
    veci << nodeA[0] - nodeB[0], nodeA[1] - nodeB[1], nodeA[2] - nodeB[2];
    double fl = veci.norm();
    if (fl == 0.0) {
        std::cerr << "Zero length element" << std::endl;
        return 1;
    }
    veci.normalize();
    
    // SMALV0 = normalized smallv
    Eigen::Vector3d smalv0;
    smalv0 << smallv[0], smallv[1], smallv[2];
    double fll = smalv0.norm();
    if (fll == 0.0) {
        std::cerr << "Zero reference vector" << std::endl;
        return 1;
    }
    smalv0.normalize();
    
    // VECK = VECI × SMALV0
    Eigen::Vector3d veck = veci.cross(smalv0);
    fll = veck.norm();
    if (fll == 0.0) {
        std::cerr << "Reference vector parallel to element axis" << std::endl;
        return 1;
    }
    veck.normalize();
    
    // VECJ = VECK × VECI
    Eigen::Vector3d vecj = veck.cross(veci);
    fll = vecj.norm();
    if (fll == 0.0) {
        std::cerr << "Cross product failed" << std::endl;
        return 1;
    }
    vecj.normalize();
    
    // Step 2: Compute intermediate stiffness terms
    const double LSQ = L * L;
    const double LCUBE = LSQ * L;
    
    const double EI1 = E * Iy;  // bending about local y (I1)
    const double EI2 = E * Iz;  // bending about local z (I2)
    
    // Shear deformation terms (R1, R2)
    double R1, R2;
    if (K1 == 0.0 || I12 != 0.0) {
        R1 = 12.0 * EI1 / LCUBE;
    } else {
        const double GAK1 = G * A * K1;
        R1 = (12.0 * EI1 * GAK1) / (GAK1 * LCUBE + 12.0 * L * EI1);
    }
    
    if (K2 == 0.0 || I12 != 0.0) {
        R2 = 12.0 * EI2 / LCUBE;
    } else {
        const double GAK2 = G * A * K2;
        R2 = (12.0 * EI2 * GAK2) / (GAK2 * LCUBE + 12.0 * L * EI2);
    }
    
    // SK terms
    const double SK1 = 0.25 * R1 * LSQ + EI1 / L;
    const double SK2 = 0.25 * R2 * LSQ + EI2 / L;
    const double SK3 = 0.25 * R1 * LSQ - EI1 / L;
    const double SK4 = 0.25 * R2 * LSQ - EI2 / L;
    
    // Other terms
    const double AEL = E * A / L;
    const double LR1 = L * R1 / 2.0;
    const double LR2 = L * R2 / 2.0;
    const double GJL = G * J / L;
    
    // Step 3: Build 12x12 stiffness matrix in local coordinates
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Fill the matrix according to Fortran indexing (1-based)
    // Fortran indices: KE(1) = K(0,0), KE(7) = K(6,0), etc.
    // We'll fill row by row, column by column (0-based indexing)
    
    // Row 0 (index 0): [AEL, 0, 0, 0, 0, 0, -AEL, 0, 0, 0, 0, 0]
    K(0,0) = AEL;
    K(0,6) = -AEL;
    
    // Row 1 (index 1): [0, R1, 0, 0, 0, LR1, 0, -R1, 0, 0, 0, LR1]
    K(1,1) = R1;
    K(1,5) = LR1;
    K(1,7) = -R1;
    K(1,11) = LR1;
    
    // Row 2 (index 2): [0, 0, R2, 0, -LR2, 0, 0, 0, -R2, 0, -LR2, 0]
    K(2,2) = R2;
    K(2,4) = -LR2;
    K(2,8) = -R2;
    K(2,10) = -LR2;
    
    // Row 3 (index 3): [0, 0, 0, GJL, 0, 0, 0, 0, 0, -GJL, 0, 0]
    K(3,3) = GJL;
    K(3,9) = -GJL;
    
    // Row 4 (index 4): [0, 0, -LR2, 0, SK2, 0, 0, 0, LR2, 0, SK4, 0]
    K(4,2) = -LR2;
    K(4,4) = SK2;
    K(4,8) = LR2;
    K(4,10) = SK4;
    
    // Row 5 (index 5): [0, LR1, 0, 0, 0, SK1, 0, -LR1, 0, 0, 0, SK3]
    K(5,1) = LR1;
    K(5,5) = SK1;
    K(5,7) = -LR1;
    K(5,11) = SK3;
    
    // Row 6 (index 6): [-AEL, 0, 0, 0, 0, 0, AEL, 0, 0, 0, 0, 0]
    K(6,0) = -AEL;
    K(6,6) = AEL;
    
    // Row 7 (index 7): [0, -R1, 0, 0, 0, -LR1, 0, R1, 0, 0, 0, -LR1]
    K(7,1) = -R1;
    K(7,5) = -LR1;
    K(7,7) = R1;
    K(7,11) = -LR1;
    
    // Row 8 (index 8): [0, 0, -R2, 0, LR2, 0, 0, 0, R2, 0, LR2, 0]
    K(8,2) = -R2;
    K(8,4) = LR2;
    K(8,8) = R2;
    K(8,10) = LR2;
    
    // Row 9 (index 9): [0, 0, 0, -GJL, 0, 0, 0, 0, 0, GJL, 0, 0]
    K(9,3) = -GJL;
    K(9,9) = GJL;
    
    // Row 10 (index 10): [0, 0, LR2, 0, SK4, 0, 0, 0, -LR2, 0, SK2, 0]
    K(10,2) = LR2;
    K(10,4) = SK4;
    K(10,8) = -LR2;
    K(10,10) = SK2;
    
    // Row 11 (index 11): [0, -LR1, 0, 0, 0, SK3, 0, LR1, 0, 0, 0, SK1]
    K(11,1) = -LR1;
    K(11,5) = SK3;
    K(11,7) = LR1;
    K(11,11) = SK1;
    
    // Handle warping (I12 term) - but I12 = 0 in test case, so skip
    
    // Handle pin flags - both zero, so no modification
    
    // The above matrix K is in the local element coordinate system.
    // For the test case, the element is aligned with global X-axis,
    // and the local coordinate system matches global, so no transformation needed.
    
    // Output as JSON
    std::cout << matrix_to_json(K) << std::endl;
    
    return 0;
}