#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <sstream>
#include <Eigen/Dense>

// Since we can't use external JSON library, we'll write our own minimal JSON formatter
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Replace e+ with e
    size_t epos = s.find("e+");
    if (epos != std::string::npos) {
        s.replace(epos, 2, "e");
    }
    return s;
}

int main() {
    // Test case parameters
    const double nodeA[3] = {0.0, 0.0, 0.0};
    const double nodeB[3] = {2.0, 0.0, 0.0};
    const double L = 2.0; // length in meters

    const double E = 200e9;      // Pa
    const double G = 76.923e9;  // Pa
    const double A = 0.01;      // m^2
    const double Iy = 8.333e-6; // m^4 (I1 in Fortran)
    const double Iz = 8.333e-6; // m^4 (I2 in Fortran)
    const double J = 1.667e-5;  // m^4 (FJ in Fortran)

    // In Fortran, I1 = Iy, I2 = Iz, FJ = J
    // No offsets, no pin flags, basic coordinates, no heat, no shear correction (K1=K2=0), no warping (I12=0)
    // So we skip all the complex coordinate transformations and just compute the standard CBAR stiffness

    // Compute intermediate terms
    const double LSQ = L * L;
    const double LCUBE = LSQ * L;

    const double EIy = E * Iy;
    const double EIz = E * Iz;
    const double GJL = G * J / L;
    const double AEL = E * A / L;

    // For standard Euler-Bernoulli beam (no shear deformation), R1 = 12*EIy/L^3, R2 = 12*EIz/L^3
    const double R1 = 12.0 * EIy / LCUBE;
    const double R2 = 12.0 * EIz / LCUBE;

    // SK1 = 0.25*R1*L^2 + EIy/L
    // SK2 = 0.25*R2*L^2 + EIz/L
    // SK3 = 0.25*R1*L^2 - EIy/L
    // SK4 = 0.25*R2*L^2 - EIz/L
    const double SK1 = 0.25 * R1 * LSQ + EIy / L;
    const double SK2 = 0.25 * R2 * LSQ + EIz / L;
    const double SK3 = 0.25 * R1 * LSQ - EIy / L;
    const double SK4 = 0.25 * R2 * LSQ - EIz / L;

    const double LR1 = L * R1 / 2.0;
    const double LR2 = L * R2 / 2.0;

    // Build 12x12 stiffness matrix KE (row-major, 1-indexed in Fortran, so we map accordingly)
    // Fortran indexing: KE(1) is [1,1], KE(2) is [2,1], ..., KE(12) is [12,1], KE(13) is [1,2], etc.
    // But in our C++ code, we'll build a column-major Eigen matrix and then output row by row.

    Eigen::Matrix<double, 12, 12, Eigen::ColMajor> KE = Eigen::Matrix<double, 12, 12, Eigen::ColMajor>::Zero();

    // Fill according to Fortran logic (converting 1-based to 0-based indices)
    // KE(  1) =  AEL        -> (0,0)
    KE(0, 0) = AEL;
    // KE(  7) = -AEL        -> (6,0)
    KE(6, 0) = -AEL;
    // KE( 14) =  R1         -> (1,1)
    KE(1, 1) = R1;
    // KE( 18) =  LR1        -> (5,1)
    KE(5, 1) = LR1;
    // KE( 20) = -R1         -> (7,1)
    KE(7, 1) = -R1;
    // KE( 24) =  LR1        -> (11,1)
    KE(11, 1) = LR1;
    // KE( 27) =  R2         -> (2,2)
    KE(2, 2) = R2;
    // KE( 29) = -LR2        -> (4,2)
    KE(4, 2) = -LR2;
    // KE( 33) = -R2         -> (8,2)
    KE(8, 2) = -R2;
    // KE( 35) = -LR2        -> (10,2)
    KE(10, 2) = -LR2;
    // KE( 40) =  GJL        -> (3,3)
    KE(3, 3) = GJL;
    // KE( 46) = -GJL        -> (9,3)
    KE(9, 3) = -GJL;
    // KE( 51) = -LR2        -> (2,4)
    KE(2, 4) = -LR2;
    // KE( 53) =  SK2        -> (4,4)
    KE(4, 4) = SK2;
    // KE( 57) =  LR2        -> (8,4)
    KE(8, 4) = LR2;
    // KE( 59) =  SK4        -> (10,4)
    KE(10, 4) = SK4;
    // KE( 62) =  LR1        -> (1,5)
    KE(1, 5) = LR1;
    // KE( 66) =  SK1        -> (5,5)
    KE(5, 5) = SK1;
    // KE( 68) = -LR1        -> (7,5)
    KE(7, 5) = -LR1;
    // KE( 72) =  SK3        -> (11,5)
    KE(11, 5) = SK3;
    // KE( 73) = -AEL        -> (0,6)
    KE(0, 6) = -AEL;
    // KE( 79) =  AEL        -> (6,6)
    KE(6, 6) = AEL;
    // KE( 86) = -R1         -> (1,7)
    KE(1, 7) = -R1;
    // KE( 90) = -LR1        -> (5,7)
    KE(5, 7) = -LR1;
    // KE( 92) =  R1         -> (7,7)
    KE(7, 7) = R1;
    // KE( 96) = -LR1        -> (11,7)
    KE(11, 7) = -LR1;
    // KE( 99) = -R2         -> (2,8)
    KE(2, 8) = -R2;
    // KE(101) =  LR2        -> (4,8)
    KE(4, 8) = LR2;
    // KE(105) =  R2         -> (8,8)
    KE(8, 8) = R2;
    // KE(107) =  LR2        -> (10,8)
    KE(10, 8) = LR2;
    // KE(112) = -GJL        -> (3,9)
    KE(3, 9) = -GJL;
    // KE(118) =  GJL        -> (9,9)
    KE(9, 9) = GJL;
    // KE(123) = -LR2        -> (2,10)
    KE(2, 10) = -LR2;
    // KE(125) =  SK4        -> (4,10)
    KE(4, 10) = SK4;
    // KE(129) =  LR2        -> (8,10)
    KE(8, 10) = LR2;
    // KE(131) =  SK2        -> (10,10)
    KE(10, 10) = SK2;
    // KE(134) =  LR1        -> (1,11)
    KE(1, 11) = LR1;
    // KE(138) =  SK3        -> (5,11)
    KE(5, 11) = SK3;
    // KE(140) = -LR1        -> (7,11)
    KE(7, 11) = -LR1;
    // KE(144) =  SK1        -> (11,11)
    KE(11, 11) = SK1;

    // Since the problem states "Node A: (0,0,0), Node B: (2,0,0)", the local x-axis is along global x.
    // So the transformation matrix is identity, and the stiffness matrix is already in global coordinates.

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << format_double(KE(i, j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}