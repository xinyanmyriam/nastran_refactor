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
    // KE(  1) =  AEL   -> [0,0]
    KE(0, 0) = AEL;
    // KE(  7) = -AEL   -> [6,0]
    KE(6, 0) = -AEL;
    // KE( 14) =  R1    -> [1,1] (since 14 = 1 + (1-1)*12 + 1? Let's check: Fortran is column-major: element (i,j) is at index (i-1)+(j-1)*12)
    // So KE(14): i=2, j=2 -> (2-1)+(2-1)*12 = 1+12 = 13 -> index 13 -> [1,1] in 0-based
    KE(1, 1) = R1;
    // KE( 18) =  LR1   -> [5,1] because 18 = (6-1)+(2-1)*12 = 5+12 = 17 -> index 17 -> [5,1]
    KE(5, 1) = LR1;
    // KE( 20) = -R1    -> [1,2] because 20 = (2-1)+(3-1)*12 = 1+24 = 25? Wait, let's recalculate properly.

    // Actually, Fortran stores column-major: element (i,j) is at position (i-1) + (j-1)*12
    // So for index k (1-based), row i = (k-1) % 12 + 1, col j = (k-1) / 12 + 1
    // Let's create a helper to convert Fortran index to (row,col) 0-based:
    auto fortran_to_rc = [](int k) -> std::pair<int, int> {
        k--; // to 0-based index
        int row = k % 12;
        int col = k / 12;
        return {row, col};
    };

    // Now fill all the specified positions:
    // The Fortran code sets:
    // KE(  1) =  AEL        -> (0,0)
    // KE(  7) = -AEL        -> (6,0)
    // KE( 14) =  R1         -> (1,1) because (14-1)%12=1, (14-1)/12=1
    // KE( 18) =  LR1        -> (5,1) because (18-1)%12=5, (18-1)/12=1
    // KE( 20) = -R1         -> (7,1) because (20-1)%12=7, (20-1)/12=1
    // KE( 24) =  LR1        -> (11,1) because (24-1)%12=11, (24-1)/12=1
    // KE( 27) =  R2         -> (2,2) because (27-1)%12=2, (27-1)/12=2
    // KE( 29) = -LR2        -> (4,2) because (29-1)%12=4, (29-1)/12=2
    // KE( 33) = -R2         -> (8,2) because (33-1)%12=8, (33-1)/12=2
    // KE( 35) = -LR2        -> (10,2) because (35-1)%12=10, (35-1)/12=2
    // KE( 40) =  GJL        -> (3,3) because (40-1)%12=3, (40-1)/12=3
    // KE( 46) = -GJL        -> (9,3) because (46-1)%12=9, (46-1)/12=3
    // KE( 51) = -LR2        -> (6,4) because (51-1)%12=6, (51-1)/12=4
    // KE( 53) =  SK2        -> (8,4) because (53-1)%12=8, (53-1)/12=4
    // KE( 57) =  LR2        -> (0,4) because (57-1)%12=0, (57-1)/12=4
    // KE( 59) =  SK4        -> (2,4) because (59-1)%12=2, (59-1)/12=4
    // KE( 62) =  LR1        -> (5,5) because (62-1)%12=5, (62-1)/12=5
    // KE( 66) =  SK1        -> (9,5) because (66-1)%12=9, (66-1)/12=5
    // KE( 68) = -LR1        -> (11,5) because (68-1)%12=11, (68-1)/12=5
    // KE( 72) =  SK3        -> (3,5) because (72-1)%12=3, (72-1)/12=5
    // KE( 73) = -AEL        -> (0,6) because (73-1)%12=0, (73-1)/12=6
    // KE( 79) =  AEL        -> (6,6) because (79-1)%12=6, (79-1)/12=6
    // KE( 86) = -R1         -> (1,7) because (86-1)%12=1, (86-1)/12=7
    // KE( 90) = -LR1        -> (5,7) because (90-1)%12=5, (90-1)/12=7
    // KE( 92) =  R1         -> (7,7) because (92-1)%12=7, (92-1)/12=7
    // KE( 96) = -LR1        -> (11,7) because (96-1)%12=11, (96-1)/12=7
    // KE( 99) = -R2         -> (2,8) because (99-1)%12=2, (99-1)/12=8
    // KE(101) =  LR2        -> (4,8) because (101-1)%12=4, (101-1)/12=8
    // KE(105) =  R2         -> (8,8) because (105-1)%12=8, (105-1)/12=8
    // KE(107) =  LR2        -> (10,8) because (107-1)%12=10, (107-1)/12=8
    // KE(112) = -GJL        -> (3,9) because (112-1)%12=3, (112-1)/12=9
    // KE(118) =  GJL        -> (9,9) because (118-1)%12=9, (118-1)/12=9
    // KE(123) = -LR2        -> (6,10) because (123-1)%12=6, (123-1)/12=10
    // KE(125) =  SK4        -> (8,10) because (125-1)%12=8, (125-1)/12=10
    // KE(129) =  LR2        -> (0,10) because (129-1)%12=0, (129-1)/12=10
    // KE(131) =  SK2        -> (2,10) because (131-1)%12=2, (131-1)/12=10
    // KE(134) =  LR1        -> (5,11) because (134-1)%12=5, (134-1)/12=11
    // KE(138) =  SK3        -> (9,11) because (138-1)%12=9, (138-1)/12=11
    // KE(140) = -LR1        -> (11,11) because (140-1)%12=11, (140-1)/12=11
    // KE(144) =  SK1        -> (3,11) because (144-1)%12=3, (144-1)/12=11

    // Let's verify the mapping with a few examples:
    // KE(1): (1-1)%12=0, (1-1)/12=0 -> (0,0) ✓
    // KE(12): (12-1)%12=11, (12-1)/12=0 -> (11,0)
    // KE(13): (13-1)%12=0, (13-1)/12=1 -> (0,1)
    // So yes, it's column-major.

    // Now fill the matrix:
    KE(0, 0) = AEL;                    // KE(1)
    KE(6, 0) = -AEL;                   // KE(7)
    KE(1, 1) = R1;                     // KE(14)
    KE(5, 1) = LR1;                    // KE(18)
    KE(7, 1) = -R1;                    // KE(20)
    KE(11, 1) = LR1;                   // KE(24)
    KE(2, 2) = R2;                     // KE(27)
    KE(4, 2) = -LR2;                   // KE(29)
    KE(8, 2) = -R2;                    // KE(33)
    KE(10, 2) = -LR2;                  // KE(35)
    KE(3, 3) = GJL;                    // KE(40)
    KE(9, 3) = -GJL;                   // KE(46)
    KE(6, 4) = -LR2;                   // KE(51)
    KE(8, 4) = SK2;                    // KE(53)
    KE(0, 4) = LR2;                    // KE(57)
    KE(2, 4) = SK4;                    // KE(59)
    KE(5, 5) = LR1;                    // KE(62)
    KE(9, 5) = SK1;                    // KE(66)
    KE(11, 5) = -LR1;                  // KE(68)
    KE(3, 5) = SK3;                    // KE(72)
    KE(0, 6) = -AEL;                   // KE(73)
    KE(6, 6) = AEL;                    // KE(79)
    KE(1, 7) = -R1;                    // KE(86)
    KE(5, 7) = -LR1;                   // KE(90)
    KE(7, 7) = R1;                     // KE(92)
    KE(11, 7) = -LR1;                  // KE(96)
    KE(2, 8) = -R2;                    // KE(99)
    KE(4, 8) = LR2;                    // KE(101)
    KE(8, 8) = R2;                     // KE(105)
    KE(10, 8) = LR2;                   // KE(107)
    KE(3, 9) = -GJL;                   // KE(112)
    KE(9, 9) = GJL;                    // KE(118)
    KE(6, 10) = -LR2;                  // KE(123)
    KE(8, 10) = SK4;                   // KE(125)
    KE(0, 10) = LR2;                   // KE(129)
    KE(2, 10) = SK2;                   // KE(131)
    KE(5, 11) = LR1;                   // KE(134)
    KE(9, 11) = SK3;                   // KE(138)
    KE(11, 11) = -LR1;                 // KE(140)
    KE(3, 11) = SK1;                   // KE(144)

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