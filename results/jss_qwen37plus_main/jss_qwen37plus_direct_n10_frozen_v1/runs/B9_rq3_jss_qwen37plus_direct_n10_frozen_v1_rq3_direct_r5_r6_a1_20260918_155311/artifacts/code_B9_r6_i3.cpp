#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Use double precision throughout
using Real = double;

// Helper function to format double in scientific notation with 15 digits
std::string format_double(Real x) {
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

int main() {
    // Test case parameters
    // Grid A: (1, 0, 0), Grid B: (0, 0, 1)
    Eigen::Vector3d GPA = {1.0, 0.0, 0.0};
    Eigen::Vector3d GPB = {0.0, 0.0, 1.0};
    // Reference vector V = (0, 1, 0)
    Eigen::Vector3d SMALLV = {0.0, 1.0, 0.0};
    // Radius of curvature R = 1.0, Angle = 90 degrees
    Real R = 1.0;
    Real BETAR = 90.0; // degrees
    // Material properties
    Real E = 200e9;
    Real nu = 0.3;
    Real G = E / (2.0 * (1.0 + nu));
    // Circular pipe: outer radius = 0.05, wall thickness = 0.005
    Real ro = 0.05;
    Real t = 0.005;
    Real ri = ro - t;
    
    // Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
    
    Real A = M_PI * (ro * ro - ri * ri);
    Real I1 = M_PI / 4.0 * (ro * ro * ro * ro - ri * ri * ri * ri);
    Real I2 = I1; // symmetric circular cross-section
    Real J = 2.0 * I1; // torsional constant for circular section

    // Initialize ECPT array (100 elements, Fortran 1-indexed -> C++ 0-indexed)
    std::vector<Real> ECPT(100, 0.0);

    // Fill ECPT with test case values (Fortran indexing: ECPT(1) is index 0)
    // ECPT(1): IELID = 1
    ECPT[0] = 1.0;
    // ECPT(2), ECPT(3): ISILNO(2) = [1, 2]
    ECPT[1] = 1.0;
    ECPT[2] = 2.0;
    // ECPT(4), ECPT(5), ECPT(6): SMALLV(3) = [0,1,0]
    ECPT[3] = 0.0;
    ECPT[4] = 1.0;
    ECPT[5] = 0.0;
    // ECPT(7): ICSSV = 0 (basic coordinate system)
    ECPT[6] = 0.0;
    // ECPT(8): IMATID = 1
    ECPT[7] = 1.0;
    // ECPT(9): A
    ECPT[8] = A;
    // ECPT(10): I1
    ECPT[9] = I1;
    // ECPT(11): I2
    ECPT[10] = I2;
    // ECPT(12): FJ (torsional constant J)
    ECPT[11] = J;
    // ECPT(13): NSM = 0
    ECPT[12] = 0.0;
    // ECPT(14): FE = 0
    ECPT[13] = 0.0;
    // ECPT(15)-ECPT(22): stress recovery coefficients = 0
    for (int i = 14; i <= 21; ++i) {
        ECPT[i] = 0.0;
    }
    // ECPT(23): K1 = 1.0 (shear factor)
    ECPT[22] = 1.0;
    // ECPT(24): K2 = 1.0
    ECPT[23] = 1.0;
    // ECPT(25): C = 0
    ECPT[24] = 0.0;
    // ECPT(26)-ECPT(28): KX, KY, KZ = 1.0
    ECPT[25] = 1.0;
    ECPT[26] = 1.0;
    ECPT[27] = 1.0;
    // ECPT(29): R = 1.0
    ECPT[28] = R;
    // ECPT(30): BETAR = 90.0
    ECPT[29] = BETAR;
    // ECPT(31): MCSIDA = 0 (basic)
    ECPT[30] = 0.0;
    // ECPT(32)-ECPT(34): GPA(3) = [1,0,0]
    ECPT[31] = GPA(0);
    ECPT[32] = GPA(1);
    ECPT[33] = GPA(2);
    // ECPT(35): MCSIDB = 0 (basic)
    ECPT[34] = 0.0;
    // ECPT(36)-ECPT(38): GPB(3) = [0,0,1]
    ECPT[35] = GPB(0);
    ECPT[36] = GPB(1);
    ECPT[37] = GPB(2);
    // ECPT(39): ELTEMP = 0
    ECPT[38] = 0.0;

    // Common block variables (simplified for this test)
    int NPVT = 1; // pivot node is node 1
    bool HEAT = false;
    bool ABASIC = true;
    bool BBASIC = true;
    int ICSIDA = 0;
    int ICSIDB = 0;
    int JCSIDA = 31; // Fortran index, so C++ index = 30
    int JCSIDB = 35; // Fortran index, so C++ index = 34

    // Set up DP array (16 elements for coordinates, 1-indexed -> indices 0-15)
    std::vector<Real> DP(16, 0.0);

    // DP(1) to DP(3): GPA in basic coordinates
    DP[0] = GPA(0);
    DP[1] = GPA(1);
    DP[2] = GPA(2);
    // DP(4) to DP(6): GPB in basic coordinates
    DP[3] = GPB(0);
    DP[4] = GPB(1);
    DP[5] = GPB(2);

    // DP(7) to DP(9): SMALLV (reference vector)
    DP[6] = SMALLV(0);
    DP[7] = SMALLV(1);
    DP[8] = SMALLV(2);

    // Normalize reference vector
    Real FLD = std::sqrt(DP[6]*DP[6] + DP[7]*DP[7] + DP[8]*DP[8]);
    if (FLD <= 1e-12) {
        std::cerr << "Error: Reference vector magnitude is zero\n";
        return 1;
    }
    DP[6] /= FLD;
    DP[7] /= FLD;
    DP[8] /= FLD;

    // Compute vector from A to B
    Eigen::Vector3d SMALV0;
    SMALV0(0) = DP[3] - DP[0];
    SMALV0(1) = DP[4] - DP[1];
    SMALV0(2) = DP[5] - DP[2];
    Real FLL = SMALV0.norm();
    if (FLL <= 1e-12) {
        std::cerr << "Error: Vector from A to B has zero length\n";
        return 1;
    }
    SMALV0.normalize();

    // Compute VECK = SMALV0 × VECJ (where VECJ = DP[6:9])
    Eigen::Vector3d VECJ(DP[6], DP[7], DP[8]);
    Eigen::Vector3d VECK = SMALV0.cross(VECJ); // Fixed: was "VE CJ" -> "VECJ"
    FLL = VECK.norm();
    if (FLL <= 1e-12) {
        std::cerr << "Error: Cross product VECK is zero\n";
        return 1;
    }
    VECK.normalize();

    // Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI = VECJ.cross(VECK);
    FLL = VECI.norm();
    if (FLL <= 1e-12) {
        std::cerr << "Error: Cross product VECI is zero\n";
        return 1;
    }
    VECI.normalize();

    // Material properties
    Real MATIDC = ECPT[7]; // IMATID
    Real ELTEMP = ECPT[38]; // ELTEMP
    Real E_val = E;
    Real G_val = G;
    Real NU_val = nu;

    // Flexibility correction factors
    Real KX = ECPT[25];
    Real KY = ECPT[26];
    Real KZ = ECPT[27];

    // Clamp KX, KY, KZ to avoid division by zero
    if (KX < 1e-8) KX = 1.0;
    if (KY < 1e-8) KY = 1.0;
    if (KZ < 1e-8) KZ = 1.0;

    Real FI1 = I1 / KZ;
    Real FI2 = I2 / KY;
    Real FJK = J / KX;

    // Area factors for shear
    Real K1 = ECPT[22];
    Real K2 = ECPT[23];
    if (K1 < 1e-8) K1 = 1.0;
    if (K2 < 1e-8) K2 = 1.0;
    if (K1 > 1.0) K1 = 1.0 / K1;
    if (K2 > 1.0) K2 = 1.0 / K2;

    // Constants for stiffness calculation
    Real DCR = 0.017453292; // degrees to radians conversion factor
    Real T = BETAR * DCR; // angle in radians
    Real RA = R / (A * E_val);
    Real RV1 = R / (2.0 * K1 * G_val * A);
    Real RV2 = (K1 / K2) * RV1;
    Real RT = R / (G_val * FJK * 2.0);
    Real RB0 = R / (E_val * FI2 * 2.0);
    Real RB1 = R / (E_val * FI1);
    Real R2 = R * R;

    // Trigonometric constants
    Real ST = std::sin(T);
    Real CT = std::cos(T);
    Real S2T = std::sin(2.0 * T);
    Real C2T = std::cos(2.0 * T);

    // Initialize 6x6 flexibility matrix F
    Eigen::Matrix<Real, 6, 6> F = Eigen::Matrix<Real, 6, 6>::Zero();

    // Axial terms
    F(0,0) += 0.25 * RA * (2.0 * T + S2T);
    F(1,1) += 0.25 * RA * (2.0 * T - S2T);
    F(0,1) += 0.50 * RA * ST * ST;
    F(1,0) = F(0,1);

    // Shear terms
    F(0,0) += 0.5 * RV1 * (2.0 * T - S2T);
    F(1,1) += 0.5 * RV1 * (2.0 * T + S2T);
    F(2,2) += 2.0 * RV2 * T;
    F(0,1) -= RV1 * ST * ST;
    F(1,0) = F(0,1);

    // Torsion terms
    F(2,2) += 0.5 * RT * R2 * (6.0 * T + S2T - 8.0 * ST);
    F(3,3) += 0.5 * RT * (2.0 * T + S2T);
    F(4,4) += 0.5 * RT * (2.0 * T - S2T);
    F(2,3) += RT * R * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) += RT * R * (2.0 - 2.0 * CT - T * ST);
    F(4,2) = F(2,4);
    F(3,4) += 0.5 * RT * (1.0 - C2T);
    F(4,3) = F(3,4);

    // Bending terms
    F(0,0) += 0.25 * RB1 * R2 * (2.0 * T * (2.0 + C2T) - 3.0 * S2T);
    F(1,1) += 0.25 * RB1 * R2 * (2.0 * T * (2.0 - C2T) + 3.0 * S2T - 8.0 * ST);
    F(2,2) += 0.50 * RB0 * R2 * (2.0 * T - S2T);
    F(3,3) += 0.50 * RB0 * (2.0 * T - S2T);
    F(4,4) += 0.50 * RB0 * (2.0 * T + S2T);
    F(5,5) += RB1 * T;
    F(0,1) += 0.25 * RB1 * R2 * (1.0 + 3.0 * C2T + 2.0 * T * S2T - 4.0 * CT);
    F(1,0) = F(0,1);
    F(0,5) -= RB1 * R * (ST - T * CT);
    F(5,0) = F(0,5);
    F(1,5) += RB1 * R * (T * ST + CT - 1.0);
    F(5,1) = F(1,5);
    F(2,3) += RB0 * R * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) -= RB0 * R * T * ST;
    F(4,2) = F(2,4);
    F(3,4) -= 0.50 * RB0 * (1.0 - C2T);
    F(4,3) = F(3,4);

    // Ensure symmetry
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (i != j) {
                F(j,i) = F(i,j);
            }
        }
    }

    // Invert flexibility matrix to get local stiffness matrix DF (6x6)
    Eigen::Matrix<Real, 6, 6> DF;
    try {
        DF = F.inverse();
    } catch (...) {
        std::cerr << "Error: Flexibility matrix is singular\n";
        return 1;
    }

    // Construct the 12x12 local stiffness matrix S
    Eigen::Matrix<Real, 12, 12> S = Eigen::Matrix<Real, 12, 12>::Zero();

    // H matrix (6x6)
    Eigen::Matrix<Real, 6, 6> H = Eigen::Matrix<Real, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i,i) = -1.0;
    }
    H(3,2) = -R * (1.0 - CT);
    H(4,2) = R * ST;
    H(5,0) = R * (1.0 - CT);
    H(5,1) = -R * ST;

    // Fill S matrix as per Fortran logic:
    // S(7:12,7:12) = DF
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i+6, j+6) = DF(i,j);
        }
    }

    // S(1:6,7:12) = H * DF
    Eigen::Matrix<Real, 6, 6> temp1 = H * DF;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j+6) = temp1(i,j);
        }
    }

    // S(1:6,1:6) = (H * DF) * H^T = temp1 * H.transpose()
    Eigen::Matrix<Real, 6, 6> temp2 = temp1 * H.transpose();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j) = temp2(i,j);
        }
    }

    // Enforce symmetry: S = (S + S^T) / 2
    S = 0.5 * (S + S.transpose());

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << format_double(S(i,j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}