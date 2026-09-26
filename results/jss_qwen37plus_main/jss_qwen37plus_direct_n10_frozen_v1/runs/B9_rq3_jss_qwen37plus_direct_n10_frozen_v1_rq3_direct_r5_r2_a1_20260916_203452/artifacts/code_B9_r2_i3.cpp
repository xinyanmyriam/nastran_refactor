#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <math.h>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Fallback in case M_PI is still not defined (e.g., some strict conforming modes)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function to format double in scientific notation with fixed precision
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos + 1);
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
    const double PI = M_PI;
    const double DCR = PI / 180.0; // degrees to radians conversion factor

    // Grid points
    Eigen::Vector3d GPA(1.0, 0.0, 0.0); // Grid A
    Eigen::Vector3d GPB(0.0, 0.0, 1.0); // Grid B
    Eigen::Vector3d SMALLV(0.0, 1.0, 0.0); // Reference vector V

    // Geometry
    double R = 1.0; // Radius of curvature
    double BETAR = 90.0; // Angle in degrees

    // Material properties
    double E = 200e9; // Pa
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));

    // Cross-section (circular pipe)
    double ro = 0.05; // outer radius
    double t = 0.005; // wall thickness
    double ri = ro - t;
    double A = PI * (ro * ro - ri * ri);
    double I1 = PI / 4.0 * (ro * ro * ro * ro - ri * ri * ri * ri); // Iy
    double I2 = I1; // Iz (symmetric)
    double FJ = 2.0 * I1; // J = polar moment

    // Initialize ECPT array (100 elements, Fortran 1-indexed -> C++ 0-indexed)
    std::vector<double> ECPT(100, 0.0);

    // Fill ECPT with test case values (using Fortran indexing convention: ECPT(1) is index 0)
    ECPT[0]  = 1.0; // IELID
    ECPT[1]  = 1.0; // ISILNO(1) - grid A ID
    ECPT[2]  = 2.0; // ISILNO(2) - grid B ID
    ECPT[3]  = SMALLV(0); // SMALLV(1)
    ECPT[4]  = SMALLV(1); // SMALLV(2)
    ECPT[5]  = SMALLV(2); // SMALLV(3)
    ECPT[6]  = 0.0; // ICSSV
    ECPT[7]  = 1.0; // IMATID
    ECPT[8]  = A;   // A
    ECPT[9]  = I1;  // I1
    ECPT[10] = I2;  // I2
    ECPT[11] = FJ;  // FJ
    ECPT[12] = 0.0; // NSM
    ECPT[13] = 0.0; // FE
    ECPT[14] = 0.0; // R1
    ECPT[15] = 0.0; // T1
    ECPT[16] = 0.0; // R2
    ECPT[17] = 0.0; // T2
    ECPT[18] = 0.0; // R3
    ECPT[19] = 0.0; // T3
    ECPT[20] = 0.0; // R4
    ECPT[21] = 0.0; // T4
    ECPT[22] = 1.0; // K1 (shear factor)
    ECPT[23] = 1.0; // K2 (shear factor)
    ECPT[24] = 0.0; // C
    ECPT[25] = 1.0; // KX
    ECPT[26] = 1.0; // KY
    ECPT[27] = 1.0; // KZ
    ECPT[28] = R;   // R
    ECPT[29] = BETAR; // BETAR
    ECPT[30] = 0.0; // MCSIDA
    ECPT[31] = GPA(0); // GPA(1)
    ECPT[32] = GPA(1); // GPA(2)
    ECPT[33] = GPA(2); // GPA(3)
    ECPT[34] = 0.0; // MCSIDB
    ECPT[35] = GPB(0); // GPB(1)
    ECPT[36] = GPB(1); // GPB(2)
    ECPT[37] = GPB(2); // GPB(3)
    ECPT[38] = 0.0; // ELTEMP

    // Simulate common blocks and variables
    int NPVT = 1; // pivot node is node 1 (grid A)
    bool HEAT = false;

    // Extract values as in Fortran
    int IPVT = 1;
    if (static_cast<int>(ECPT[1]) != NPVT) {
        IPVT = 2;
    }

    // Set up coordinate system IDs
    int JCSIDA = 31; // Fortran index 31 -> C++ index 30
    int JCSIDB = 35; // Fortran index 35 -> C++ index 34
    int ICSIDA = static_cast<int>(ECPT[30]); // ECPT(31)
    int ICSIDB = static_cast<int>(ECPT[34]); // ECPT(35)

    // DP array: store coordinates and vectors
    std::vector<double> DP(17, 0.0); // DP(1) to DP(16), 1-indexed -> 0-indexed

    // Define location of end A, end B
    DP[0] = ECPT[JCSIDA];     // DP(1) = ECPT(JCSIDA+1) -> ECPT[31]
    DP[1] = ECPT[JCSIDA+1];   // DP(2) = ECPT(JCSIDA+2) -> ECPT[32]
    DP[2] = ECPT[JCSIDA+2];   // DP(3) = ECPT(JCSIDA+3) -> ECPT[33]
    DP[3] = ECPT[JCSIDB];     // DP(4) = ECPT(JCSIDB+1) -> ECPT[35]
    DP[4] = ECPT[JCSIDB+1];   // DP(5) = ECPT(JCSIDB+2) -> ECPT[36]
    DP[5] = ECPT[JCSIDB+2];   // DP(6) = ECPT(JCSIDB+3) -> ECPT[37]

    // Vector from end A to center of curvature (SMALLV)
    DP[6] = ECPT[3]; // SMALLV(1)
    DP[7] = ECPT[4]; // SMALLV(2)
    DP[8] = ECPT[5]; // SMALLV(3)
    double FLD = std::sqrt(DP[6]*DP[6] + DP[7]*DP[7] + DP[8]*DP[8]);
    if (FLD <= 0.0) {
        std::cerr << "Error: Zero length vector for center of curvature\n";
        return 1;
    }
    DP[6] /= FLD;
    DP[7] /= FLD;
    DP[8] /= FLD;

    // Determine if points are in basic coordinates
    bool ABASIC = (ICSIDA == 0);
    bool BBASIC = (ICSIDB == 0);

    // VECJ = normalized SMALLV (center direction)
    Eigen::Vector3d VECJ(DP[6], DP[7], DP[8]);

    // Compute true length of elbow
    double T = BETAR * DCR; // angle in radians
    double FL = R * T;
    if (FL == 0.0) {
        std::cerr << "Error: Zero elbow length\n";
        return 1;
    }

    // Vector from A to B
    Eigen::Vector3d SMALV0(DP[3]-DP[0], DP[4]-DP[1], DP[5]-DP[2]);
    double FLL = SMALV0.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero vector from A to B\n";
        return 1;
    }
    SMALV0.normalize();

    // Compute VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK = SMALV0.cross(VECJ);
    FLL = VECK.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero cross product for VECK\n";
        return 1;
    }
    VECK.normalize();

    // Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI = VECJ.cross(VECK);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero cross product for VECI\n";
        return 1;
    }
    VECI.normalize();

    // Material properties
    double E_val = E;
    double G_val = G;
    double nu_val = nu;

    // Flexibility correction factors
    double KX = ECPT[24];
    double KY = ECPT[25];
    double KZ = ECPT[26];
    if (KX < 1.0e-8) KX = 1.0;
    if (KY < 1.0e-8) KY = 1.0;
    if (KZ < 1.0e-8) KZ = 1.0;
    double FI1 = I1 / KZ;
    double FI2 = I2 / KY;
    double FJK = FJ / KX;

    // Shear area factors (Fortran logic: if >1.0, invert)
    double K1 = ECPT[21];
    double K2 = ECPT[22];
    if (K1 < 1.0e-8) K1 = 1.0;
    if (K2 < 1.0e-8) K2 = 1.0;
    if (K1 > 1.0) K1 = 1.0 / K1;
    if (K2 > 1.0) K2 = 1.0 / K2;

    // Section property constants
    double RA  = R / (A * E_val);
    double RV1 = R / (2.0 * K1 * G_val * A);
    double RV2 = (K1 / K2) * RV1;
    double RT  = R / (G_val * FJK * 2.0);
    double RB0 = R / (E_val * FI2 * 2.0);
    double RB1 = R / (E_val * FI1);
    double R2  = R * R;

    // Trigonometric constants
    double ST  = std::sin(T);
    double CT  = std::cos(T);
    double S2T = std::sin(2.0 * T);
    double C2T = std::cos(2.0 * T);

    // Build 6x6 flexibility matrix F (Fortran F(6,6))
    Eigen::Matrix<double, 6, 6> F = Eigen::Matrix<double, 6, 6>::Zero();

    // Axial
    F(0,0) += 0.25 * RA * (2.0*T + S2T);
    F(1,1) += 0.25 * RA * (2.0*T - S2T);
    F(0,1) += 0.50 * RA * ST * ST;
    F(1,0) = F(0,1);

    // Shear
    F(0,0) += 0.5 * RV1 * (2.0*T - S2T);
    F(1,1) += 0.5 * RV1 * (2.0*T + S2T);
    F(2,2) += 2.0 * RV2 * T;
    F(0,1) -= RV1 * ST * ST;
    F(1,0) = F(0,1);

    // Torsion
    F(2,2) += 0.5 * RT * R2 * (6.0*T + S2T - 8.0*ST);
    F(3,3) += 0.5 * RT * (2.0*T + S2T);
    F(4,4) += 0.5 * RT * (2.0*T - S2T);
    F(2,3) += RT * R * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) += RT * R * (2.0 - 2.0*CT - T*ST);
    F(4,2) = F(2,4);
    F(3,4) += 0.5 * RT * (1.0 - C2T);
    F(4,3) = F(3,4);

    // Bending
    F(0,0) += 0.25 * RB1 * R2 * (2.0*T*(2.0+C2T) - 3.0*S2T);
    F(1,1) += 0.25 * RB1 * R2 * (2.0*T*(2.0-C2T) + 3.0*S2T - 8.0*ST);
    F(2,2) += 0.50 * RB0 * R2 * (2.0*T - S2T);
    F(3,3) += 0.50 * RB0 * (2.0*T - S2T);
    F(4,4) += 0.50 * RB0 * (2.0*T + S2T);
    F(5,5) += RB1 * T;
    F(0,1) += 0.25 * RB1 * R2 * (1.0 + 3.0*C2T + 2.0*T*S2T - 4.0*CT);
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
            if (i != j && std::abs(F(i,j) - F(j,i)) > 1e-15) {
                F(i,j) = 0.5 * (F(i,j) + F(j,i));
                F(j,i) = F(i,j);
            }
        }
    }

    // Invert flexibility matrix to get local stiffness submatrix DF
    // Use Eigen's partialPivLu for robust inversion
    Eigen::Matrix<double, 6, 6> DF;
    try {
        DF = F.inverse();
    } catch (...) {
        std::cerr << "Error: Failed to invert flexibility matrix\n";
        return 1;
    }

    // Build 12x12 stiffness matrix S
    Eigen::Matrix<double, 12, 12> S = Eigen::Matrix<double, 12, 12>::Zero();

    // Set bottom-right 6x6 block: S(7:12,7:12) = DF
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i+6, j+6) = DF(i, j);
        }
    }

    // Build H matrix (6x6)
    Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i, i) = -1.0;
    }
    H(3,2) = -R * (1.0 - CT); // H(4,3) in Fortran (1-indexed)
    H(4,2) = R * ST;          // H(5,3) in Fortran
    H(5,0) = -H(3,2);         // H(6,1) in Fortran
    H(5,1) = -H(4,2);         // H(6,2) in Fortran

    // Compute S(1:6,7:12) = H * DF
    Eigen::Matrix<double, 6, 6> temp1 = H * DF;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j+6) = temp1(i, j);
        }
    }

    // Compute S(1:6,1:6) = (H * DF) * H^T = temp1 * H.transpose()
    Eigen::Matrix<double, 6, 6> temp2 = temp1 * H.transpose();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j) = temp2(i, j);
        }
    }

    // Enforce symmetry: S = (S + S^T) / 2
    Eigen::Matrix<double, 12, 12> S_sym = 0.5 * (S + S.transpose());

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << format_double(S_sym(i, j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}