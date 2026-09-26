#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function to format double in scientific notation with fixed precision
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos) {
        if (s[pos] == '.') {
            s.erase(pos); // remove '.' only if it's the last character
        } else {
            s.erase(pos + 1);
        }
    }
    // Remove '+' from exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos && epos + 1 < s.length() && s[epos + 1] == '+') {
        s.erase(epos + 1, 1);
    }
    return s;
}

// Invert a 6x6 matrix using Eigen's partialPivLu
Eigen::Matrix<double, 6, 6> invert_6x6(const Eigen::Matrix<double, 6, 6>& A) {
    Eigen::Matrix<double, 6, 6> A_inv;
    Eigen::PartialPivLU<Eigen::Matrix<double, 6, 6>> lu(A);
    A_inv = lu.inverse();
    return A_inv;
}

int main() {
    // Test case parameters
    const double PI = M_PI;
    const double DCR = PI / 180.0; // degrees to radians conversion factor

    // Grid points
    Eigen::Vector3d GA(1.0, 0.0, 0.0); // Grid A
    Eigen::Vector3d GB(0.0, 0.0, 1.0); // Grid B
    Eigen::Vector3d SMALLV(0.0, 1.0, 0.0); // Reference vector

    // Geometry
    double R = 1.0; // radius of curvature
    double BETAR = 90.0; // angle in degrees

    // Material properties
    double E = 200e9; // Pa
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));

    // Cross-section: circular pipe
    double ro = 0.05; // outer radius
    double t = 0.005; // wall thickness
    double ri = ro - t; // inner radius
    double A = PI * (ro * ro - ri * ri); // area
    double I1 = PI / 4.0 * (ro * ro * ro * ro - ri * ri * ri * ri); // moment of inertia about y-axis
    double I2 = I1; // symmetric, same about z-axis
    double J = 2.0 * I1; // torsional constant

    // Initialize ECPT array (100 elements, double precision)
    std::vector<double> ECPT(100, 0.0);

    // Fill ECPT with test case values
    // ECPT(1) = IELID (element ID) - unused in computation but needed for indexing
    ECPT[0] = 1.0;

    // ECPT(2), ECPT(3) = ISILNO(1), ISILNO(2) - grid point IDs
    ECPT[1] = 1.0; // grid A ID
    ECPT[2] = 2.0; // grid B ID

    // ECPT(4), ECPT(5), ECPT(6) = SMALLV(3) - reference vector
    ECPT[3] = SMALLV(0);
    ECPT[4] = SMALLV(1);
    ECPT[5] = SMALLV(2);

    // ECPT(7) = ICSSV - coordinate system ID for SMALLV (0 = basic)
    ECPT[6] = 0.0;

    // ECPT(8) = IMATID - material ID
    ECPT[7] = 1.0;

    // ECPT(9) = A - cross-sectional area
    ECPT[8] = A;

    // ECPT(10), ECPT(11) = I1, I2 - area moments of inertia
    ECPT[9] = I1;
    ECPT[10] = I2;

    // ECPT(12) = FJ - torsional constant
    ECPT[11] = J;

    // ECPT(13) = NSM - non-structural mass (not used)
    ECPT[12] = 0.0;

    // ECPT(14) = FE - force element description (not used)
    ECPT[13] = 0.0;

    // ECPT(15)-ECPT(22) = stress recovery coefficients (not used)
    // ECPT(23), ECPT(24) = K1, K2 - shear area factors
    ECPT[22] = 1.0; // K1
    ECPT[23] = 1.0; // K2

    // ECPT(25) = C - stress intensification factor (not used)
    ECPT[24] = 0.0;

    // ECPT(26)-ECPT(28) = KX, KY, KZ - flexibility correction factors
    ECPT[25] = 1.0; // KX
    ECPT[26] = 1.0; // KY
    ECPT[27] = 1.0; // KZ

    // ECPT(29) = R - radius of curvature
    ECPT[28] = R;

    // ECPT(30) = BETAR - angle from GA to GB in degrees
    ECPT[29] = BETAR;

    // ECPT(31) = MCSIDA - coordinate system ID for grid A (0 = basic)
    ECPT[30] = 0.0;

    // ECPT(32)-ECPT(34) = GPA(3) - basic coordinates for grid A
    ECPT[31] = GA(0);
    ECPT[32] = GA(1);
    ECPT[33] = GA(2);

    // ECPT(35) = MCSIDB - coordinate system ID for grid B (0 = basic)
    ECPT[34] = 0.0;

    // ECPT(36)-ECPT(38) = GPB(3) - basic coordinates for grid B
    ECPT[35] = GB(0);
    ECPT[36] = GB(1);
    ECPT[37] = GB(2);

    // ECPT(39) = ELTEMP - average element temperature (not used)
    ECPT[38] = 0.0;

    // NPVT is the pivot node - we'll compute for both nodes, but per NASTRAN logic,
    // we need to determine which is the pivot. For this test, we'll assume NPVT = 1 (grid A)
    int NPVT = 1;

    // Set up common block variables (simplified for this test)
    int ISILNO[2] = {1, 2}; // grid point IDs
    double A_val = A;
    double I1_val = I1;
    double I2_val = I2;
    double FJ_val = J;
    double K1_val = 1.0;
    double K2_val = 1.0;
    double KX_val = 1.0;
    double KY_val = 1.0;
    double KZ_val = 1.0;
    double R_val = R;
    double BETAR_val = BETAR;
    double MCSIDA_val = 0.0;
    double MCSIDB_val = 0.0;
    double GPA[3] = {GA(0), GA(1), GA(2)};
    double GPB[3] = {GB(0), GB(1), GB(2)};
    double TEMPEL_val = 0.0;
    int IMATID_val = 1;

    // Determine pivot point
    int IPVT = 1;
    if (ISILNO[0] != NPVT) {
        IPVT = 2;
        if (ISILNO[1] != NPVT) {
            // Error case - should not happen in test
            std::cerr << "Error: NPVT not found in ISILNO" << std::endl;
            return 1;
        }
    }

    // Set up pointers to coordinate system IDs
    int JCSIDA = 31; // index 30 in 0-based
    int JCSIDB = 35; // index 34 in 0-based

    // Bounds check
    if (JCSIDA - 1 >= static_cast<int>(ECPT.size()) || JCSIDB - 1 >= static_cast<int>(ECPT.size())) {
        std::cerr << "Error: ECPT index out of bounds" << std::endl;
        return 1;
    }

    int ICSIDA = static_cast<int>(ECPT[JCSIDA - 1]); // ECPT(31) -> index 30
    int ICSIDB = static_cast<int>(ECPT[JCSIDB - 1]); // ECPT(35) -> index 34

    // Define location of end A, end B in terms of DP(1) through DP(6)
    std::vector<double> DP(16, 0.0);

    // Bounds checks for DP assignments
    auto safe_get = [&](size_t idx) -> double {
        return (idx < ECPT.size()) ? ECPT[idx] : 0.0;
    };

    DP[0] = safe_get(JCSIDA);     // DP(1) = ECPT(JCSIDA+1) -> index JCSIDA (31) -> ECPT[31]
    DP[1] = safe_get(JCSIDA + 1); // DP(2) = ECPT(JCSIDA+2) -> ECPT[32]
    DP[2] = safe_get(JCSIDA + 2); // DP(3) = ECPT(JCSIDA+3) -> ECPT[33]
    DP[3] = safe_get(JCSIDB);     // DP(4) = ECPT(JCSIDB+1) -> ECPT[34]
    DP[4] = safe_get(JCSIDB + 1); // DP(5) = ECPT(JCSIDB+2) -> ECPT[35]
    DP[5] = safe_get(JCSIDB + 2); // DP(6) = ECPT(JCSIDB+3) -> ECPT[36]

    // Define components of vector from end A to center of curvature
    // In NASTRAN, SMALLV(3) contains the vector from A to center of curvature
    DP[6] = safe_get(3); // SMALLV(1) -> ECPT[3] (index 3)
    DP[7] = safe_get(4); // SMALLV(2) -> ECPT[4] (index 4)
    DP[8] = safe_get(5); // SMALLV(3) -> ECPT[5] (index 5)

    double FLD = std::sqrt(DP[6] * DP[6] + DP[7] * DP[7] + DP[8] * DP[8]);
    if (FLD <= 0.0) {
        std::cerr << "Error: Zero length vector to center of curvature" << std::endl;
        return 1;
    }
    DP[6] /= FLD;
    DP[7] /= FLD;
    DP[8] /= FLD;

    // Determine if points are in basic coordinates
    bool ABASIC = (ICSIDA == 0);
    bool BBASIC = (ICSIDB == 0);

    // Compute transformation matrices - for basic coordinates, identity
    Eigen::Matrix3d TA = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d TB = Eigen::Matrix3d::Identity();

    // VECJ = unit vector from A to center of curvature
    Eigen::Vector3d VECJ(DP[6], DP[7], DP[8]);

    // Transform coordinates if needed (not needed for basic)
    if (!ABASIC) {
        // Not implemented - assume basic
    }
    if (!BBASIC) {
        // Not implemented - assume basic
    }

    // Calculate true length of elbow
    double T = BETAR_val * DCR; // angle in radians
    double FL = R_val * T;
    if (FL == 0.0) {
        std::cerr << "Error: Zero elbow length" << std::endl;
        return 1;
    }

    // Construct vector from A to B
    Eigen::Vector3d SMALV0;
    SMALV0 << DP[3] - DP[0], DP[4] - DP[1], DP[5] - DP[2];
    double FLL = SMALV0.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero vector from A to B" << std::endl;
        return 1;
    }
    SMALV0.normalize();

    // Compute VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK = SMALV0.cross(VECJ);
    FLL = VECK.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero cross product for VECK" << std::endl;
        return 1;
    }
    VECK.normalize();

    // Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI = VECJ.cross(VECK);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero cross product for VECI" << std::endl;
        return 1;
    }
    VECI.normalize();

    // Material properties - use given E, G, nu
    double E_mat = E;
    double G_mat = G;
    double nu_mat = nu;

    // Set up intermediate variables
    double KX = KX_val;
    double KY = KY_val;
    double KZ = KZ_val;
    if (KX < 1.0e-8) KX = 1.0;
    if (KY < 1.0e-8) KY = 1.0;
    if (KZ < 1.0e-8) KZ = 1.0;
    double FI1 = I1_val / KZ;
    double FI2 = I2_val / KY;
    double FJK = FJ_val / KX;

    // Area factors for shear
    double K1 = K1_val;
    double K2 = K2_val;
    if (K1 < 1.0e-8) K1 = 1.0;
    if (K2 < 1.0e-8) K2 = 1.0;
    if (K1 > 1.0) K1 = 1.0 / K1;
    if (K2 > 1.0) K2 = 1.0 / K2;

    // Section property constants
    double RA = R_val / (A_val * E_mat);
    double RV1 = R_val / (2.0 * K1 * G_mat * A_val);
    double RV2 = (K1 / K2) * RV1;
    double RT = R_val / (G_mat * FJK * 2.0);
    double RB0 = R_val / (E_mat * FI2 * 2.0);
    double RB1 = R_val / (E_mat * FI1);
    double R2 = R_val * R_val;

    // Trigonometric constants
    double ST = std::sin(T);
    double CT = std::cos(T);
    double S2T = std::sin(2.0 * T);
    double C2T = std::cos(2.0 * T);

    // Form the 6x6 flexibility matrix F
    Eigen::Matrix<double, 6, 6> F = Eigen::Matrix<double, 6, 6>::Zero();

    // Axial
    F(0,0) += 0.25 * RA * (2.0 * T + S2T);
    F(1,1) += 0.25 * RA * (2.0 * T - S2T);
    F(0,1) += 0.50 * RA * ST * ST;
    F(1,0) = F(0,1);

    // Shear
    F(0,0) += 0.5 * RV1 * (2.0 * T - S2T);
    F(1,1) += 0.5 * RV1 * (2.0 * T + S2T);
    F(2,2) += 2.0 * RV2 * T;
    F(0,1) -= RV1 * ST * ST;
    F(1,0) = F(0,1);

    // Torsion
    F(2,2) += 0.5 * RT * R2 * (6.0 * T + S2T - 8.0 * ST);
    F(3,3) += 0.5 * RT * (2.0 * T + S2T);
    F(4,4) += 0.5 * RT * (2.0 * T - S2T);
    F(2,3) += RT * R_val * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) += RT * R_val * (2.0 - 2.0 * CT - T * ST);
    F(4,2) = F(2,4);
    F(3,4) += 0.5 * RT * (1.0 - C2T);
    F(4,3) = F(3,4);

    // Bending
    F(0,0) += 0.25 * RB1 * R2 * (2.0 * T * (2.0 + C2T) - 3.0 * S2T);
    F(1,1) += 0.25 * RB1 * R2 * (2.0 * T * (2.0 - C2T) + 3.0 * S2T - 8.0 * ST);
    F(2,2) += 0.50 * RB0 * R2 * (2.0 * T - S2T);
    F(3,3) += 0.50 * RB0 * (2.0 * T - S2T);
    F(4,4) += 0.50 * RB0 * (2.0 * T + S2T);
    F(5,5) += RB1 * T;
    F(0,1) += 0.25 * RB1 * R2 * (1.0 + 3.0 * C2T + 2.0 * T * S2T - 4.0 * CT);
    F(1,0) = F(0,1);
    F(0,5) -= RB1 * R_val * (ST - T * CT);
    F(5,0) = F(0,5);
    F(1,5) += RB1 * R_val * (T * ST + CT - 1.0);
    F(5,1) = F(1,5);
    F(2,3) += RB0 * R_val * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) -= RB0 * R_val * T * ST;
    F(4,2) = F(2,4);
    F(3,4) -= 0.50 * RB0 * (1.0 - C2T);
    F(4,3) = F(3,4);

    // Ensure symmetry
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (i != j && std::abs(F(i,j) - F(j,i)) > 1e-15) {
                F(j,i) = F(i,j);
            }
        }
    }

    // Invert flexibility matrix to get local stiffness matrix DF
    Eigen::Matrix<double, 6, 6> DF = invert_6x6(F);

    // Set up H matrix (force transformation)
    Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i,i) = -1.0;
    }
    H(3,2) = -R_val * (1.0 - CT); // H(4,3) in Fortran (1-indexed)
    H(4,2) = R_val * ST;          // H(5,3) in Fortran
    H(5,0) = -H(3,2);            // H(6,1) in Fortran
    H(5,1) = -H(4,2);            // H(6,2) in Fortran

    // Form the 12x12 stiffness matrix S
    Eigen::Matrix<double, 12, 12> S = Eigen::Matrix<double, 12, 12>::Zero();

    // S(7:12,7:12) = DF
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i+6, j+6) = DF(i,j);
        }
    }

    // S(1:6,7:12) = H * DF
    Eigen::Matrix<double, 6, 6> temp1 = H * DF;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j+6) = temp1(i,j);
        }
    }

    // S(1:6,1:6) = (H * DF) * H^T
    Eigen::Matrix<double, 6, 6> temp2 = temp1 * H.transpose();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j) = temp2(i,j);
        }
    }

    // Reflect for symmetry: S = S^T
    for (int i = 0; i < 12; ++i) {
        for (int j = i+1; j < 12; ++j) {
            S(j,i) = S(i,j);
        }
    }

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << format_double(S(i,j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}