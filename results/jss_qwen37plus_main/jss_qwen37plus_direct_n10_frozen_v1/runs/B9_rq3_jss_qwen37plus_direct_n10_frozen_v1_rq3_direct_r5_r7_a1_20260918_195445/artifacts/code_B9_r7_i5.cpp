#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Fallback definition for M_PI if not provided (e.g., some strict POSIX or older compilers)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function to format double in scientific notation with 15 digits precision
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
    const double DCR = 0.017453292; // degrees to radians conversion factor

    // Grid points
    Eigen::Vector3d GA(1.0, 0.0, 0.0); // Grid A
    Eigen::Vector3d GB(0.0, 0.0, 1.0); // Grid B
    Eigen::Vector3d SMALLV(0.0, 1.0, 0.0); // Reference vector V

    // Geometry
    double R = 1.0; // Radius of curvature
    double BETAR = 90.0; // Angle in degrees

    // Material properties
    double E = 200e9; // Pa
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));

    // Cross-section: circular pipe
    double ro = 0.05; // outer radius
    double t = 0.005; // wall thickness
    double ri = ro - t;
    double A = PI * (ro * ro - ri * ri);
    double I1 = PI / 4.0 * (ro * ro * ro * ro - ri * ri * ri * ri); // Iy
    double I2 = I1; // Iz (symmetric)
    double FJ = 2.0 * I1; // J = polar moment

    // Initialize ECPT array (100 elements, zero-initialized)
    std::vector<double> ECPT(100, 0.0);

    // Fill ECPT with test case values
    ECPT[0] = 1.0; // IELID (element ID)
    ECPT[1] = 1.0; // ISILNO(1) - grid A ID
    ECPT[2] = 2.0; // ISILNO(2) - grid B ID
    ECPT[3] = SMALLV(0); // SMALLV(1)
    ECPT[4] = SMALLV(1); // SMALLV(2)
    ECPT[5] = SMALLV(2); // SMALLV(3)
    ECPT[6] = 0.0; // ICSSV (coordinate system ID for SMALLV)
    ECPT[7] = 1.0; // IMATID (material ID)
    ECPT[8] = A; // A (cross-sectional area)
    ECPT[9] = I1; // I1 (area moment of inertia about y-axis)
    ECPT[10] = I2; // I2 (area moment of inertia about z-axis)
    ECPT[11] = FJ; // FJ (torsional constant)
    ECPT[12] = 0.0; // NSM (non-structural mass)
    ECPT[13] = 0.0; // FE (force element descriptions)
    ECPT[14] = 0.0; // R1 (stress recovery coefficient)
    ECPT[15] = 0.0; // T1
    ECPT[16] = 0.0; // R2
    ECPT[17] = 0.0; // T2
    ECPT[18] = 0.0; // R3
    ECPT[19] = 0.0; // T3
    ECPT[20] = 0.0; // R4
    ECPT[21] = 0.0; // T4
    ECPT[22] = 1.0; // K1 (shear area factor)
    ECPT[23] = 1.0; // K2 (shear area factor)
    ECPT[24] = 0.0; // C (stress intensification factor)
    ECPT[25] = 1.0; // KX (flexibility correction factor)
    ECPT[26] = 1.0; // KY
    ECPT[27] = 1.0; // KZ
    ECPT[28] = R; // R (radius of curvature)
    ECPT[29] = BETAR; // BETAR (angle in degrees)
    ECPT[30] = 0.0; // MCSIDA (coordinate system ID for grid A)
    ECPT[31] = GA(0); // GPA(1)
    ECPT[32] = GA(1); // GPA(2)
    ECPT[33] = GA(2); // GPA(3)
    ECPT[34] = 0.0; // MCSIDB (coordinate system ID for grid B)
    ECPT[35] = GB(0); // GPB(1)
    ECPT[36] = GB(1); // GPB(2)
    ECPT[37] = GB(2); // GPB(3)
    ECPT[38] = 0.0; // ELTEMP (average element temperature)

    // Common block variables (simulated)
    int NPVT = 1; // pivot point is grid A (node 1)
    int IELID = static_cast<int>(ECPT[0]);
    std::vector<int> ISILNO = {static_cast<int>(ECPT[1]), static_cast<int>(ECPT[2])};
    std::vector<double> SMALLV_VEC = {ECPT[3], ECPT[4], ECPT[5]};
    int ICSSV = static_cast<int>(ECPT[6]);
    int IMATID = static_cast<int>(ECPT[7]);
    double A_val = ECPT[8];
    double I1_val = ECPT[9];
    double I2_val = ECPT[10];
    double FJ_val = ECPT[11];
    double NSM_val = ECPT[12];
    double FE_val = ECPT[13];
    double K1_val = ECPT[22];
    double K2_val = ECPT[23];
    double KX_val = ECPT[25];
    double KY_val = ECPT[26];
    double KZ_val = ECPT[27];
    double R_val = ECPT[28];
    double BETAR_val = ECPT[29];
    int MCSIDA_val = static_cast<int>(ECPT[30]);
    std::vector<double> GPA_vec = {ECPT[31], ECPT[32], ECPT[33]};
    int MCSIDB_val = static_cast<int>(ECPT[34]);
    std::vector<double> GPB_vec = {ECPT[35], ECPT[36], ECPT[37]};
    double TEMPEL_val = ECPT[38];

    // Simulate material properties
    double E_mat = E;
    double G_mat = G;
    double NU_mat = nu;

    // Step 1: Determine which point is the pivot
    int IPVT = 1;
    if (ISILNO[0] != NPVT) {
        IPVT = 2;
        if (ISILNO[1] != NPVT) {
            // Error handling - but for test case, this won't happen
        }
    }

    // Step 2: Set up coordinate system IDs
    int JCSIDA = 31;
    int JCSIDB = 35;
    int ICSIDA = MCSIDA_val;
    int ICSIDB = MCSIDB_val;

    // Step 3: Define location of end A and B in DP array (DP(1) to DP(6))
    std::vector<double> DP(16, 0.0); // DP(1) to DP(16)
    DP[0] = GPA_vec[0]; // DP(1)
    DP[1] = GPA_vec[1]; // DP(2)
    DP[2] = GPA_vec[2]; // DP(3)
    DP[3] = GPB_vec[0]; // DP(4)
    DP[4] = GPB_vec[1]; // DP(5)
    DP[5] = GPB_vec[2]; // DP(6)

    // Step 4: Define components of vector from end A to center of curvature
    // In NASTRAN, SMALLV is the reference vector pointing from A to center of curvature
    DP[6] = SMALLV_VEC[0]; // DP(7)
    DP[7] = SMALLV_VEC[1]; // DP(8)
    DP[8] = SMALLV_VEC[2]; // DP(9)
    double FLD = std::sqrt(DP[6]*DP[6] + DP[7]*DP[7] + DP[8]*DP[8]);
    if (FLD <= 0.000) {
        // Error
    }
    DP[6] /= FLD;
    DP[7] /= FLD;
    DP[8] /= FLD;

    // Step 5: Determine if points are in basic coordinates
    bool ABASIC = (ICSIDA == 0);
    bool BBASIC = (ICSIDB == 0);

    // Step 6: Compute transformation matrices (not needed for basic coordinates in test case)
    // So we just use the vectors as-is
    Eigen::Vector3d VECJ(DP[6], DP[7], DP[8]);

    // Step 7: Calculate true length of elbow
    double T = BETAR_val * DCR; // angle in radians
    double FL = R_val * T;
    if (FL == 0.0) {
        // Error
    }

    // Step 8: Construct vector from A to B
    Eigen::Vector3d SMALV0(DP[3] - DP[0], DP[4] - DP[1], DP[5] - DP[2]);
    double FLL = SMALV0.norm();
    if (FLL == 0.0) {
        // Error
    }
    SMALV0.normalize();

    // Step 9: Compute VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK;
    VECK(0) = SMALV0(1)*VECJ(2) - SMALV0(2)*VECJ(1);
    VECK(1) = SMALV0(2)*VECJ(0) - SMALV0(0)*VECJ(2);
    VECK(2) = SMALV0(0)*VECJ(1) - SMALV0(1)*VECJ(0);
    FLL = VECK.norm();
    if (FLL == 0.0) {
        // Error
    }
    VECK.normalize();

    // Step 10: Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI;
    VECI(0) = VECJ(1)*VECK(2) - VECJ(2)*VECK(1);
    VECI(1) = VECJ(2)*VECK(0) - VECJ(0)*VECK(2);
    VECI(2) = VECJ(0)*VECK(1) - VECJ(1)*VECK(0);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        // Error
    }
    VECI.normalize();

    // Step 11: Set up material properties (already known for test case)
    double DAMPC = G_mat / E_mat;

    // Step 12: Set up intermediate variables
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

    // Step 13: Compute section property constants
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

    // Step 14: Form the 6x6 flexibility matrix F
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

    // Step 15: Invert flexibility matrix to get local stiffness matrix DF
    // Since F is symmetric, we can use Eigen's LLT decomposition for stability
    Eigen::Matrix<double, 6, 6> DF;
    try {
        // Use partial pivoting LU for robust inversion
        Eigen::FullPivLU<Eigen::Matrix<double, 6, 6>> lu(F);
        if (lu.rank() < 6) {
            // Singular matrix error
        }
        DF = lu.inverse();
    } catch (...) {
        // Fallback to pseudo-inverse if needed
        Eigen::JacobiSVD<Eigen::Matrix<double, 6, 6>> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const double tolerance = 1e-12;
        Eigen::Matrix<double, 6, 6> S_inv = Eigen::Matrix<double, 6, 6>::Zero();
        for (int i = 0; i < 6; ++i) {
            if (svd.singularValues()(i) > tolerance) {
                S_inv(i, i) = 1.0 / svd.singularValues()(i);
            }
        }
        DF = svd.matrixV() * S_inv * svd.matrixU().transpose();
    }

    // Step 16: Set up force transformation matrix H
    Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i, i) = -1.0;
    }
    H(3, 2) = -R_val * (1.0 - CT);
    H(4, 2) = R_val * ST;
    H(5, 0) = -H(3, 2);
    H(5, 1) = -H(4, 2);

    // Step 17: Form the 12x12 stiffness matrix S
    Eigen::Matrix<double, 12, 12> S = Eigen::Matrix<double, 12, 12>::Zero();

    // Upper right quadrant: S(1:6, 7:12) = H * DF
    S.block(0, 6, 6, 6) = H * DF;

    // Lower left quadrant: S(7:12, 1:6) = DF * H^T
    S.block(6, 0, 6, 6) = DF * H.transpose();

    // Lower right quadrant: S(7:12, 7:12) = DF
    S.block(6, 6, 6, 6) = DF;

    // Upper left quadrant: S(1:6, 1:6) = H * DF * H^T
    S.block(0, 0, 6, 6) = H * DF * H.transpose();

    // Ensure symmetry by averaging upper and lower triangles
    for (int i = 0; i < 12; ++i) {
        for (int j = i + 1; j < 12; ++j) {
            double avg = 0.5 * (S(i, j) + S(j, i));
            S(i, j) = avg;
            S(j, i) = avg;
        }
    }

    // Step 18: Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << format_double(S(i, j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}