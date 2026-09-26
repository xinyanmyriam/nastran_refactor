#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Helper function to convert double to scientific notation string with fixed precision
std::string to_scientific(double x, int precision = 6) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Replace +e+ with e
    size_t epos = s.find("e+");
    if (epos != std::string::npos) {
        s.replace(epos, 2, "e");
    }
    return s;
}

int main() {
    // Test case parameters
    const double PI = 3.14159265358979323846;
    const double DCR = 0.017453292; // degrees to radians conversion factor

    // Grid points
    Eigen::Vector3d GA(1.0, 0.0, 0.0); // Grid A
    Eigen::Vector3d GB(0.0, 0.0, 1.0); // Grid B
    Eigen::Vector3d V(0.0, 1.0, 0.0);  // Reference vector

    // Geometry
    double R = 1.0;           // Radius of curvature
    double BETAR = 90.0;      // Angle in degrees

    // Material properties
    double E = 200e9;         // Pa
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));

    // Cross-section (circular pipe)
    double ro = 0.05;         // outer radius
    double t = 0.005;         // wall thickness
    double ri = ro - t;       // inner radius
    double A = PI * (ro*ro - ri*ri);
    double I1 = PI/4.0 * (ro*ro*ro*ro - ri*ri*ri*ri); // Iy = Iz for circular
    double I2 = I1;           // I1 and I2 are symmetric for circular section
    double J = 2.0 * I1;      // torsional constant for circular section

    // Initialize ECPT array (100 elements, double precision)
    std::vector<double> ECPT(100, 0.0);

    // Fill ECPT with test case values
    ECPT[0] = 1.0; // IELID (element ID)
    
    // ISILNO(2): grid point IDs (we'll use 1 and 2 for simplicity)
    ECPT[1] = 1.0; // ISILNO(1)
    ECPT[2] = 2.0; // ISILNO(2)
    
    // SMALLV(3): reference vector components
    ECPT[3] = V(0);
    ECPT[4] = V(1);
    ECPT[5] = V(2);
    
    // ICSSV: coordinate system ID for SMALLV (0 for basic)
    ECPT[6] = 0.0;
    
    // IMATID: material ID
    ECPT[7] = 1.0;
    
    // A, I1, I2, FJ (torsional constant), NSM
    ECPT[8] = A;
    ECPT[9] = I1;
    ECPT[10] = I2;
    ECPT[11] = J;
    ECPT[12] = 0.0; // NSM
    
    // FE, C1, C2, D1, D2, F1, F2, G1, G2, K1, K2, C
    // Set K1=K2=1.0 for circular section (shear correction factors)
    ECPT[13] = 0.0; // FE
    ECPT[14] = 0.0; // C1
    ECPT[15] = 0.0; // C2
    ECPT[16] = 0.0; // D1
    ECPT[17] = 0.0; // D2
    ECPT[18] = 0.0; // F1
    ECPT[19] = 0.0; // F2
    ECPT[20] = 0.0; // G1
    ECPT[21] = 0.0; // G2
    ECPT[22] = 1.0; // K1 (shear correction)
    ECPT[23] = 1.0; // K2 (shear correction)
    ECPT[24] = 0.0; // C (stress intensification)
    
    // KX, KY, KZ: flexibility correction factors
    ECPT[25] = 1.0; // KX
    ECPT[26] = 1.0; // KY
    ECPT[27] = 1.0; // KZ
    
    // R, BETAR, MCSIDA, GPA(3), MCSIDB, GPB(3), ELTEMP
    ECPT[28] = R;
    ECPT[29] = BETAR;
    ECPT[30] = 0.0; // MCSIDA (basic coordinate system)
    ECPT[31] = GA(0);
    ECPT[32] = GA(1);
    ECPT[33] = GA(2);
    ECPT[34] = 0.0; // MCSIDB (basic coordinate system)
    ECPT[35] = GB(0);
    ECPT[36] = GB(1);
    ECPT[37] = GB(2);
    ECPT[38] = 0.0; // ELTEMP

    // NPVT: pivot node (we'll use node 1 as pivot)
    int NPVT = 1;

    // ISILNO(2) from ECPT
    int ISILNO[2] = {static_cast<int>(ECPT[1]), static_cast<int>(ECPT[2])};

    // Determine which point is the pivot
    int IPVT = 1;
    if (ISILNO[0] != NPVT) {
        IPVT = 2;
        if (ISILNO[1] != NPVT) {
            std::cerr << "Error: NPVT not found in ISILNO" << std::endl;
            return 1;
        }
    }

    // JCSIDA = 31, JCSIDB = 35 (Fortran indexing, but we use 0-based)
    int JCSIDA = 30; // 31-1
    int JCSIDB = 34; // 35-1

    // ICSIDA, ICSIDB
    int ICSIDA = static_cast<int>(ECPT[JCSIDA]);
    int ICSIDB = static_cast<int>(ECPT[JCSIDB]);

    // DP(1) through DP(6): coordinates of A and B
    std::vector<double> DP(16, 0.0);
    DP[0] = ECPT[JCSIDA + 1]; // GA_x
    DP[1] = ECPT[JCSIDA + 2]; // GA_y
    DP[2] = ECPT[JCSIDA + 3]; // GA_z
    DP[3] = ECPT[JCSIDB + 1]; // GB_x
    DP[4] = ECPT[JCSIDB + 2]; // GB_y
    DP[5] = ECPT[JCSIDB + 3]; // GB_z

    // DP(7) through DP(9): components of vector from A to center of curvature
    // In this test case, center of curvature is at origin (0,0,0) since GA=(1,0,0), GB=(0,0,1), R=1
    // So vector from A to center is (-1,0,0)
    DP[6] = -1.0;
    DP[7] = 0.0;
    DP[8] = 0.0;

    // Normalize DP(7-9)
    double FLD = std::sqrt(DP[6]*DP[6] + DP[7]*DP[7] + DP[8]*DP[8]);
    if (FLD <= 0.0) {
        std::cerr << "Error: FLD <= 0.0" << std::endl;
        return 1;
    }
    DP[6] /= FLD;
    DP[7] /= FLD;
    DP[8] /= FLD;

    // ABASIC, BBASIC
    bool ABASIC = (ICSIDA == 0);
    bool BBASIC = (ICSIDB == 0);

    // VECJ: normalized vector from A to center
    Eigen::Vector3d VECJ(DP[6], DP[7], DP[8]);

    // SMALV0: vector from A to B
    Eigen::Vector3d SMALV0(DP[3]-DP[0], DP[4]-DP[1], DP[5]-DP[2]);
    double FLL = SMALV0.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: FLL == 0.0" << std::endl;
        return 1;
    }
    SMALV0.normalize();

    // VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK = SMALV0.cross(VECJ);
    FLL = VECK.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: VECK norm == 0.0" << std::endl;
        return 1;
    }
    VECK.normalize();

    // VECI = VECJ × VECK
    Eigen::Vector3d VECI = VECJ.cross(VECK);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: VECI norm == 0.0" << std::endl;
        return 1;
    }
    VECI.normalize();

    // Material properties
    double MATIDC = ECPT[7];
    double ELTEMP = ECPT[38];
    double E_val = E;
    double G_val = G;
    double NU_val = nu;

    // Flexibility correction factors
    double KX = ECPT[25];
    double KY = ECPT[26];
    double KZ = ECPT[27];
    if (KX < 1.0e-8) KX = 1.0;
    if (KY < 1.0e-8) KY = 1.0;
    if (KZ < 1.0e-8) KZ = 1.0;

    double FI1 = I1 / KZ;
    double FI2 = I2 / KY;
    double FJK = J / KX;

    // Shear correction factors K1, K2
    double K1 = ECPT[22];
    double K2 = ECPT[23];
    if (K1 < 1.0e-8) K1 = 1.0;
    if (K2 < 1.0e-8) K2 = 1.0;
    if (K1 > 1.0) K1 = 1.0 / K1;
    if (K2 > 1.0) K2 = 1.0 / K2;

    // Section property constants
    double T = BETAR * DCR; // angle in radians
    double RA = R / (A * E_val);
    double RV1 = R / (2.0 * K1 * G_val * A);
    double RV2 = (K1 / K2) * RV1;
    double RT = R / (G_val * FJK * 2.0);
    double RB0 = R / (E_val * FI2 * 2.0);
    double RB1 = R / (E_val * FI1);
    double R2 = R * R;

    // Trigonometric constants
    double ST = std::sin(T);
    double CT = std::cos(T);
    double S2T = std::sin(2.0 * T);
    double C2T = std::cos(2.0 * T);

    // Initialize 6x6 flexibility matrix F
    Eigen::Matrix<double, 6, 6> F = Eigen::Matrix<double, 6, 6>::Zero();

    // Axial terms
    F(0,0) += 0.25 * RA * (2.0*T + S2T);
    F(1,1) += 0.25 * RA * (2.0*T - S2T);
    F(0,1) += 0.50 * RA * ST * ST;
    F(1,0) = F(0,1);

    // Shear terms
    F(0,0) += 0.5 * RV1 * (2.0*T - S2T);
    F(1,1) += 0.5 * RV1 * (2.0*T + S2T);
    F(2,2) += 2.0 * RV2 * T;
    F(0,1) -= RV1 * ST * ST;
    F(1,0) = F(0,1);

    // Torsion terms
    F(2,2) += 0.5 * RT * R2 * (6.0*T + S2T - 8.0*ST);
    F(3,3) += 0.5 * RT * (2.0*T + S2T);
    F(4,4) += 0.5 * RT * (2.0*T - S2T);
    F(2,3) += RT * R * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) += RT * R * (2.0 - 2.0*CT - T*ST);
    F(4,2) = F(2,4);
    F(3,4) += 0.5 * RT * (1.0 - C2T);
    F(4,3) = F(3,4);

    // Bending terms
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
            if (i != j) {
                F(j,i) = F(i,j);
            }
        }
    }

    // Invert flexibility matrix to get local stiffness matrix DF (6x6)
    Eigen::Matrix<double, 6, 6> DF;
    try {
        DF = F.inverse();
    } catch (...) {
        std::cerr << "Error: Flexibility matrix is singular" << std::endl;
        return 1;
    }

    // H matrix (6x6)
    Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i,i) = -1.0;
    }
    H(3,2) = -R * (1.0 - CT);
    H(4,2) = R * ST;
    H(5,0) = R * (1.0 - CT);
    H(5,1) = -R * ST;

    // Build 12x12 stiffness matrix S
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

    // Make symmetric
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
            std::cout << to_scientific(S(i,j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}