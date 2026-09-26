#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Use double precision throughout
using Real = double;

// Helper function to format double in scientific notation
std::string format_scientific(Real x) {
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
    // Grid A: (1, 0, 0), Grid B: (0, 0, 1)
    // Reference vector V = (0, 1, 0)
    // Radius of curvature R = 1.0, Angle = 90 degrees
    // E = 200e9, nu = 0.3, G = E/(2*(1+nu))
    // Circular pipe: outer radius = 0.05, wall thickness = 0.005

    const Real R = 1.0;
    const Real BETAR = 90.0; // degrees
    const Real E = 200e9;
    const Real nu = 0.3;
    const Real G = E / (2.0 * (1.0 + nu));
    
    const Real ro = 0.05;
    const Real t = 0.005;
    const Real ri = ro - t;
    
    // Ensure M_PI is available — _USE_MATH_DEFINES above guarantees it for <cmath>
    const Real A = M_PI * (ro * ro - ri * ri);
    const Real I1 = M_PI / 4.0 * (ro * ro * ro * ro - ri * ri * ri * ri);
    const Real I2 = I1; // circular symmetry
    const Real J = 2.0 * I1;

    // Convert angle to radians
    const Real DCR = 0.017453292; // degrees to radians
    const Real T = BETAR * DCR;

    // Define geometry
    Eigen::Vector3d GPA = Eigen::Vector3d(1.0, 0.0, 0.0); // Grid A
    Eigen::Vector3d GPB = Eigen::Vector3d(0.0, 0.0, 1.0); // Grid B
    Eigen::Vector3d SMALLV = Eigen::Vector3d(0.0, 1.0, 0.0); // reference vector

    // Set up ECPT array (Fortran common block equivalent)
    std::vector<Real> ECPT(100, 0.0);

    // Fill ECPT with test case values
    // ECPT(1)  - IELID         ELEMENT ID
    ECPT[0] = 1.0;
    // ECPT(2-3) - ISILNO(2)     grid point IDs (we'll use 1 and 2)
    ECPT[1] = 1.0;
    ECPT[2] = 2.0;
    // ECPT(4-6) - SMALLV(3)     reference vector
    ECPT[3] = SMALLV(0);
    ECPT[4] = SMALLV(1);
    ECPT[5] = SMALLV(2);
    // ECPT(7)  - ICSSV         coordinate system ID for SMALLV (0 for basic)
    ECPT[6] = 0.0;
    // ECPT(8)  - IMATID        material ID
    ECPT[7] = 1.0;
    // ECPT(9)  - A             cross-sectional area
    ECPT[8] = A;
    // ECPT(10-11) - I1, I2     area moments of inertia
    ECPT[9] = I1;
    ECPT[10] = I2;
    // ECPT(12) - FJ            torsional constant (J)
    ECPT[11] = J;
    // ECPT(13) - NSM           non-structural mass (0)
    ECPT[12] = 0.0;
    // ECPT(14) - FE            force element descriptions (0)
    ECPT[13] = 0.0;
    // ECPT(15-22) - stress recovery coefficients (not used)
    // ECPT(23-24) - K1, K2    area factors for shear (use 1.0)
    ECPT[22] = 1.0;
    ECPT[23] = 1.0;
    // ECPT(25) - C             stress intensification factor (0)
    ECPT[24] = 0.0;
    // ECPT(26-28) - KX, KY, KZ flexibility correction factors (1.0)
    ECPT[25] = 1.0;
    ECPT[26] = 1.0;
    ECPT[27] = 1.0;
    // ECPT(29) - R             radius of curvature
    ECPT[28] = R;
    // ECPT(30) - BETAR         angle from GA to GB (in degrees)
    ECPT[29] = BETAR;
    // ECPT(31) - MCSIDA        coord sys ID for grid point A (0 for basic)
    ECPT[30] = 0.0;
    // ECPT(32-34) - GPA(3)     basic coordinates for grid point A
    ECPT[31] = GPA(0);
    ECPT[32] = GPA(1);
    ECPT[33] = GPA(2);
    // ECPT(35) - MCSIDB        coord sys ID for grid point B (0 for basic)
    ECPT[34] = 0.0;
    // ECPT(36-38) - GPB(3)     basic coordinates for grid point B
    ECPT[35] = GPB(0);
    ECPT[36] = GPB(1);
    ECPT[37] = GPB(2);
    // ECPT(39) - ELTEMP        avg element temperature (0)
    ECPT[38] = 0.0;

    // Set up common block variables (simplified for this test)
    int NPVT = 1; // pivot node is node 1 (grid point 1)
    int ISILNO[2] = {1, 2};
    Real SMALLV_arr[3] = {SMALLV(0), SMALLV(1), SMALLV(2)};
    int ICSSV = 0;
    int IMATID = 1;
    Real A_val = A;
    Real I1_val = I1;
    Real I2_val = I2;
    Real FJ_val = J;
    Real NSM_val = 0.0;
    Real FE_val = 0.0;
    Real K1_val = 1.0;
    Real K2_val = 1.0;
    Real C_val = 0.0;
    Real KX_val = 1.0;
    Real KY_val = 1.0;
    Real KZ_val = 1.0;
    Real R_val = R;
    Real BETAR_val = BETAR;
    int MCSIDA_val = 0;
    Real GPA_arr[3] = {GPA(0), GPA(1), GPA(2)};
    int MCSIDB_val = 0;
    Real GPB_arr[3] = {GPB(0), GPB(1), GPB(2)};
    Real TEMPEL_val = 0.0;

    // Material properties
    Real E_mat = E;
    Real G_mat = G;
    Real NU_mat = nu;

    // Compute trigonometric values
    auto SID = [&](Real x) -> Real { return std::sin(x * DCR); };
    auto COD = [&](Real x) -> Real { return std::cos(x * DCR); };
    auto DTR = [&](Real x) -> Real { return x * DCR; };

    Real ST = SID(BETAR_val);
    Real CT = COD(BETAR_val);
    Real S2T = SID(2.0 * BETAR_val);
    Real C2T = COD(2.0 * BETAR_val);

    // Determine which point is the pivot
    int IPVT = 1;
    if (ISILNO[0] != NPVT) {
        IPVT = 2;
        if (ISILNO[1] != NPVT) {
            // Error case - not handled in test
        }
    }

    // Set up pointers to coordinate system IDs
    int JCSIDA = 31; // Fortran indexing, so 31 means index 30 in 0-based
    int JCSIDB = 35; // Fortran indexing, so 35 means index 34 in 0-based
    int ICSIDA = static_cast<int>(ECPT[30]); // MCSIDA at index 30
    int ICSIDB = static_cast<int>(ECPT[34]); // MCSIDB at index 34

    // Define location of end A, end B in terms of DP(1) through DP(6)
    std::vector<Real> DP(16, 0.0);
    DP[0] = ECPT[JCSIDA]; // GPA(1) - but in our case, ECPT[31] is GPA(1)
    DP[1] = ECPT[JCSIDA + 1]; // GPA(2) - ECPT[32]
    DP[2] = ECPT[JCSIDA + 2]; // GPA(3) - ECPT[33]
    DP[3] = ECPT[JCSIDB]; // GPB(1) - ECPT[35]
    DP[4] = ECPT[JCSIDB + 1]; // GPB(2) - ECPT[36]
    DP[5] = ECPT[JCSIDB + 2]; // GPB(3) - ECPT[37]

    // In our test case, we have basic coordinates, so:
    DP[0] = GPA(0);
    DP[1] = GPA(1);
    DP[2] = GPA(2);
    DP[3] = GPB(0);
    DP[4] = GPB(1);
    DP[5] = GPB(2);

    // Define components of vector from end A to center of curvature, C
    // In the Fortran code, this is SMALLV(3) at ECPT(4-6), which is our reference vector
    DP[6] = ECPT[3]; // SMALLV(1)
    DP[7] = ECPT[4]; // SMALLV(2)
    DP[8] = ECPT[5]; // SMALLV(3)

    Real FLD = std::sqrt(DP[6] * DP[6] + DP[7] * DP[7] + DP[8] * DP[8]);
    if (FLD <= 0.0) {
        // Error case
    }
    DP[6] /= FLD;
    DP[7] /= FLD;
    DP[8] /= FLD;

    // Determine if points are in basic coordinates
    bool ABASIC = (ICSIDA == 0);
    bool BBASIC = (ICSIDB == 0);

    // Since both are basic in our test, we skip transformation
    Eigen::Vector3d VECJ(DP[6], DP[7], DP[8]);

    // Calculate true length of elbow
    Real FL = R_val * DTR(BETAR_val);
    if (FL == 0.0) {
        // Error case
    }

    // Construct vector from A to B
    Eigen::Vector3d SMALV0(DP[3] - DP[0], DP[4] - DP[1], DP[5] - DP[2]);
    Real FLL = SMALV0.norm();
    if (FLL == 0.0) {
        // Error case
    }
    SMALV0.normalize();

    // Compute VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK;
    VECK(0) = SMALV0(1) * VECJ(2) - SMALV0(2) * VECJ(1);
    VECK(1) = SMALV0(2) * VECJ(0) - SMALV0(0) * VECJ(2);
    VECK(2) = SMALV0(0) * VECJ(1) - SMALV0(1) * VECJ(0);
    FLL = VECK.norm();
    if (FLL == 0.0) {
        // Error case
    }
    VECK.normalize();

    // Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI;
    VECI(0) = VECJ(1) * VECK(2) - VECJ(2) * VECK(1);
    VECI(1) = VECJ(2) * VECK(0) - VECJ(0) * VECK(2);
    VECI(2) = VECJ(0) * VECK(1) - VECJ(1) * VECK(0);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        // Error case
    }
    VECI.normalize();

    // Set up intermediate variables for stiffness matrix calculation
    Real KX = KX_val;
    Real KY = KY_val;
    Real KZ = KZ_val;
    if (KX < 1.0e-8) KX = 1.0;
    if (KY < 1.0e-8) KY = 1.0;
    if (KZ < 1.0e-8) KZ = 1.0;
    Real FI1 = I1_val / KZ;
    Real FI2 = I2_val / KY;
    Real FJK = FJ_val / KX;

    // Area factors for shear
    Real K1 = K1_val;
    Real K2 = K2_val;
    if (K1 < 1.0e-8) K1 = 1.0;
    if (K2 < 1.0e-8) K2 = 1.0;
    if (K1 > 1.0) K1 = 1.0 / K1;
    if (K2 > 1.0) K2 = 1.0 / K2;

    // Compute section property constants
    Real RA = R_val / (A_val * E_mat);
    Real RV1 = R_val / (2.0 * K1 * G_mat * A_val);
    Real RV2 = (K1 / K2) * RV1;
    Real RT = R_val / (G_mat * FJK * 2.0);
    Real RB0 = R_val / (E_mat * FI2 * 2.0);
    Real RB1 = R_val / (E_mat * FI1);
    Real R2 = R_val * R_val;

    // Form the node flexibility matrix F(6,6)
    Eigen::Matrix<Real, 6, 6> F = Eigen::Matrix<Real, 6, 6>::Zero();

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

    // Make F symmetric
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < i; ++j) {
            F(j,i) = F(i,j);
        }
    }

    // Invert flexibility matrix to get local stiffness matrix DF
    Eigen::Matrix<Real, 6, 6> DF;
    try {
        DF = F.inverse();
    } catch (...) {
        // In practice, we'd handle singularity, but for this test case it should be fine
        DF = F.completeOrthogonalDecomposition().pseudoInverse();
    }

    // Set up H matrix (force transformation)
    Eigen::Matrix<Real, 6, 6> H = Eigen::Matrix<Real, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i,i) = -1.0;
    }
    H(3,2) = -R_val * (1.0 - CT);
    H(4,2) = R_val * ST;
    H(5,0) = -H(3,2);
    H(5,1) = -H(4,2);

    // Form the 12x12 stiffness matrix S
    Eigen::Matrix<Real, 12, 12> S = Eigen::Matrix<Real, 12, 12>::Zero();

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

    // S(1:6,1:6) = (H * DF) * H^T
    Eigen::Matrix<Real, 6, 6> temp2 = temp1 * H.transpose();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j) = temp2(i,j);
        }
    }

    // Reflect for symmetry (S is already symmetric due to construction, but ensure)
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
            std::cout << format_scientific(S(i,j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}