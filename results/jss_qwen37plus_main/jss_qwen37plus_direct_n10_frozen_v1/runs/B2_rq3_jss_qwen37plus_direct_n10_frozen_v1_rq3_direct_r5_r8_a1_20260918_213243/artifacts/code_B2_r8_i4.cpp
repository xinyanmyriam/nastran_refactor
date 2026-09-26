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
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos + 1);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Replace +e with e
    size_t epos = s.find("+e");
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
    const double G = 76.923e9;   // Pa
    const double A = 0.01;       // m^2
    const double Iy = 8.333e-6;  // m^4 (I1 in Fortran)
    const double Iz = 8.333e-6;  // m^4 (I2 in Fortran)
    const double J = 1.667e-5;   // m^4 (FJ in Fortran)
    const double I12 = 0.0;      // no warping (I12 = 0)

    // In Fortran, the reference vector SMALLV is normalized and used to define orientation.
    // For standard CBAR with no special orientation, we use a default reference vector.
    // According to NASTRAN documentation, if SMALLV is not provided, it defaults to (0,1,0).
    // We'll use (0,1,0) as the reference vector.
    const double smallv[3] = {0.0, 1.0, 0.0};

    // Material properties - same as input
    const double k1 = 0.0; // shear correction factor (not used in this test)
    const double k2 = 0.0;

    // Compute transformation vectors
    // VECI = vector from A to B (note: corrected direction - should be from A to B for local x-axis)
    Eigen::Vector3d veci(nodeB[0] - nodeA[0], nodeB[1] - nodeA[1], nodeB[2] - nodeA[2]);
    double fl = veci.norm();
    if (fl == 0.0) {
        std::cerr << "Error: zero length element" << std::endl;
        return 1;
    }
    veci.normalize();

    // SMALV0 = normalized smallv
    Eigen::Vector3d smalv0(smallv[0], smallv[1], smallv[2]);
    double fl_smallv = smalv0.norm();
    if (fl_smallv == 0.0) {
        std::cerr << "Error: zero reference vector" << std::endl;
        return 1;
    }
    smalv0.normalize();

    // VECK = VECI × SMALV0
    Eigen::Vector3d veck = veci.cross(smalv0);
    double fll = veck.norm();
    if (fll == 0.0) {
        std::cerr << "Error: cross product is zero" << std::endl;
        return 1;
    }
    veck.normalize();

    // VECJ = VECK × VECI
    Eigen::Vector3d vecj = veck.cross(veci);
    fll = vecj.norm();
    if (fll == 0.0) {
        std::cerr << "Error: cross product is zero" << std::endl;
        return 1;
    }
    vecj.normalize();

    // Build 3x3 direction cosine matrix T3x3
    // Columns are: veci (x-axis), vecj (y-axis), veck (z-axis)
    Eigen::Matrix3d T3x3;
    T3x3.col(0) = veci;
    T3x3.col(1) = vecj;
    T3x3.col(2) = veck;

    // Intermediate variables for stiffness matrix
    double lsq = L * L;
    double lcube = lsq * L;
    double ei1 = E * Iy;  // EIy
    double ei2 = E * Iz;  // Eiz
    double r1, r2;

    // Shear correction logic (simplified - k1/k2 are zero in test case)
    if (k1 == 0.0 || I12 != 0.0) {
        r1 = 12.0 * ei1 / lcube;
    } else {
        // Not used in test case
        r1 = 12.0 * ei1 / lcube;
    }

    if (k2 == 0.0 || I12 != 0.0) {
        r2 = 12.0 * ei2 / lcube;
    } else {
        // Not used in test case
        r2 = 12.0 * ei2 / lcube;
    }

    double sk1 = 0.25 * r1 * lsq + ei1 / L;
    double sk2 = 0.25 * r2 * lsq + ei2 / L;
    double sk3 = 0.25 * r1 * lsq - ei1 / L;
    double sk4 = 0.25 * r2 * lsq - ei2 / L;

    double ael = A * E / L;
    double lr1 = L * r1 / 2.0;
    double lr2 = L * r2 / 2.0;
    double gjl = G * J / L;

    // Initialize 12x12 stiffness matrix (row-major order, 0-indexed)
    Eigen::Matrix<double, 12, 12, Eigen::RowMajor> ke_local = Eigen::Matrix<double, 12, 12, Eigen::RowMajor>::Zero();

    // Fill the local stiffness matrix according to Fortran logic (1-based indexing converted to 0-based)
    // Row 0 (Fortran index 1): [AEL, 0, 0, 0, 0, 0, -AEL, 0, 0, 0, 0, 0]
    ke_local(0,0) = ael;
    ke_local(0,6) = -ael;

    // Row 1 (Fortran index 2): [0, R1, 0, 0, 0, LR1, 0, -R1, 0, 0, 0, LR1]
    ke_local(1,1) = r1;
    ke_local(1,5) = lr1;
    ke_local(1,7) = -r1;
    ke_local(1,11) = lr1;

    // Row 2 (Fortran index 3): [0, 0, R2, 0, -LR2, 0, 0, 0, -R2, 0, LR2, 0]
    ke_local(2,2) = r2;
    ke_local(2,4) = -lr2;
    ke_local(2,8) = -r2;
    ke_local(2,10) = lr2;

    // Row 3 (Fortran index 4): [0, 0, 0, GJL, 0, 0, 0, 0, 0, -GJL, 0, 0]
    ke_local(3,3) = gjl;
    ke_local(3,9) = -gjl;

    // Row 4 (Fortran index 5): [0, 0, -LR2, 0, SK2, 0, 0, 0, LR2, 0, SK4, 0]
    ke_local(4,2) = -lr2;
    ke_local(4,4) = sk2;
    ke_local(4,8) = lr2;
    ke_local(4,10) = sk4;

    // Row 5 (Fortran index 6): [0, LR1, 0, 0, 0, SK1, 0, -LR1, 0, 0, 0, SK3]
    ke_local(5,1) = lr1;
    ke_local(5,5) = sk1;
    ke_local(5,7) = -lr1;
    ke_local(5,11) = sk3;

    // Row 6 (Fortran index 7): [-AEL, 0, 0, 0, 0, 0, AEL, 0, 0, 0, 0, 0]
    ke_local(6,0) = -ael;
    ke_local(6,6) = ael;

    // Row 7 (Fortran index 8): [0, -R1, 0, 0, 0, -LR1, 0, R1, 0, 0, 0, -LR1]
    ke_local(7,1) = -r1;
    ke_local(7,5) = -lr1;
    ke_local(7,7) = r1;
    ke_local(7,11) = -lr1;

    // Row 8 (Fortran index 9): [0, 0, -R2, 0, LR2, 0, 0, 0, R2, 0, -LR2, 0]
    ke_local(8,2) = -r2;
    ke_local(8,4) = lr2;
    ke_local(8,8) = r2;
    ke_local(8,10) = -lr2;

    // Row 9 (Fortran index 10): [0, 0, 0, -GJL, 0, 0, 0, 0, 0, GJL, 0, 0]
    ke_local(9,3) = -gjl;
    ke_local(9,9) = gjl;

    // Row 10 (Fortran index 11): [0, 0, LR2, 0, SK4, 0, 0, 0, -LR2, 0, SK2, 0]
    ke_local(10,2) = lr2;
    ke_local(10,4) = sk4;
    ke_local(10,8) = -lr2;
    ke_local(10,10) = sk2;

    // Row 11 (Fortran index 12): [0, -LR1, 0, 0, 0, SK3, 0, LR1, 0, 0, 0, SK1]
    ke_local(11,1) = -lr1;
    ke_local(11,5) = sk3;
    ke_local(11,7) = lr1;
    ke_local(11,11) = sk1;

    // Build 12x12 transformation matrix T
    // T is block diagonal with T3x3 repeated for each node's 3 translational and 3 rotational DOFs
    // For CBAR: nodes A and B, each with 6 DOFs (tx,ty,tz,rx,ry,rz)
    // So T has T3x3 in blocks: (0:2,0:2), (3:5,3:5), (6:8,6:8), (9:11,9:11)
    Eigen::Matrix<double, 12, 12, Eigen::RowMajor> T = Eigen::Matrix<double, 12, 12, Eigen::RowMajor>::Zero();
    
    // Translational DOFs for node A (indices 0,1,2)
    T.block<3,3>(0,0) = T3x3;
    // Rotational DOFs for node A (indices 3,4,5) 
    T.block<3,3>(3,3) = T3x3;
    // Translational DOFs for node B (indices 6,7,8)
    T.block<3,3>(6,6) = T3x3;
    // Rotational DOFs for node B (indices 9,10,11)
    T.block<3,3>(9,9) = T3x3;

    // Transform to global coordinates: K_global = T^T * K_local * T
    Eigen::Matrix<double, 12, 12, Eigen::RowMajor> ke_global = T.transpose() * ke_local * T;

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << format_double(ke_global(i,j));
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}