#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <json/json.h>
#include <sstream>
#include <string>

// Since we can't use external JSON library in strict C++17 without dependencies,
// we'll implement minimal JSON serialization for the required format.
// We'll use a simple string-based approach.

std::string double_to_scientific(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t pos = s.find_last_not_of('0');
    if (pos != std::string::npos && s[pos] == '.') {
        s.erase(pos);
    } else if (pos != std::string::npos) {
        s.erase(pos + 1);
    }
    // Remove '+' from exponent
    pos = s.find("e+");
    if (pos != std::string::npos) {
        s.replace(pos, 2, "e");
    }
    return s;
}

std::string matrix_to_json(const Eigen::MatrixXd& mat) {
    std::ostringstream json;
    json << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < mat.rows(); ++i) {
        if (i > 0) json << ",";
        json << "[";
        for (int j = 0; j < mat.cols(); ++j) {
            if (j > 0) json << ",";
            json << double_to_scientific(mat(i, j));
        }
        json << "]";
    }
    
    json << "]}";
    return json.str();
}

int main() {
    // Test case parameters
    const double PI = M_PI;
    const double DEG_TO_RAD = PI / 180.0;

    // Grid A: (1, 0, 0), Grid B: (0, 0, 1)
    Eigen::Vector3d GPA(1.0, 0.0, 0.0);
    Eigen::Vector3d GPB(0.0, 0.0, 1.0);
    
    // Reference vector V = (0, 1, 0)
    Eigen::Vector3d SMALLV(0.0, 1.0, 0.0);
    
    // Geometry
    double R = 1.0;                    // Radius of curvature
    double BETAR = 90.0;               // Angle in degrees
    
    // Material properties
    double E = 200e9;                  // Pa
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu)); // Shear modulus
    
    // Cross-section: circular pipe
    double ro = 0.05;                  // outer radius
    double t = 0.005;                  // wall thickness
    double ri = ro - t;                // inner radius
    
    double A = PI * (ro*ro - ri*ri);   // Area
    double I1 = PI/4.0 * (ro*ro*ro*ro - ri*ri*ri*ri); // Moment of inertia about y-axis
    double I2 = I1;                    // Moment of inertia about z-axis (symmetric)
    double J = 2.0 * I1;               // Torsional constant (polar moment)
    
    // Flexibility correction factors (KX, KY, KZ) - set to 1.0 as per typical usage
    double KX = 1.0;
    double KY = 1.0;
    double KZ = 1.0;
    
    // Shear area factors (K1, K2) - for circular hollow section, typically 0.5
    // According to Roark: "0.5 FOR THIN-WALLED HOOLOW CIRCULAR SECTION"
    double K1 = 0.5;
    double K2 = 0.5;
    
    // Convert angle to radians
    double T = BETAR * DEG_TO_RAD;
    double ST = std::sin(T);
    double CT = std::cos(T);
    double S2T = std::sin(2.0*T);
    double C2T = std::cos(2.0*T);
    
    // Intermediate variables
    double RA = R / (A * E);
    double RV1 = R / (2.0 * K1 * G * A);
    double RV2 = (K1/K2) * RV1;
    double RT = R / (G * J * 2.0);
    double RB0 = R / (E * I2 * 2.0);
    double RB1 = R / (E * I1);
    double R2 = R * R;
    
    // Compute local coordinate system vectors
    // Vector from A to center of curvature: same as reference vector SMALLV
    Eigen::Vector3d VECJ = SMALLV.normalized();
    
    // Vector from A to B
    Eigen::Vector3d SMALV0 = GPB - GPA;
    double FLL = SMALV0.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Zero length between nodes" << std::endl;
        return 1;
    }
    SMALV0.normalize();
    
    // Compute VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK = SMALV0.cross(VECJ);
    FLL = VECK.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Cross product is zero" << std::endl;
        return 1;
    }
    VECK.normalize();
    
    // Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI = VECJ.cross(VECK);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        std::cerr << "Error: Cross product is zero" << std::endl;
        return 1;
    }
    VECI.normalize();
    
    // Build 6x6 flexibility matrix F
    Eigen::MatrixXd F = Eigen::MatrixXd::Zero(6, 6);
    
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
    F(2,3) += RT * R * (ST - T*CT);
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
    F(0,5) -= RB1 * R * (ST - T*CT);
    F(5,0) = F(0,5);
    F(1,5) += RB1 * R * (T*ST + CT - 1.0);
    F(5,1) = F(1,5);
    F(2,3) += RB0 * R * (ST - T*CT);
    F(3,2) = F(2,3);
    F(2,4) -= RB0 * R * T * ST;
    F(4,2) = F(2,4);
    F(3,4) -= 0.50 * RB0 * (1.0 - C2T);
    F(4,3) = F(3,4);
    
    // Make F symmetric
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < i; ++j) {
            F(i,j) = F(j,i);
        }
    }
    
    // Invert flexibility matrix to get stiffness submatrix
    // Use Eigen's LU decomposition for inversion
    Eigen::FullPivLU<Eigen::MatrixXd> lu(F);
    if (!lu.isInvertible()) {
        std::cerr << "Error: Flexibility matrix is singular" << std::endl;
        return 1;
    }
    Eigen::MatrixXd DF = lu.inverse();
    
    // Build 12x12 stiffness matrix S
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(12, 12);
    
    // Set bottom-right 6x6 block
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i+6, j+6) = DF(i, j);
        }
    }
    
    // Build H matrix (6x6)
    Eigen::MatrixXd H = Eigen::MatrixXd::Zero(6, 6);
    for (int i = 0; i < 6; ++i) {
        H(i, i) = -1.0;
    }
    H(3, 2) = -R * (1.0 - CT);  // H(4,3) in Fortran (1-indexed)
    H(4, 2) = R * ST;          // H(5,3) in Fortran
    H(5, 0) = -H(3, 2);        // H(6,1) in Fortran
    H(5, 1) = -H(4, 2);        // H(6,2) in Fortran
    
    // Compute S(1:6,7:12) = H * DF
    Eigen::MatrixXd S_top_right = H * DF;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j+6) = S_top_right(i, j);
        }
    }
    
    // Compute S(1:6,1:6) = S(1:6,7:12) * H^T = (H * DF) * H^T
    Eigen::MatrixXd S_top_left = S_top_right * H.transpose();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j) = S_top_left(i, j);
        }
    }
    
    // Make symmetric: S(7:12,1:6) = S(1:6,7:12)^T
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i+6, j) = S(j, i+6);
        }
    }
    
    // Output as JSON
    std::cout << matrix_to_json(S) << std::endl;
    
    return 0;
}