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

// Helper function to format double in scientific notation with fixed precision
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

// Invert a 6x6 matrix using Eigen's LU decomposition
bool invert_6x6(const Eigen::Matrix<Real, 6, 6>& A, Eigen::Matrix<Real, 6, 6>& invA) {
    Eigen::FullPivLU<Eigen::Matrix<Real, 6, 6>> lu(A);
    if (lu.rank() < 6) {
        return false;
    }
    invA = lu.inverse();
    return true;
}

// Compute the 12x12 stiffness matrix for the curved pipe elbow element
Eigen::Matrix<Real, 12, 12> compute_elbow_stiffness() {
    // Test case parameters
    // Grid A: (1, 0, 0), Grid B: (0, 0, 1)
    Eigen::Vector3d GPA = {1.0, 0.0, 0.0};
    Eigen::Vector3d GPB = {0.0, 0.0, 1.0};
    
    // Reference vector V = (0, 1, 0)
    Eigen::Vector3d SMALLV = {0.0, 1.0, 0.0};
    
    // Radius of curvature R = 1.0, Angle = 90 degrees = π/2 radians
    Real R = 1.0;
    Real BETAR = 90.0; // degrees
    
    // Material properties
    Real E = 200e9;      // Pa
    Real nu = 0.3;
    Real G = E / (2.0 * (1.0 + nu));
    
    // Circular pipe geometry
    Real ro = 0.05;      // outer radius
    Real t = 0.005;      // wall thickness
    Real ri = ro - t;    // inner radius
    
    // Cross-section properties
    Real A = M_PI * (ro * ro - ri * ri);
    Real I1 = M_PI / 4.0 * (ro * ro * ro * ro - ri * ri * ri * ri); // Iy
    Real I2 = I1; // Iz (symmetric circular section)
    Real FJ = 2.0 * I1; // J = polar moment = Iy + Iz = 2*I1
    
    // Flexibility correction factors (KX, KY, KZ) - use default 1.0
    Real KX = 1.0;
    Real KY = 1.0;
    Real KZ = 1.0;
    
    // Shear area factors (K1, K2) - for thin-walled circular: K1=K2=0.5
    // But per NASTRAN doc: "K1,K2 IN ELBOW ARE EQUIVALENT TO 1./K1 AND 1./K2 IN BAR"
    // So for thin-walled circular, use K1=K2=0.5 -> corrected to 1/0.5 = 2.0
    Real K1 = 0.5;
    Real K2 = 0.5;
    
    // Convert angle to radians
    Real T = BETAR * M_PI / 180.0;
    
    // Compute unit vectors for local coordinate system
    // Vector from A to B
    Eigen::Vector3d AB = GPB - GPA;
    Real FLL = AB.norm();
    if (FLL == 0.0) {
        throw std::runtime_error("Zero length element");
    }
    Eigen::Vector3d SMALV0 = AB / FLL;
    
    // Reference vector (SMALLV) - already given as (0,1,0)
    Eigen::Vector3d VECJ = SMALLV;
    Real FLD = VECJ.norm();
    if (FLD == 0.0) {
        throw std::runtime_error("Zero reference vector");
    }
    VECJ /= FLD;
    
    // Compute VECK = SMALV0 × VECJ
    Eigen::Vector3d VECK;
    VECK << SMALV0(1)*VECJ(2) - SMALV0(2)*VECJ(1),
              SMALV0(2)*VECJ(0) - SMALV0(0)*VECJ(2),
              SMALV0(0)*VECJ(1) - SMALV0(1)*VECJ(0);
    
    FLL = VECK.norm();
    if (FLL == 0.0) {
        throw std::runtime_error("Cross product zero - collinear vectors");
    }
    VECK /= FLL;
    
    // Compute VECI = VECJ × VECK
    Eigen::Vector3d VECI;
    VECI << VECJ(1)*VECK(2) - VECJ(2)*VECK(1),
             VECJ(2)*VECK(0) - VECJ(0)*VECK(2),
             VECJ(0)*VECK(1) - VECJ(1)*VECK(0);
    FLL = VECI.norm();
    if (FLL == 0.0) {
        throw std::runtime_error("VECJ × VECK zero");
    }
    VECI /= FLL;
    
    // Set up intermediate variables for stiffness calculation
    if (KX < 1e-8) KX = 1.0;
    if (KY < 1e-8) KY = 1.0;
    if (KZ < 1e-8) KZ = 1.0;
    
    Real FI1 = I1 / KZ;
    Real FI2 = I2 / KY;
    Real FJK = FJ / KX;
    
    // Area factors for shear: convert to stiffness correction factors
    if (K1 < 1e-8) K1 = 1.0;
    if (K2 < 1e-8) K2 = 1.0;
    if (K1 > 1.0) K1 = 1.0 / K1;
    if (K2 > 1.0) K2 = 1.0 / K2;
    
    // Section property constants
    Real RA  = R / (A * E);
    Real RV1 = R / (2.0 * K1 * G * A);
    Real RV2 = (K1 / K2) * RV1;
    Real RT  = R / (G * FJK * 2.0);
    Real RB0 = R / (E * FI2 * 2.0);
    Real RB1 = R / (E * FI1);
    Real R2  = R * R;
    
    // Trigonometric constants
    Real ST  = std::sin(T);
    Real CT  = std::cos(T);
    Real S2T = std::sin(2.0 * T);
    Real C2T = std::cos(2.0 * T);
    
    // Build 6x6 flexibility matrix F (local coordinates at node I)
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
    F(2,3) += RT * R * (ST - T * CT);
    F(3,2) = F(2,3);
    F(2,4) += RT * R * (2.0 - 2.0 * CT - T * ST);
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
        for (int j = 0; j < i; ++j) {
            F(j,i) = F(i,j);
        }
    }
    
    // Invert flexibility matrix to get local stiffness matrix at node I
    Eigen::Matrix<Real, 6, 6> DF;
    if (!invert_6x6(F, DF)) {
        throw std::runtime_error("Flexibility matrix is singular");
    }
    
    // Build transformation matrix H
    Eigen::Matrix<Real, 6, 6> H = Eigen::Matrix<Real, 6, 6>::Zero();
    for (int i = 0; i < 6; ++i) {
        H(i,i) = -1.0;
    }
    H(3,2) = -R * (1.0 - CT); // H(4,3) in Fortran (1-indexed)
    H(4,2) = R * ST;          // H(5,3) in Fortran
    H(5,0) = -H(3,2);         // H(6,1) in Fortran
    H(5,1) = -H(4,2);         // H(6,2) in Fortran
    
    // Build 12x12 stiffness matrix S in local coordinates
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
    
    // S(1:6,1:6) = (H * DF) * H^T = temp1 * H^T
    Eigen::Matrix<Real, 6, 6> temp2 = temp1 * H.transpose();
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            S(i, j) = temp2(i,j);
        }
    }
    
    // Enforce symmetry: S = (S + S^T) / 2
    for (int i = 0; i < 12; ++i) {
        for (int j = i+1; j < 12; ++j) {
            Real avg = 0.5 * (S(i,j) + S(j,i));
            S(i,j) = avg;
            S(j,i) = avg;
        }
    }
    
    return S;
}

int main() {
    try {
        // Compute the 12x12 stiffness matrix
        Eigen::Matrix<Real, 12, 12> K = compute_elbow_stiffness();
        
        // Output as JSON
        std::cout << "{\"stiffness_matrix\":[";
        for (int i = 0; i < 12; ++i) {
            std::cout << "[";
            for (int j = 0; j < 12; ++j) {
                std::cout << format_double(K(i,j));
                if (j < 11) std::cout << ",";
            }
            std::cout << "]";
            if (i < 11) std::cout << ",";
        }
        std::cout << "]}" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}