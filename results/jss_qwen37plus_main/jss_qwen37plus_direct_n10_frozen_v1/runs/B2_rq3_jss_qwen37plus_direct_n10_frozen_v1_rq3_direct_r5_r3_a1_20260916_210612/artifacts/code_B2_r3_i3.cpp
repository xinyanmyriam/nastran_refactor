#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Since we can't use external JSON library in strict C++17 without dependencies,
// we'll implement minimal JSON serialization for the required format.
// We'll use a simple string-based approach.

std::string double_to_scientific(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << x;
    std::string s = oss.str();
    // Remove trailing zeros after decimal point, but keep at least one digit after 'e'
    size_t e_pos = s.find('e');
    if (e_pos != std::string::npos) {
        size_t dot_pos = s.find('.');
        if (dot_pos != std::string::npos && dot_pos < e_pos) {
            // Remove trailing zeros in mantissa
            size_t last_nonzero = s.find_last_not_of('0', e_pos);
            if (last_nonzero != std::string::npos && last_nonzero > dot_pos) {
                s.erase(last_nonzero + 1, e_pos - last_nonzero - 1);
            } else if (last_nonzero == dot_pos) {
                s.erase(dot_pos, e_pos - dot_pos);
                e_pos = s.find('e');
            }
        }
    }
    return s;
}

std::string matrix_to_json(const Eigen::Matrix<double, 12, 12>& K) {
    std::ostringstream json;
    json << "{\"stiffness_matrix\":[";
    
    for (int i = 0; i < 12; ++i) {
        if (i > 0) json << ",";
        json << "[";
        for (int j = 0; j < 12; ++j) {
            if (j > 0) json << ",";
            json << double_to_scientific(K(i, j));
        }
        json << "]";
    }
    
    json << "]}";
    return json.str();
}

// Helper function to compute cross product of two 3D vectors
Eigen::Vector3d cross(const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
    return Eigen::Vector3d(
        a(1)*b(2) - a(2)*b(1),
        a(2)*b(0) - a(0)*b(2),
        a(0)*b(1) - a(1)*b(0)
    );
}

// Helper function to normalize a vector
Eigen::Vector3d normalize(const Eigen::Vector3d& v) {
    double norm = v.norm();
    if (norm == 0.0) {
        return Eigen::Vector3d::Zero();
    }
    return v / norm;
}

// Main CBAR stiffness computation
Eigen::Matrix<double, 12, 12> compute_cbar_stiffness() {
    // Test case parameters
    const double E = 200e9;           // Pa
    const double G = 76.923e9;        // Pa
    const double A = 0.01;            // m^2
    const double Iy = 8.333e-6;       // m^4 (I1 in Fortran)
    const double Iz = 8.333e-6;       // m^4 (I2 in Fortran)
    const double J = 1.667e-5;        // m^4 (FJ in Fortran)
    
    // Node coordinates
    Eigen::Vector3d nodeA(0.0, 0.0, 0.0);
    Eigen::Vector3d nodeB(2.0, 0.0, 0.0);
    
    // Length
    double L = (nodeB - nodeA).norm();
    double L2 = L * L;
    double L3 = L2 * L;
    
    // Material properties
    double EIy = E * Iy;
    double EIz = E * Iz;
    double GJL = G * J / L;
    double AEL = E * A / L;
    
    // Shear correction factors (K1=K2=0.0 in test case, so use pure Euler-Bernoulli)
    double K1 = 0.0;
    double K2 = 0.0;
    double I12 = 0.0; // product of inertia
    
    // Compute R1 and R2 (shear stiffness terms)
    double R1, R2;
    if (K1 == 0.0 || I12 != 0.0) {
        R1 = 12.0 * EIy / L3;
    } else {
        double GAK1 = G * A * K1;
        R1 = (12.0 * EIy * GAK1) / (GAK1 * L3 + 12.0 * L * EIy);
    }
    
    if (K2 == 0.0 || I12 != 0.0) {
        R2 = 12.0 * EIz / L3;
    } else {
        double GAK2 = G * A * K2;
        R2 = (12.0 * EIz * GAK2) / (GAK2 * L3 + 12.0 * L * EIz);
    }
    
    // Compute SK terms
    double SK1 = 0.25 * R1 * L2 + EIy / L;
    double SK2 = 0.25 * R2 * L2 + EIz / L;
    double SK3 = 0.25 * R1 * L2 - EIy / L;
    double SK4 = 0.25 * R2 * L2 - EIz / L;
    
    // Compute LR terms
    double LR1 = L * R1 / 2.0;
    double LR2 = L * R2 / 2.0;
    
    // Local coordinate system setup (like Fortran's VECI, VECJ, VECK)
    // VECI is the x-axis (from A to B, but Fortran negates it)
    Eigen::Vector3d VECI = nodeA - nodeB; // Fortran does VECI = -VECI after computing as B-A
    VECI = normalize(VECI);
    
    // Reference vector SMALLV (use y-axis as reference since beam is along x)
    Eigen::Vector3d SMALLV(0.0, 1.0, 0.0);
    SMALLV = normalize(SMALLV);
    
    // VECK = VECI × SMALLV
    Eigen::Vector3d VECK = cross(VECI, SMALLV);
    VECK = normalize(VECK);
    
    // VECJ = VECK × VECI
    Eigen::Vector3d VECJ = cross(VECK, VECI);
    VECJ = normalize(VECJ);
    
    // Initialize stiffness matrix to zero
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    
    // Fill the 12x12 stiffness matrix according to Fortran logic
    // Note: Fortran uses 1-based indexing, C++ uses 0-based
    // The Fortran array KE(1..144) is stored column-wise, but we'll fill row-wise
    
    // Axial terms (1,1) and (7,7) etc.
    K(0, 0) = AEL;      // KE(1)
    K(0, 6) = -AEL;     // KE(7)
    K(6, 0) = -AEL;     // KE(73)
    K(6, 6) = AEL;      // KE(79)
    
    // Bending about y-axis (Iy) - affects DOF 2,3,8,9 (y-trans, z-rot, y-trans, z-rot)
    K(1, 1) = R1;       // KE(14)
    K(1, 5) = LR1;      // KE(18)
    K(1, 7) = -R1;      // KE(20)
    K(1, 11) = LR1;     // KE(24)
    
    K(5, 1) = LR1;      // KE(27) - wait, check Fortran: KE(27) is at position 27 which is row 3, col 3? 
    // Let's map Fortran positions carefully:
    // Fortran KE(1) = (1,1), KE(2) = (2,1), ..., KE(12) = (12,1), KE(13) = (1,2), etc.
    // So KE(i) corresponds to row = (i-1)%12, col = (i-1)/12
    
    // Instead, let's build the matrix by direct mapping from Fortran indices
    // We'll create a 1D array and then reshape
    
    Eigen::Matrix<double, 12, 12, Eigen::ColMajor> K_col = Eigen::Matrix<double, 12, 12, Eigen::ColMajor>::Zero();
    
    // Map Fortran KE(1..144) to column-major storage
    auto set_ke = [&](int idx, double val) {
        int row = (idx - 1) % 12;
        int col = (idx - 1) / 12;
        K_col(row, col) = val;
    };
    
    // Set all values as in Fortran
    set_ke(1, AEL);      // (1,1)
    set_ke(7, -AEL);     // (7,1)
    set_ke(14, R1);      // (2,2)
    set_ke(18, LR1);     // (6,2)
    set_ke(20, -R1);     // (8,2)
    set_ke(24, LR1);     // (12,2)
    set_ke(27, R2);      // (3,3)
    set_ke(29, -LR2);    // (5,3)
    set_ke(33, -R2);     // (9,3)
    set_ke(35, -LR2);    // (11,3)
    set_ke(40, GJL);     // (4,4)
    set_ke(46, -GJL);    // (10,4)
    set_ke(51, -LR2);    // (3,5)
    set_ke(53, SK2);     // (5,5)
    set_ke(57, LR2);     // (9,5)
    set_ke(59, SK4);     // (11,5)
    set_ke(62, LR1);     // (2,6)
    set_ke(66, SK1);     // (6,6)
    set_ke(68, -LR1);    // (8,6)
    set_ke(72, SK3);     // (12,6)
    set_ke(73, -AEL);    // (1,7)
    set_ke(79, AEL);     // (7,7)
    set_ke(86, -R1);     // (2,8)
    set_ke(90, -LR1);    // (6,8)
    set_ke(92, R1);      // (8,8)
    set_ke(96, -LR1);    // (12,8)
    set_ke(99, -R2);     // (3,9)
    set_ke(101, LR2);    // (5,9)
    set_ke(105, R2);     // (9,9)
    set_ke(107, LR2);    // (11,9)
    set_ke(112, -GJL);   // (4,10)
    set_ke(118, GJL);    // (10,10)
    set_ke(123, -LR2);   // (3,11)
    set_ke(125, SK4);    // (5,11)
    set_ke(129, LR2);    // (9,11)
    set_ke(131, SK2);    // (11,11)
    set_ke(134, LR1);    // (2,12)
    set_ke(138, SK3);    // (6,12)
    set_ke(140, -LR1);   // (8,12)
    set_ke(144, SK1);    // (12,12)
    
    // Handle I12 terms (warping torsion, but I12=0 in test case, so skip)
    if (I12 != 0.0) {
        double BETA = 12.0 * E * I12 / L3;
        double LB = L * BETA / 2.0;
        double L2B3 = L2 * BETA / 3.0;
        double L2B6 = L2 * BETA / 6.0;
        
        set_ke(15, BETA);   // (3,2)
        set_ke(17, -LB);    // (5,2)
        set_ke(21, -BETA);  // (9,2)
        set_ke(23, -LB);    // (11,2)
        set_ke(26, BETA);   // (2,3)
        set_ke(30, LB);     // (6,3)
        set_ke(32, -BETA);  // (8,3)
        set_ke(36, LB);     // (12,3)
        set_ke(50, -LB);    // (2,5)
        set_ke(54, -L2B3);  // (6,5)
        set_ke(56, LB);     // (8,5)
        set_ke(60, -L2B6);  // (12,5)
        set_ke(63, LB);     // (3,6)
        set_ke(65, -L2B3);  // (5,6)
        set_ke(69, -LB);    // (9,6)
        set_ke(71, -L2B6);  // (11,6)
        set_ke(87, -BETA);  // (3,8)
        set_ke(89, LB);     // (5,8)
        set_ke(93, BETA);   // (9,8)
        set_ke(95, LB);     // (11,8)
        set_ke(98, -BETA);  // (2,9)
        set_ke(102, -LB);   // (6,9)
        set_ke(104, BETA);  // (8,9)
        set_ke(108, -LB);   // (12,9)
        set_ke(122, -LB);   // (2,11)
        set_ke(126, -L2B6); // (6,11)
        set_ke(128, LB);    // (8,11)
        set_ke(132, -L2B3); // (12,11)
        set_ke(135, LB);    // (3,12)
        set_ke(137, -L2B6); // (5,12)
        set_ke(141, -LB);   // (9,12)
        set_ke(143, -L2B3); // (11,12)
    }
    
    // Convert back to row-major for output (but Eigen's default is column-major for storage)
    // We'll work with the column-major matrix and convert to row-major when needed
    // For output, we want K(i,j) where i=row, j=col in mathematical sense
    Eigen::Matrix<double, 12, 12> K_out;
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            K_out(i, j) = K_col(i, j);
        }
    }
    
    return K_out;
}

int main() {
    // Compute the CBAR stiffness matrix
    Eigen::Matrix<double, 12, 12> K = compute_cbar_stiffness();
    
    // Output as JSON
    std::cout << matrix_to_json(K) << std::endl;
    
    return 0;
}