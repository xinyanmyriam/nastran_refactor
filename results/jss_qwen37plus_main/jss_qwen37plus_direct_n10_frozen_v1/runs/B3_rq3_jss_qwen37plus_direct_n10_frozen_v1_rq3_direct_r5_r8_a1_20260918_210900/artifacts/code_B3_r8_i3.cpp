#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Ensure M_PI is defined (not standard in C++)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// For JSON output formatting
std::string double_to_scientific(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t e_pos = s.find('e');
    if (e_pos != std::string::npos) {
        // Remove trailing zeros after decimal
        size_t dot_pos = s.find('.');
        if (dot_pos != std::string::npos && dot_pos < e_pos) {
            // Remove zeros from end of mantissa
            size_t last_nonzero = e_pos - 1;
            while (last_nonzero > dot_pos && s[last_nonzero] == '0') {
                --last_nonzero;
            }
            if (last_nonzero > dot_pos && s[last_nonzero] == '.') {
                ++last_nonzero; // keep the dot if it's the last
            }
            s.erase(last_nonzero + 1, e_pos - last_nonzero - 1);
            // Reconstruct
            std::string mantissa = s.substr(0, e_pos);
            std::string exp = s.substr(e_pos);
            // Remove '+' from exponent if present
            if (exp.length() >= 2 && exp[1] == '+') {
                exp = exp.substr(0, 1) + exp.substr(2);
            }
            s = mantissa + exp;
        }
    }
    return s;
}

int main() {
    // Test case: Nodes at (0,0), (2,0), (1,1.5) meters
    const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    const double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    const double x3 = 1.0, y3 = 1.5, z3 = 0.0;

    // Material properties
    const double E = 2.1e11;   // Pa
    const double nu = 0.3;
    const double t = 0.01;      // thickness in meters

    // Compute element geometry
    // E matrix: 3x3 for coordinate system construction
    Eigen::Matrix<double, 3, 3, Eigen::ColMajor> E_mat;
    
    // First vector (I-vector) = node2 - node1
    double e1 = x2 - x1;
    double e3 = y2 - y1;
    double e5 = z2 - z1;
    
    double xsubb = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    if (xsubb < 1.0e-6) {
        std::cerr << "Error: Degenerate element - zero length I-vector" << std::endl;
        return 1;
    }
    
    // Normalize I-vector
    e1 /= xsubb;
    e3 /= xsubb;
    e5 /= xsubb;
    
    // Second temporary vector = node3 - node1
    double e2 = x3 - x1;
    double e4 = y3 - y1;
    double e6 = z3 - z1;
    
    // XSUBC = I . (node3 - node1)
    double xsubc = e1*e2 + e3*e4 + e5*e6;
    
    // K-vector = I cross (node3 - node1)
    double e7 = e3*e6 - e5*e4;
    double e8 = e5*e2 - e1*e6;
    double e9 = e1*e4 - e3*e2;
    
    double ysubc = std::sqrt(e7*e7 + e8*e8 + e9*e9);
    if (ysubc < 1.0e-6) {
        std::cerr << "Error: Degenerate element - zero length K-vector" << std::endl;
        return 1;
    }
    
    // Normalize K-vector
    e7 /= ysubc;
    e8 /= ysubc;
    e9 /= ysubc;
    
    // J-vector = K cross I
    double j1 = e5*e8 - e3*e9;
    double j2 = e1*e9 - e5*e7;
    double j3 = e3*e7 - e1*e8;
    
    double temp = std::sqrt(j1*j1 + j2*j2 + j3*j3);
    if (temp == 0.0) {
        std::cerr << "Error: Degenerate element - zero length J-vector" << std::endl;
        return 1;
    }
    
    j1 /= temp;
    j2 /= temp;
    j3 /= temp;
    
    // Fill E matrix: columns are I, J, K vectors
    E_mat << e1, j1, e7,
             e3, j2, e8,
             e5, j3, e9;

    // Volume of element
    double vol = xsubb * ysubc * t / 2.0;
    
    // Material constants for isotropic plane stress
    double reelmu = 1.0 / xsubb;
    double flambda = 1.0 / ysubc;
    double delta = xsubc / xsubb - 1.0;
    
    // C matrix: 3x6, partitioned as [C_A | C_B | C_C] where each C_X is 3x2
    Eigen::Matrix<double, 3, 6, Eigen::ColMajor> C_mat;
    
    // C_A (columns 0-1)
    C_mat(0,0) = -reelmu; C_mat(0,1) = 0.0;
    C_mat(1,0) = 0.0;     C_mat(1,1) = flambda * delta;
    C_mat(2,0) = flambda * delta; C_mat(2,1) = -reelmu;
    
    // C_B (columns 2-3)
    C_mat(0,2) = reelmu; C_mat(0,3) = 0.0;
    C_mat(1,2) = 0.0;      C_mat(1,3) = -flambda * reelmu * xsubc;
    C_mat(2,2) = -flambda * reelmu * xsubc; C_mat(2,3) = reelmu;
    
    // C_C (columns 4-5)
    C_mat(0,4) = 0.0; C_mat(0,5) = 0.0;
    C_mat(1,4) = 0.0; C_mat(1,5) = flambda;
    C_mat(2,4) = flambda; C_mat(2,5) = 0.0;
    
    // Material matrix G for isotropic plane stress
    // G = [G11 G12 0; G12 G22 0; 0 0 G33]
    // For isotropic: G11 = G22 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G33 = E/(2*(1+nu))
    double denom = 1.0 - nu*nu;
    double g11 = E / denom;
    double g12 = nu * E / denom;
    double g33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3, Eigen::ColMajor> G_mat;
    G_mat << g11, g12, 0.0,
             g12, g11, 0.0,
             0.0, 0.0, g33;
    
    // Compute stiffness matrix: K = vol * (E * C)^T * G * (E * C)
    // First compute E * C (3x6 matrix)
    Eigen::Matrix<double, 3, 6, Eigen::ColMajor> EC_mat = E_mat * C_mat;
    
    // Then compute (E*C)^T * G * (E*C) which is 6x6
    Eigen::Matrix<double, 6, 6, Eigen::ColMajor> K_temp = EC_mat.transpose() * G_mat * EC_mat;
    
    // Multiply by volume
    Eigen::Matrix<double, 6, 6, Eigen::ColMajor> K_mat = vol * K_temp;
    
    // Since this is a membrane element with 2 DOF per node (ux, uy),
    // the ordering is: node1_ux, node1_uy, node2_ux, node2_uy, node3_ux, node3_uy
    
    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << double_to_scientific(K_mat(i,j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}