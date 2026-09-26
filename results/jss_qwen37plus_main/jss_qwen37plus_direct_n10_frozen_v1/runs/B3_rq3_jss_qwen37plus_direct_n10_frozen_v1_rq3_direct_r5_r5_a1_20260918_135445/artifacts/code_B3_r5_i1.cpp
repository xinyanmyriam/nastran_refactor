#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

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
                ++last_nonzero; // Keep the decimal point if it's the last
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
    const double E_modulus = 2.1e11;   // Pa (renamed to avoid conflict with matrix E)
    const double nu = 0.3;
    const double t = 0.01;     // thickness in meters

    // Constants - define M_PI if not available
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
    const double DEGRA = M_PI / 180.0; // degrees to radians

    // Standard CST element:
    // Area = 0.5 * |(x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)|
    double area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    
    // Strain-displacement matrix B (3x6)
    // B = [b1 0 b2 0 b3 0;
    //      0 c1 0 c2 0 c3;
    //      c1 b1 c2 b2 c3 b3] / (2*area)
    // where b_i = y_{j} - y_{k}, c_i = x_{k} - x_{j} for cyclic permutations
    
    // Node numbering: 1,2,3
    // For node 1: j=2, k=3
    double b1 = y2 - y3;
    double c1 = x3 - x2;
    // For node 2: j=3, k=1
    double b2 = y3 - y1;
    double c2 = x1 - x3;
    // For node 3: j=1, k=2
    double b3 = y1 - y2;
    double c3 = x2 - x1;
    
    Eigen::Matrix<double, 3, 6> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    B /= (2.0 * area);

    // Constitutive matrix D (3x3) for plane stress
    Eigen::Matrix<double, 3, 3> D;
    double denom = 1.0 - nu*nu;
    D << E_modulus/denom, nu*E_modulus/denom, 0.0,
         nu*E_modulus/denom, E_modulus/denom, 0.0,
         0.0, 0.0, E_modulus/(2.0*(1.0+nu));

    // Stiffness matrix K = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6> K = t * area * B.transpose() * D * B;

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << double_to_scientific(K(i,j));
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}