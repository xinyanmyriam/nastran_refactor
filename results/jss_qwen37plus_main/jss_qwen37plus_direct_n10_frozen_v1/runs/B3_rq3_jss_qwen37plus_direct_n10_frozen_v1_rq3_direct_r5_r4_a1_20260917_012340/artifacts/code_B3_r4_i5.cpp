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
    oss << std::scientific << std::setprecision(12) << x;
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
            if (s[last_nonzero] == '.') {
                ++last_nonzero;
            }
            s.erase(last_nonzero + 1, e_pos - last_nonzero - 1);
            // Ensure at least one digit after decimal
            if (s[dot_pos + 1] == 'e') {
                s.insert(dot_pos + 1, "0");
            }
        }
        // Fix exponent sign: ensure + is present for positive exponents
        if (e_pos + 1 < s.length() && s[e_pos + 1] == '-') {
            // already has minus, ok
        } else if (e_pos + 1 < s.length()) {
            // has digit, need to insert +
            s.insert(e_pos + 1, "+");
        }
    }
    return s;
}

int main() {
    // Test case parameters
    const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    const double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    const double x3 = 1.0, y3 = 1.5, z3 = 0.0;
    
    const double E_modulus = 2.1e11;  // Pa
    const double nu = 0.3;
    const double t = 0.01;    // m

    // Material properties: isotropic plane stress
    // G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G22 = E/(1-nu^2)
    // G33 = E/(2*(1+nu)) for shear modulus
    const double denom = 1.0 - nu * nu;
    const double G11 = E_modulus / denom;
    const double G12 = nu * E_modulus / denom;
    const double G22 = E_modulus / denom;
    const double G33 = E_modulus / (2.0 * (1.0 + nu)); // Fixed: shear modulus
    const double G13 = 0.0, G23 = 0.0;

    // Build G matrix (3x3) as in Fortran: [G11, G12, G13; G12, G22, G23; G13, G23, G33]
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> G;
    G << G11, G12, G13,
         G12, G22, G23,
         G13, G23, G33;

    // Compute element geometry - area of triangle
    double area = 0.5 * std::abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
    if (area < 1.0e-12) {
        std::cerr << "Error: degenerate element - zero area" << std::endl;
        return 1;
    }

    // Standard CTRIA3 B-matrix coefficients (strain-displacement)
    // For cyclic ordering (1,2,3): 
    // beta1 = y2 - y3, gamma1 = x3 - x2
    // beta2 = y3 - y1, gamma2 = x1 - x3
    // beta3 = y1 - y2, gamma3 = x2 - x1
    double beta1 = y2 - y3;
    double gamma1 = x3 - x2;  // Fixed: was x2 - x3
    double beta2 = y3 - y1;
    double gamma2 = x1 - x3;  // Fixed: was x3 - x1
    double beta3 = y1 - y2;
    double gamma3 = x2 - x1;  // Fixed: was x1 - x2

    // Build C matrix (3x6) where C = 2*area*B
    // Row 0: [beta1, 0, beta2, 0, beta3, 0] for epsilon_x
    // Row 1: [0, gamma1, 0, gamma2, 0, gamma3] for epsilon_y  
    // Row 2: [gamma1, beta1, gamma2, beta2, gamma3, beta3] for gamma_xy
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C;
    C << beta1, 0.0, beta2, 0.0, beta3, 0.0,
         0.0, gamma1, 0.0, gamma2, 0.0, gamma3,
         gamma1, beta1, gamma2, beta2, gamma3, beta3;

    // Volume = area * thickness
    double vol = area * t;

    // For CTRIA3, the strain-displacement matrix B = C / (2*area)
    // So stiffness K = vol * B^T * G * B = vol * (C/(2*area))^T * G * (C/(2*area))
    // = vol * C^T * G * C / (4 * area^2) = t * area * C^T * G * C / (4 * area^2)
    // = t * C^T * G * C / (4 * area)
    double scale_factor = t / (4.0 * area);

    // Compute K = scale_factor * C^T * G * C
    Eigen::Matrix<double, 6, 6> K = scale_factor * C.transpose() * G * C;

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