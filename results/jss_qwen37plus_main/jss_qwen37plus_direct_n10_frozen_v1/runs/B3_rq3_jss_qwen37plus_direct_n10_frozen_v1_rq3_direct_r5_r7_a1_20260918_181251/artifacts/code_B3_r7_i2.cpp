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
    size_t epos = s.find('e');
    if (epos != std::string::npos) {
        // Remove trailing zeros after decimal point before 'e'
        size_t dotpos = s.find('.');
        if (dotpos != std::string::npos && dotpos < epos) {
            // Find last non-zero digit before 'e'
            size_t last_nonzero = epos - 1;
            while (last_nonzero > dotpos && s[last_nonzero] == '0') {
                --last_nonzero;
            }
            if (last_nonzero > dotpos && s[last_nonzero] == '.') {
                ++last_nonzero; // keep the decimal point
            }
            s.erase(last_nonzero + 1, epos - last_nonzero - 1);
        }
    }
    return s;
}

// JSON-safe string output for matrix
void print_json_stiffness(const Eigen::Matrix<double, 6, 6>& K) {
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << format_double(K(i, j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
}

int main() {
    // Test case data
    const std::vector<std::vector<double>> nodes = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {1.0, 1.5, 0.0}};
    const double E = 2.1e11;      // Pa
    const double nu = 0.3;
    const double t = 0.01;       // m

    // Extract node coordinates
    const double x1 = nodes[0][0], y1 = nodes[0][1], z1 = nodes[0][2];
    const double x2 = nodes[1][0], y2 = nodes[1][1], z2 = nodes[1][2];
    const double x3 = nodes[2][0], y3 = nodes[2][1], z3 = nodes[2][2];

    // Compute element geometry
    // Vector AB = B - A
    double abx = x2 - x1;
    double aby = y2 - y1;
    double abz = z2 - z1;
    
    // Vector AC = C - A
    double acx = x3 - x1;
    double acy = y3 - y1;
    double acz = z3 - z1;
    
    // Cross product AB × AC gives normal vector magnitude = 2 * area
    double cross_x = aby * acz - abz * acy;
    double cross_y = abz * acx - abx * acz;
    double cross_z = abx * acy - aby * acx;
    
    double area_times_2 = std::sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
    if (area_times_2 < 1.0e-12) {
        std::cerr << "Error: Degenerate element - zero area" << std::endl;
        return 1;
    }
    
    double area = area_times_2 / 2.0;
    
    // Compute shape function derivatives
    // For linear triangle: Ni = ai + bi*x + ci*y, where
    // b1 = (y2-y3), c1 = (x3-x2)
    // b2 = (y3-y1), c2 = (x1-x3)  
    // b3 = (y1-y2), c3 = (x2-x1)
    // and all divided by (2*area)
    
    double two_area = area_times_2; // 2 * area
    
    double b1 = y2 - y3;
    double c1 = x3 - x2;
    double b2 = y3 - y1;
    double c2 = x1 - x3;
    double b3 = y1 - y2;
    double c3 = x2 - x1;
    
    // Strain-displacement matrix B (3x6)
    // Row 0: [∂N1/∂x, 0, ∂N2/∂x, 0, ∂N3/∂x, 0]
    // Row 1: [0, ∂N1/∂y, 0, ∂N2/∂y, 0, ∂N3/∂y]
    // Row 2: [∂N1/∂y, ∂N1/∂x, ∂N2/∂y, ∂N2/∂x, ∂N3/∂y, ∂N3/∂x]
    
    Eigen::Matrix<double, 3, 6> B;
    
    // ∂Ni/∂x = bi / (2*area), ∂Ni/∂y = ci / (2*area)
    double inv_two_area = 1.0 / two_area;
    
    B(0, 0) = b1 * inv_two_area;  B(0, 1) = 0.0;
    B(0, 2) = b2 * inv_two_area;  B(0, 3) = 0.0;
    B(0, 4) = b3 * inv_two_area;  B(0, 5) = 0.0;
    
    B(1, 0) = 0.0;                B(1, 1) = c1 * inv_two_area;
    B(1, 2) = 0.0;                B(1, 3) = c2 * inv_two_area;
    B(1, 4) = 0.0;                B(1, 5) = c3 * inv_two_area;
    
    B(2, 0) = c1 * inv_two_area;  B(2, 1) = b1 * inv_two_area;
    B(2, 2) = c2 * inv_two_area;  B(2, 3) = b2 * inv_two_area;
    B(2, 4) = c3 * inv_two_area;  B(2, 5) = b3 * inv_two_area;
    
    // Material matrix D (3x3) for isotropic plane stress
    // D11 = E/(1-nu^2), D12 = nu*E/(1-nu^2), D22 = E/(1-nu^2), D33 = E/(2*(1+nu))
    double denom = 1.0 - nu * nu;
    double d11 = E / denom;
    double d12 = nu * E / denom;
    double d22 = E / denom;
    double d33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3> D;
    D << d11, d12, 0.0,
         d12, d22, 0.0,
         0.0, 0.0, d33;
    
    // Stiffness matrix K = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6> K = B.transpose() * D * B;
    K *= t * area;
    
    // Print as JSON
    print_json_stiffness(K);
    
    return 0;
}