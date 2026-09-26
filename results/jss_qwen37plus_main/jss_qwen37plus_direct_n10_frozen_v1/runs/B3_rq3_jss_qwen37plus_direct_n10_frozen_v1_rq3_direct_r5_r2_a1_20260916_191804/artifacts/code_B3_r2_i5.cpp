#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// For JSON output formatting
std::string double_to_scientific(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos) {
        // Remove trailing zeros after decimal
        size_t dot = s.find('.');
        if (dot != std::string::npos && dot < epos) {
            while (s.back() == '0' && s.size() > dot + 2) {
                s.pop_back();
            }
            if (s.back() == '.') {
                s.pop_back();
            }
        }
        // Remove '+' from exponent if present
        size_t plus = s.find('+', epos);
        if (plus != std::string::npos) {
            s.erase(plus, 1);
        }
    }
    return s;
}

int main() {
    // Test case data
    const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    const double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    const double x3 = 1.0, y3 = 1.5, z3 = 0.0;
    
    const double E = 2.1e11;  // Pa
    const double nu = 0.3;
    const double t = 0.01;    // m
    
    // Compute signed area of triangle
    double area_signed = 0.5 * ((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    double area = std::abs(area_signed);
    if (std::abs(area) < 1.0e-12) {
        std::cerr << "Error: degenerate element - zero area" << std::endl;
        return 1;
    }
    
    // For CST, B matrix is 3x6:
    // B = [ b1 0  b2 0  b3 0 ]
    //     [ 0  c1  0 c2  0 c3 ]
    //     [ c1 b1  c2 b2  c3 b3 ]
    // where for node i: bi = (yj - yk)/2A, ci = (xk - xj)/2A (with j,k cyclic)
    // Use absolute area for B-matrix coefficients to ensure proper scaling
    double twoA = 2.0 * area;
    double b1 = (y2 - y3) / twoA;
    double b2 = (y3 - y1) / twoA;
    double b3 = (y1 - y2) / twoA;
    // Fixed: use (xk - xj)/2A for ci, where (i,j,k) is cyclic order
    // For node 1: j=2, k=3 -> c1 = (x3 - x2)/2A
    // For node 2: j=3, k=1 -> c2 = (x1 - x3)/2A  
    // For node 3: j=1, k=2 -> c3 = (x2 - x1)/2A
    double c1 = (x3 - x2) / twoA;
    double c2 = (x1 - x3) / twoA;
    double c3 = (x2 - x1) / twoA;
    
    // Build B matrix (3x6) - using RowMajor for consistency with Eigen's default
    Eigen::Matrix<double, 3, 6> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    
    // D matrix for plane stress (3x3)
    Eigen::Matrix<double, 3, 3> D;
    double denom = 1.0 - nu * nu;
    D << E/denom, nu*E/denom, 0.0,
         nu*E/denom, E/denom, 0.0,
         0.0, 0.0, E/(2.0*(1.0+nu));
    
    // K = t * |area| * B^T * D * B
    // B.transpose() is 6x3, D is 3x3, B is 3x6
    // So B.transpose() * D is 6x3, then (B.transpose() * D) * B is 6x6
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