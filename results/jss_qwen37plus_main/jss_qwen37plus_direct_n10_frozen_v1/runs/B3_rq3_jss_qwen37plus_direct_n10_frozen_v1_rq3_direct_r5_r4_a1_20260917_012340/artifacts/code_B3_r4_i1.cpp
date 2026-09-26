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
    // G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G22 = E/(1-nu^2), others zero for plane stress
    const double denom = 1.0 - nu * nu;
    const double G11 = E_modulus / denom;
    const double G12 = nu * E_modulus / denom;
    const double G22 = E_modulus / denom;
    const double G13 = 0.0, G23 = 0.0, G33 = 0.0;

    // Build G matrix (3x3) as in Fortran: [G11, G12, G13; G12, G22, G23; G13, G23, G33]
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> G;
    G << G11, G12, G13,
         G12, G22, G23,
         G13, G23, G33;

    // Compute element geometry
    // E matrix: 3x3 where columns are i, j, k vectors (in Fortran, stored row-wise)
    // First compute i-vector = (x2-x1, y2-y1, z2-z1)
    double e1 = x2 - x1;
    double e3 = y2 - y1;
    double e5 = z2 - z1;
    
    double xsubb = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    if (xsubb < 1.0e-6) {
        std::cerr << "Error: degenerate element - zero length i-vector" << std::endl;
        return 1;
    }
    
    // Normalize i-vector
    e1 /= xsubb;
    e3 /= xsubb;
    e5 /= xsubb;
    
    // j-temp = (x3-x1, y3-y1, z3-z1)
    double e2 = x3 - x1;
    double e4 = y3 - y1;
    double e6 = z3 - z1;
    
    // xsubc = i . (rC - rA)
    double xsubc = e1*e2 + e3*e4 + e5*e6;
    
    // k-vector = i cross (rC - rA)
    double e7 = e3*e6 - e5*e4;
    double e8 = e5*e2 - e1*e6;
    double e9 = e1*e4 - e3*e2;
    
    double ysubc = std::sqrt(e7*e7 + e8*e8 + e9*e9);
    if (ysubc < 1.0e-6) {
        std::cerr << "Error: degenerate element - zero length k-vector" << std::endl;
        return 1;
    }
    
    // Normalize k-vector
    e7 /= ysubc;
    e8 /= ysubc;
    e9 /= ysubc;
    
    // j-vector = k cross i
    double j1 = e5*e8 - e3*e9;
    double j3 = e1*e9 - e5*e7;
    double j5 = e3*e7 - e1*e8;
    
    double temp = std::sqrt(j1*j1 + j3*j3 + j5*j5);
    if (temp == 0.0) {
        std::cerr << "Error: degenerate element - zero length j-vector" << std::endl;
        return 1;
    }
    
    // Normalize j-vector
    j1 /= temp;
    j3 /= temp;
    j5 /= temp;
    
    // Now build E matrix (3x3) with columns: i, j, k
    // In Fortran, E is stored row-wise: E(1),E(3),E(5) = i; E(2),E(4),E(6) = j; E(7),E(8),E(9) = k
    // So in column-major storage for Eigen, we want:
    // Column 0 (i): [e1, e3, e5]
    // Column 1 (j): [j1, j3, j5]
    // Column 2 (k): [e7, e8, e9]
    Eigen::Matrix<double, 3, 3, Eigen::ColMajor> E_matrix;
    E_matrix.col(0) << e1, e3, e5;
    E_matrix.col(1) << j1, j3, j5;
    E_matrix.col(2) << e7, e8, e9;
    
    // Volume = xsubb * ysubc * t / 2.0
    double vol = xsubb * ysubc * t / 2.0;
    
    // reelmU = 1/xsubb, flambda = 1/ysubc, delta = xsubc/xsubb - 1.0
    double reelmU = 1.0 / xsubb;
    double flambda = 1.0 / ysubc;
    double delta = xsubc / xsubb - 1.0;
    
    // Build C matrix (3x6) partitioned as CSUBA (3x2), CSUBB (3x2), CSUBC (3x2)
    // Stored row-wise in Fortran: C(1..18)
    // We'll build a 3x6 matrix in row-major order
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C;
    
    // CSUBA (first 2 columns)
    C(0,0) = -reelmU;      C(0,1) = 0.0;
    C(1,0) = 0.0;          C(1,1) = flambda * delta;
    C(2,0) = 0.0;          C(2,1) = -reelmU;
    
    // CSUBB (next 2 columns)
    C(0,2) = reelmU;       C(0,3) = 0.0;
    C(1,2) = 0.0;          C(1,3) = -flambda * reelmU * xsubc;
    C(2,2) = 0.0;          C(2,3) = reelmU;
    
    // CSUBC (last 2 columns)
    C(0,4) = 0.0;          C(0,5) = 0.0;
    C(1,4) = 0.0;          C(1,5) = flambda;
    C(2,4) = 0.0;          C(2,5) = flambda;
    
    // Since the test case has no coordinate system transformations (all CSID = 0),
    // we skip the transformation matrices TI and TJ (they are identity)
    // So we can directly compute: K = vol * (E * C)^T * G * (E * C)
    
    // First compute E * C (3x6)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> EC = E_matrix * C;
    
    // Then compute (E*C)^T * G * (E*C) which is 6x6
    // Let M = E*C, then result = M^T * G * M
    Eigen::Matrix<double, 6, 6> K = EC.transpose() * G * EC;
    
    // Scale by volume
    K *= vol;
    
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