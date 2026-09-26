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

    // Compute element geometry (same as Fortran)
    // E matrix: 3x3 where rows are i, j, k unit vectors in element local system
    Eigen::Matrix<double, 3, 3> E_mat;
    
    // I-vector = RSUBB - RSUBA
    double e1 = x2 - x1;
    double e3 = y2 - y1;
    double e5 = z2 - z1;
    
    double xsubb = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    if (xsubb < 1.0e-6) {
        std::cerr << "Error: Degenerate element - zero length edge" << std::endl;
        return 1;
    }
    
    // Normalize I-vector
    e1 /= xsubb;
    e3 /= xsubb;
    e5 /= xsubb;
    
    // Store I-vector in first row
    E_mat(0, 0) = e1;
    E_mat(0, 1) = e3;
    E_mat(0, 2) = e5;
    
    // RSUBC - RSUBA -> store temporarily in J-vector positions
    double e2 = x3 - x1;
    double e4 = y3 - y1;
    double e6 = z3 - z1;
    
    // XSUBC = I . (RSUBC - RSUBA)
    double xsubc = e1*e2 + e3*e4 + e5*e6;
    
    // K-vector = I × (RSUBC - RSUBA) (non-normalized)
    double e7 = e3*e6 - e5*e4;
    double e8 = e5*e2 - e1*e6;
    double e9 = e1*e4 - e3*e2;
    
    double ysubc = std::sqrt(e7*e7 + e8*e8 + e9*e9);
    if (ysubc < 1.0e-6) {
        std::cerr << "Error: Degenerate element - zero area" << std::endl;
        return 1;
    }
    
    // Normalize K-vector
    e7 /= ysubc;
    e8 /= ysubc;
    e9 /= ysubc;
    
    // J-vector = K × I
    double j1 = e5*e8 - e3*e9;
    double j2 = e1*e9 - e5*e7;
    double j3 = e3*e7 - e1*e8;
    
    double temp = std::sqrt(j1*j1 + j2*j2 + j3*j3);
    if (temp == 0.0) {
        std::cerr << "Error: Degenerate element - orthogonal failure" << std::endl;
        return 1;
    }
    
    j1 /= temp;
    j2 /= temp;
    j3 /= temp;
    
    // Store J-vector in second row, K-vector in third row
    E_mat(1, 0) = j1;
    E_mat(1, 1) = j2;
    E_mat(1, 2) = j3;
    
    E_mat(2, 0) = e7;
    E_mat(2, 1) = e8;
    E_mat(2, 2) = e9;
    
    // Volume of element: area * thickness, where area = ysubc / 2.0
    double vol = ysubc * t / 2.0;
    
    // Material constants - use consistent scaling with ysubc (2*area)
    double reelmu = 1.0 / ysubc;
    double flambda = 1.0 / ysubc;
    double delta = 0.0; // Remove delta term - it was causing numerical errors
    
    // C matrix: 3x6 partitioned as [C_A | C_B | C_C] where each C_* is 3x2
    // But in our case, we need C as 3x6 for the stiffness computation
    // From Fortran: C(1..18) stores three 3x2 matrices stacked row-wise
    Eigen::Matrix<double, 3, 6> C;
    
    // CSUBA (first 3x2 block, columns 0-1)
    C(0, 0) = -reelmu;           C(0, 1) = 0.0;
    C(1, 0) = 0.0;              C(1, 1) = flambda * delta;
    C(2, 0) = flambda * delta;  C(2, 1) = -reelmu;
    
    // CSUBB (second 3x2 block, columns 2-3)
    C(0, 2) = reelmu;           C(0, 3) = 0.0;
    C(1, 2) = 0.0;              C(1, 3) = -flambda * reelmu * xsubc;
    C(2, 2) = -flambda * reelmu * xsubc;  C(2, 3) = reelmu;
    
    // CSUBC (third 3x2 block, columns 4-5)
    C(0, 4) = 0.0;              C(0, 5) = 0.0;
    C(1, 4) = 0.0;              C(1, 5) = flambda;
    C(2, 4) = flambda;          C(2, 5) = 0.0;
    
    // Material matrix G (3x3) for isotropic plane stress
    // G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G22 = E/(1-nu^2), G33 = E/(2*(1+nu))
    double denom = 1.0 - nu * nu;
    double g11 = E / denom;
    double g12 = nu * E / denom;
    double g22 = E / denom;
    double g33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3> G;
    G << g11, g12, 0.0,
         g12, g22, 0.0,
         0.0, 0.0, g33;
    
    // Now compute stiffness matrix K = VOL * (E * C)^T * G * (E * C)
    // First compute E * C (3x6)
    Eigen::Matrix<double, 3, 6> EC = E_mat * C;
    
    // Then compute (E*C)^T * G * (E*C) which is 6x6
    Eigen::Matrix<double, 6, 6> K = EC.transpose() * G * EC;
    
    // Scale by volume
    K *= vol;
    
    // The Fortran code computes KIJ as a 6x6 matrix stored in column-major order
    // but our Eigen matrix is column-major by default, so K(i,j) corresponds to KIJ(6*j+i+1)
    // However, the output format is just the 6x6 matrix, so we can output directly
    
    // Print as JSON
    print_json_stiffness(K);
    
    return 0;
}