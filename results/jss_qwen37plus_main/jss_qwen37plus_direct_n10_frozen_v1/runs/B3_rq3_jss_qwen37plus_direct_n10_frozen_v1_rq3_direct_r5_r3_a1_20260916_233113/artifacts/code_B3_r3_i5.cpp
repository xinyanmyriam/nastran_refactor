#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Fallback M_PI definition if not provided
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// For JSON output formatting
std::string double_to_scientific(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos) {
        // Remove trailing zeros after decimal
        size_t dot = s.find('.');
        if (dot != std::string::npos && dot < epos) {
            while (s.length() > dot + 1 && s.back() == '0') {
                s.pop_back();
            }
            if (s.back() == '.') {
                s.pop_back();
            }
        }
        // Ensure exponent has sign
        std::string exp = s.substr(epos+1);
        if (exp.length() > 0 && exp[0] != '+' && exp[0] != '-') {
            s.replace(epos+1, 0, "+");
        }
    }
    return s;
}

int main() {
    // Test case data
    const std::vector<std::vector<double>> nodes = {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {1.0, 1.5, 0.0}};
    const double E = 2.1e11;  // Pa
    const double nu = 0.3;
    const double t = 0.01;   // m

    // Extract node coordinates
    const double X1 = nodes[0][0], Y1 = nodes[0][1], Z1 = nodes[0][2];
    const double X2 = nodes[1][0], Y2 = nodes[1][1], Z2 = nodes[1][2];
    const double X3 = nodes[2][0], Y3 = nodes[2][1], Z3 = nodes[2][2];

    // Compute element geometry
    // E matrix: 3x3, but only first two columns used for membrane (I, J, K vectors)
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> E_mat;
    E_mat.setZero();

    // I-vector = RSUBB - RSUBA
    E_mat(0,0) = X2 - X1;
    E_mat(1,0) = Y2 - Y1;
    E_mat(2,0) = Z2 - Z1;

    // XSUBB = length of I-vector
    double XSUBB = std::sqrt(E_mat(0,0)*E_mat(0,0) + E_mat(1,0)*E_mat(1,0) + E_mat(2,0)*E_mat(2,0));
    if (XSUBB < 1.0e-6) {
        std::cerr << "Error: Zero length I-vector" << std::endl;
        return 1;
    }

    // Normalize I-vector
    E_mat(0,0) /= XSUBB;
    E_mat(1,0) /= XSUBB;
    E_mat(2,0) /= XSUBB;

    // Store RSUBC - RSUBA temporarily in column 1 (J-vector location)
    E_mat(0,1) = X3 - X1;
    E_mat(1,1) = Y3 - Y1;
    E_mat(2,1) = Z3 - Z1;

    // XSUBC = I . (RSUBC - RSUBA)
    double XSUBC = E_mat(0,0)*E_mat(0,1) + E_mat(1,0)*E_mat(1,1) + E_mat(2,0)*E_mat(2,1);

    // K-vector = I cross (RSUBC - RSUBA) [non-normalized]
    E_mat(0,2) = E_mat(1,0)*E_mat(2,1) - E_mat(2,0)*E_mat(1,1);
    E_mat(1,2) = E_mat(2,0)*E_mat(0,1) - E_mat(0,0)*E_mat(2,1);
    E_mat(2,2) = E_mat(0,0)*E_mat(1,1) - E_mat(1,0)*E_mat(0,1);

    // YSUBC = length of K-vector
    double YSUBC = std::sqrt(E_mat(0,2)*E_mat(0,2) + E_mat(1,2)*E_mat(1,2) + E_mat(2,2)*E_mat(2,2));
    if (YSUBC < 1.0e-6) {
        std::cerr << "Error: Zero length K-vector" << std::endl;
        return 1;
    }

    // Normalize K-vector
    E_mat(0,2) /= YSUBC;
    E_mat(1,2) /= YSUBC;
    E_mat(2,2) /= YSUBC;

    // J-vector = K cross I
    E_mat(0,1) = E_mat(1,2)*E_mat(2,0) - E_mat(2,2)*E_mat(1,0);
    E_mat(1,1) = E_mat(2,2)*E_mat(0,0) - E_mat(0,2)*E_mat(2,0);
    E_mat(2,1) = E_mat(0,2)*E_mat(1,0) - E_mat(1,2)*E_mat(0,0);

    // Normalize J-vector
    double TEMP = std::sqrt(E_mat(0,1)*E_mat(0,1) + E_mat(1,1)*E_mat(1,1) + E_mat(2,1)*E_mat(2,1));
    if (TEMP == 0.0) {
        std::cerr << "Error: Zero length J-vector" << std::endl;
        return 1;
    }
    E_mat(0,1) /= TEMP;
    E_mat(1,1) /= TEMP;
    E_mat(2,1) /= TEMP;

    // Volume, REELMU, FLAMDA, DELTA
    double VOL = XSUBB * YSUBC * t / 2.0;
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB - 1.0;

    // Build C matrix (3x6), partitioned as CSUBA, CSUBB, CSUBC (each 3x2)
    // Stored row-wise: C(1..18) = [CSUBA(1,1), CSUBA(1,2), CSUBA(2,1), CSUBA(2,2), CSUBA(3,1), CSUBA(3,2),
    //                               CSUBB(1,1), ... , CSUBC(3,2)]
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C_mat;
    C_mat.setZero();

    // CSUBA (first 2 columns)
    C_mat(0,0) = -REELMU;           // C(1)
    C_mat(0,1) = 0.0;              // C(2)
    C_mat(1,0) = 0.0;              // C(3)
    C_mat(1,1) = FLAMDA * DELTA;   // C(4)
    C_mat(2,0) = FLAMDA * DELTA;   // C(5)
    C_mat(2,1) = -REELMU;          // C(6)

    // CSUBB (next 2 columns)
    C_mat(0,2) = REELMU;           // C(7)
    C_mat(0,3) = 0.0;              // C(8)
    C_mat(1,2) = 0.0;              // C(9)
    C_mat(1,3) = -FLAMDA * REELMU * XSUBC;  // C(10)
    C_mat(2,2) = -FLAMDA * REELMU * XSUBC;  // C(11)
    C_mat(2,3) = REELMU;           // C(12)

    // CSUBC (last 2 columns)
    C_mat(0,4) = 0.0;              // C(13)
    C_mat(0,5) = 0.0;              // C(14)
    C_mat(1,4) = 0.0;              // C(15)
    C_mat(1,5) = FLAMDA;           // C(16)
    C_mat(2,4) = FLAMDA;           // C(17)
    C_mat(2,5) = 0.0;              // C(18)

    // Material properties: isotropic plane stress
    // G matrix (3x3) for plane stress: [G11, G12, G13; G12, G22, G23; G13, G23, G33]
    // For isotropic material: 
    // G11 = E/(1-nu^2), G22 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G33 = G/2 = E/(2*(1+nu))
    double G11 = E / (1.0 - nu*nu);
    double G12 = nu * E / (1.0 - nu*nu);
    double G22 = G11;
    double G33 = E / (2.0 * (1.0 + nu));

    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> G_mat;
    G_mat << G11, G12, 0.0,
             G12, G22, 0.0,
             0.0, 0.0, G33;

    // Compute stiffness matrix: K = VOL * (C * E^T)^T * G * (C * E^T)
    // C is 3x6, E^T is 3x3 => C * E^T is 3x6
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C_E_T = C_mat * E_mat.transpose();

    // Then K = VOL * (C * E^T)^T * G * (C * E^T) → 6x6
    Eigen::Matrix<double, 6, 3, Eigen::RowMajor> C_E_T_T = C_E_T.transpose();
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> K_full = VOL * C_E_T_T * G_mat * C_E_T;

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << double_to_scientific(K_full(i,j));
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}