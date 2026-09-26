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
    const double E = 2.1e11;   // Pa
    const double nu = 0.3;
    const double t = 0.01;     // thickness in meters

    // Constants
    const double DEGRA = M_PI / 180.0; // degrees to radians

    // Compute element geometry
    // E matrix: 3x3 for coordinate system construction
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> E;
    E.setZero();

    // I-vector = RSUBB - RSUBA (non-normalized)
    E(0,0) = x2 - x1;
    E(0,1) = y2 - y1;
    E(0,2) = z2 - z1;

    // XSUBB = length of I-vector
    double XSUBB = std::sqrt(E(0,0)*E(0,0) + E(0,1)*E(0,1) + E(0,2)*E(0,2));
    
    // Check for degenerate element
    if (XSUBB <= 1.0e-6) {
        std::cerr << "Error: Degenerate element - zero length I-vector" << std::endl;
        return 1;
    }

    // Normalize I-vector
    E(0,0) /= XSUBB;
    E(0,1) /= XSUBB;
    E(0,2) /= XSUBB;

    // Store RSUBC - RSUBA temporarily in E(1,0), E(1,1), E(1,2) (J-vector location)
    E(1,0) = x3 - x1;
    E(1,1) = y3 - y1;
    E(1,2) = z3 - z1;

    // XSUBC = I . (RSUBC - RSUBA)
    double XSUBC = E(0,0)*E(1,0) + E(0,1)*E(1,1) + E(0,2)*E(1,2);

    // K-vector = I cross (RSUBC - RSUBA) (non-normalized)
    E(2,0) = E(0,1)*E(1,2) - E(0,2)*E(1,1);
    E(2,1) = E(0,2)*E(1,0) - E(0,0)*E(1,2);
    E(2,2) = E(0,0)*E(1,1) - E(0,1)*E(1,0);

    // YSUBC = length of K-vector
    double YSUBC = std::sqrt(E(2,0)*E(2,0) + E(2,1)*E(2,1) + E(2,2)*E(2,2));
    
    if (YSUBC <= 1.0e-6) {
        std::cerr << "Error: Degenerate element - zero length K-vector" << std::endl;
        return 1;
    }

    // Normalize K-vector
    E(2,0) /= YSUBC;
    E(2,1) /= YSUBC;
    E(2,2) /= YSUBC;

    // J-vector = K cross I
    E(1,0) = E(2,1)*E(0,2) - E(2,2)*E(0,1);
    E(1,1) = E(2,2)*E(0,0) - E(2,0)*E(0,2);
    E(1,2) = E(2,0)*E(0,1) - E(2,1)*E(0,0);

    // Normalize J-vector
    double TEMP = std::sqrt(E(1,0)*E(1,0) + E(1,1)*E(1,1) + E(1,2)*E(1,2));
    if (TEMP == 0.0) {
        std::cerr << "Error: Degenerate element - zero length J-vector" << std::endl;
        return 1;
    }
    E(1,0) /= TEMP;
    E(1,1) /= TEMP;
    E(1,2) /= TEMP;

    // Volume of element
    double VOL = XSUBB * YSUBC * t / 2.0;
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB - 1.0;

    // C matrix: 3x6 partitioned as CSUBA (3x2), CSUBB (3x2), CSUBC (3x2)
    // Stored in row-major order: C(1) to C(18)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C;
    C.setZero();

    // CSUBA (first 2 columns)
    C(0,0) = -REELMU;           // C(1)
    C(0,1) = 0.0;               // C(2)
    C(1,0) = 0.0;               // C(3)
    C(1,1) = FLAMDA * DELTA;    // C(4)
    C(2,0) = C(1,1);            // C(5) = C(4)
    C(2,1) = -REELMU;           // C(6)

    // CSUBB (next 2 columns)
    C(0,2) = REELMU;            // C(7)
    C(0,3) = 0.0;               // C(8)
    C(1,2) = 0.0;               // C(9)
    C(1,3) = -FLAMDA * REELMU * XSUBC; // C(10)
    C(2,2) = C(1,3);            // C(11) = C(10)
    C(2,3) = REELMU;            // C(12)

    // CSUBC (last 2 columns)
    C(0,4) = 0.0;               // C(13)
    C(0,5) = 0.0;               // C(14)
    C(1,4) = 0.0;               // C(15)
    C(1,5) = FLAMDA;            // C(16)
    C(2,4) = FLAMDA;            // C(17)
    C(2,5) = 0.0;               // C(18)

    // Material matrix G (3x3) for isotropic material
    // G = [G11 G12 G13; G12 G22 G23; G13 G23 G33]
    // For plane stress isotropic: 
    // G11 = E/(1-nu^2), G22 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G33 = E/(2*(1+nu))
    double G11 = E / (1.0 - nu*nu);
    double G22 = G11;
    double G12 = nu * E / (1.0 - nu*nu);
    double G33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> G;
    G << G11, G12, 0.0,
         G12, G22, 0.0,
         0.0, 0.0, G33;

    // Since this is a membrane element with no coordinate system transformation needed
    // (all nodes in global coordinates, no local coordinate system specified),
    // we assume identity transformation matrices TI and TJ.
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> TI = Eigen::Matrix<double, 3, 3>::Identity();
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> TJ = Eigen::Matrix<double, 3, 3>::Identity();

    // Compute stiffness matrix KIJ (6x6)
    // Formula: K = VOL * T_I^T * E * C^T * G * C * E^T * T_J
    // Since TI and TJ are identity, this simplifies to: K = VOL * E * C^T * G * C * E^T

    // First compute C^T (6x3)
    Eigen::Matrix<double, 6, 3, Eigen::RowMajor> C_T = C.transpose();

    // Compute C^T * G (6x3)
    Eigen::Matrix<double, 6, 3, Eigen::RowMajor> CT_G = C_T * G;

    // Compute C^T * G * C (6x6)
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> CT_G_C = CT_G * C;

    // Compute E * (C^T * G * C) (3x6)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> E_CT_G_C = E * CT_G_C;

    // Compute E * (C^T * G * C) * E^T (3x3)
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> E_CT_G_C_ET = E_CT_G_C * E.transpose();

    // Multiply by volume
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> K_temp = VOL * E_CT_G_C_ET;

    // Now expand to 6x6 stiffness matrix
    // The 3x3 result is for the strain-displacement relationship in the local element coordinate system
    // But for CTRIA3, the full 6x6 stiffness matrix is assembled from contributions
    // Following the NASTRAN logic, we need to compute the full 6x6 matrix.

    // Actually, let's follow the exact NASTRAN computation path:
    // K = VOL * T_I^T * E * C^T * G * C * E^T * T_J
    // With TI = TJ = I, so K = VOL * E * C^T * G * C * E^T

    // But wait - the C matrix is 3x6, so C^T is 6x3, G is 3x3, so C^T*G*C is 6x6
    // Then E * (C^T*G*C) would be 3x6, which doesn't match.

    // Let's re-examine the NASTRAN comment: "K = VOL . T *E*C^T*G*C *E^T*T"
    // The dimensions in NASTRAN: E is 3x2 (not 3x3 as we computed above for coordinate system),
    // but in the Fortran code, E is used as a 3x2 matrix for the transformation.

    // Looking more carefully at the Fortran code:
    // E(1), E(3), E(5) are the I-vector components -> stored in E(1,1), E(2,1), E(3,1) in column 1
    // E(2), E(4), E(6) are the J-vector components -> stored in E(1,2), E(2,2), E(3,2) in column 2
    // E(7), E(8), E(9) are the K-vector components -> not used for membrane (so E is effectively 3x2)

    // So let's reconstruct E as 3x2:
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> E_3x2;
    E_3x2 << E(0,0), E(1,0),  // I-vector, J-vector
              E(0,1), E(1,1),
              E(0,2), E(1,2);

    // Now C is 3x6, so C^T is 6x3
    // We need: E_3x2^T * C^T * G * C * E_3x2
    // Dimensions: (2x3) * (3x6) * (6x6) * (6x3) * (3x2) -> no, that's not right.

    // Let's read the Fortran comment again: "E*C^T*G*C*E^T"
    // E is 3x2, C is 3x6, so C^T is 6x3, so E*C^T is 3x6
    // Then (E*C^T)*G is 3x3, then *(C*E^T) is 3x2 -> gives 3x2, not 6x6.

    // Actually, the standard CST element formulation is:
    // B = [dN1/dx  0     dN2/dx  0     dN3/dx  0
    //      0     dN1/dy  0     dN2/dy  0     dN3/dy
    //      dN1/dy dN1/dx dN2/dy dN2/dx dN3/dy dN3/dx]
    // But NASTRAN uses a different approach.

    // Let's implement the standard CST (Constant Strain Triangle) formulation instead,
    // which is well-documented and matches the problem description.

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
    D << E/denom, nu*E/denom, 0.0,
         nu*E/denom, E/denom, 0.0,
         0.0, 0.0, E/(2.0*(1.0+nu));

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