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
            // Trim zeros after decimal before exponent
            size_t last_nonzero = epos - 1;
            while (last_nonzero > dot && s[last_nonzero] == '0') {
                --last_nonzero;
            }
            if (last_nonzero > dot && s[last_nonzero] == '.') {
                ++last_nonzero; // keep the decimal point if it's the last
            }
            s.erase(last_nonzero + 1, epos - last_nonzero - 1);
            // Now fix exponent sign: replace "e+0" with "e+"
            std::string exp_part = s.substr(epos);
            if (exp_part.length() >= 3 && exp_part[1] == '+') {
                if (exp_part[2] == '0') {
                    // remove the zero after +
                    s.replace(epos + 1, 2, "+");
                }
            }
        }
    }
    return s;
}

int main() {
    // Test case parameters
    const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    const double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    const double x3 = 1.0, y3 = 1.5, z3 = 0.0;
    const double E = 2.1e11;      // Pa
    const double nu = 0.3;
    const double t = 0.01;        // m

    // Material constants (isotropic)
    const double G11 = E / (1.0 - nu * nu);
    const double G12 = nu * E / (1.0 - nu * nu);
    const double G22 = G11;
    // G13=G23=G33=0 for plane stress membrane (not needed in 2D)

    // Constants from Fortran COMMON blocks (simplified for this test)
    const double DEGRA = M_PI / 180.0;  // degrees to radians
    const double ANGLE = 0.0;            // material angle = 0
    const double THETA = ANGLE * DEGRA;
    const double SINTH = std::sin(THETA);
    const double COSTH = std::cos(THETA);

    // Node IDs (for indexing, not used in computation but needed for structure)
    const int nodeA = 1, nodeB = 2, nodeC = 3;
    const int NPVT = nodeA; // pivot is node A (first node)

    // Compute geometry
    double E_vec[9]; // E matrix: 3x3, but only first two columns used (I, J, K vectors)
    
    // I-vector = B - A
    E_vec[0] = x2 - x1;
    E_vec[2] = y2 - y1;
    E_vec[4] = z2 - z1;
    
    double XSUBB = std::sqrt(E_vec[0]*E_vec[0] + E_vec[2]*E_vec[2] + E_vec[4]*E_vec[4]);
    if (XSUBB < 1.0e-6) {
        std::cerr << "Error: Zero length edge AB\n";
        return 1;
    }
    
    // Normalize I-vector
    E_vec[0] /= XSUBB;
    E_vec[2] /= XSUBB;
    E_vec[4] /= XSUBB;
    
    // Temporarily store C - A in E(2),E(4),E(6) (J-vector location)
    E_vec[1] = x3 - x1;
    E_vec[3] = y3 - y1;
    E_vec[5] = z3 - z1;
    
    // XSUBC = I . (C-A)
    double XSUBC = E_vec[0]*E_vec[1] + E_vec[2]*E_vec[3] + E_vec[4]*E_vec[5];
    
    // K-vector = I × (C-A) (non-normalized)
    E_vec[6] = E_vec[2]*E_vec[5] - E_vec[4]*E_vec[3];
    E_vec[7] = E_vec[4]*E_vec[1] - E_vec[0]*E_vec[5];
    E_vec[8] = E_vec[0]*E_vec[3] - E_vec[2]*E_vec[1];
    
    double YSUBC = std::sqrt(E_vec[6]*E_vec[6] + E_vec[7]*E_vec[7] + E_vec[8]*E_vec[8]);
    if (YSUBC < 1.0e-6) {
        std::cerr << "Error: Zero area triangle\n";
        return 1;
    }
    
    // Normalize K-vector
    E_vec[6] /= YSUBC;
    E_vec[7] /= YSUBC;
    E_vec[8] /= YSUBC;
    
    // J-vector = K × I
    E_vec[1] = E_vec[4]*E_vec[7] - E_vec[2]*E_vec[8];
    E_vec[3] = E_vec[0]*E_vec[8] - E_vec[4]*E_vec[6];
    E_vec[5] = E_vec[2]*E_vec[6] - E_vec[0]*E_vec[7];
    
    // Normalize J-vector
    double TEMP = std::sqrt(E_vec[1]*E_vec[1] + E_vec[3]*E_vec[3] + E_vec[5]*E_vec[5]);
    if (TEMP == 0.0) {
        std::cerr << "Error: Zero J-vector\n";
        return 1;
    }
    E_vec[1] /= TEMP;
    E_vec[3] /= TEMP;
    E_vec[5] /= TEMP;
    
    // Volume (area * thickness)
    double VOL = XSUBB * YSUBC * t / 2.0;
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB - 1.0;
    
    // Build C matrix (3x6): partitioned as [C_A | C_B | C_C], each 3x2
    // C is stored column-wise in Fortran, but we'll build as Eigen matrix
    Eigen::Matrix<double, 3, 6> C;
    C.setZero();
    
    // C_A (columns 0,1)
    C(0,0) = -REELMU;
    C(1,0) = 0.0;
    C(2,0) = 0.0;
    C(0,1) = FLAMDA * DELTA;
    C(1,1) = FLAMDA * DELTA; // Note: Fortran has C(4)=C(5)=FLAMDA*DELTA, so row1,col1 and row2,col1
    C(2,1) = -REELMU;
    
    // C_B (columns 2,3)
    C(0,2) = REELMU;
    C(1,2) = 0.0;
    C(2,2) = 0.0;
    C(0,3) = -FLAMDA * REELMU * XSUBC;
    C(1,3) = -FLAMDA * REELMU * XSUBC; // same as C(10) in Fortran
    C(2,3) = REELMU;
    
    // C_C (columns 4,5)
    C(0,4) = 0.0;
    C(1,4) = 0.0;
    C(2,4) = 0.0;
    C(0,5) = FLAMDA;
    C(1,5) = FLAMDA;
    C(2,5) = 0.0;
    
    // Build G matrix (3x3) for isotropic material in local coordinates
    // From Fortran: G(1)=G11, G(2)=G12, G(3)=G13, G(4)=G12, G(5)=G22, G(6)=G23, G(7)=G13, G(8)=G23, G(9)=G33
    // For plane stress membrane: G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G22 = G11, others zero
    Eigen::Matrix<double, 3, 3> G;
    G << G11, G12, 0.0,
         G12, G22, 0.0,
         0.0, 0.0, 0.0;
    
    // Build E matrix (3x2) from I and J vectors (first two columns of E_vec)
    // E = [I J] where I = [E_vec[0], E_vec[2], E_vec[4]]^T, J = [E_vec[1], E_vec[3], E_vec[5]]^T
    Eigen::Matrix<double, 3, 2> E_mat;
    E_mat << E_vec[0], E_vec[1],
              E_vec[2], E_vec[3],
              E_vec[4], E_vec[5];
    
    // Compute stiffness: K = VOL * (E * C)^T * G * (E * C)
    // First compute E*C (3x6)
    Eigen::Matrix<double, 3, 6> EC = E_mat * C;
    
    // Then compute (EC)^T * G * EC (6x6)
    Eigen::Matrix<double, 6, 6> K = VOL * EC.transpose() * G * EC;
    
    // Since the problem states "2 DOF per node: ux, uy", and nodes are A,B,C,
    // the ordering is: [u1, v1, u2, v2, u3, v3]
    // This matches the Fortran's KIJ(1..36) storage (column-major in Fortran, but we have row-major in Eigen)
    // Our K is already in the correct order: row i corresponds to DOF i, column j to DOF j
    
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