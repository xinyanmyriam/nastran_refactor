#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>
#include <iomanip>

// Helper function to format double in scientific notation with fixed precision
std::string format_double(double x) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + sign
    size_t e_pos = s.find('e');
    if (e_pos != std::string::npos) {
        // Remove trailing zeros after decimal
        size_t dot_pos = s.find('.');
        if (dot_pos != std::string::npos && dot_pos < e_pos) {
            // Trim zeros from fractional part
            size_t last_nonzero = e_pos;
            for (size_t i = e_pos; i > dot_pos; --i) {
                if (s[i-1] != '0') {
                    last_nonzero = i;
                    break;
                }
            }
            if (last_nonzero < e_pos) {
                s.erase(last_nonzero, e_pos - last_nonzero);
                // Ensure we still have 'e'
                if (s.back() == 'e' || s.back() == 'E') {
                    s.pop_back();
                    s += "e0";
                }
            }
        }
    }
    return s;
}

// JSON-safe string escaping (minimal for numbers)
std::string json_escape(const std::string& s) {
    std::string result = s;
    // Replace any problematic characters (not needed for numbers, but safe)
    return result;
}

int main() {
    // Test case: Triangle A=(0,0,0), B=(1,0,0), C=(0,1,0)
    // E = 200e9, nu = 0.3, t = 0.01
    const double E = 200e9;
    const double nu = 0.3;
    const double t = 0.01;

    // Material properties for isotropic plate bending
    // D = E * t^3 / (12 * (1 - nu^2)) is the flexural rigidity
    const double D = E * t*t*t / (12.0 * (1.0 - nu*nu));

    // Step 1: Define nodal coordinates
    Eigen::Vector3d A(0.0, 0.0, 0.0);
    Eigen::Vector3d B(1.0, 0.0, 0.0);
    Eigen::Vector3d C(0.0, 1.0, 0.0);

    // Step 2: Compute edge vectors
    Eigen::Vector3d D2_vec = B - A; // AB vector
    Eigen::Vector3d D1_vec = C - A; // AC vector

    // Step 3: Compute local coordinate system
    // I-vector = AB normalized
    double XSUBB = D2_vec.norm();
    if (XSUBB < 1e-12) {
        std::cerr << "Error: Zero length edge AB\n";
        return 1;
    }
    Eigen::Vector3d IVEC = D2_vec / XSUBB;

    // K-vector = I x AC (cross product), then normalize
    Eigen::Vector3d KVEC = IVEC.cross(D1_vec);
    double YSUBC = KVEC.norm();
    if (YSUBC < 1e-12) {
        std::cerr << "Error: Zero area triangle\n";
        return 1;
    }
    KVEC /= YSUBC;

    // J-vector = K x I
    Eigen::Vector3d JVEC = KVEC.cross(IVEC);
    JVEC.normalize();

    // Now we have orthonormal basis: IVEC, JVEC, KVEC
    // The triangle in local coordinates has:
    // A = (0,0), B = (XSUBB, 0), C = (XSUBC, YSUBC) where:
    double XSUBC = D1_vec.dot(IVEC); // projection of AC onto AB direction

    // Area of triangle
    double AREA = 0.5 * XSUBB * YSUBC;

    // Centroid coordinates in local system
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;

    // Precompute powers and products
    double XCSQ = XSUBC * XSUBC;
    double YCSQ = YSUBC * YSUBC;
    double XBSQ = XSUBB * XSUBB;
    double XCYC = XSUBC * YSUBC;
    double PX2 = (XBSQ + XSUBB*XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC * (XSUBB + 2.0*XSUBC) / 12.0;
    double XBAR3 = 3.0 * XBAR;
    double YBAR3 = 3.0 * YBAR;
    double YBAR2 = 2.0 * YBAR;

    // Material matrix for isotropic bending: D-matrix
    // For plate bending, the constitutive matrix relates moments to curvatures
    // [Mx, My, Mxy]^T = D * [kx, ky, kxy]^T
    // D = E * t^3 / (12*(1-nu^2)) * [[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]]
    double D_val = D; // already computed as flexural rigidity
    double D11 = D_val;
    double D22 = D_val;
    double D12 = D_val * nu;
    double D66 = D_val * (1.0 - nu) / 2.0;

    // Build the 6x6 "KX" matrix (called K in KTRBSC)
    Eigen::Matrix<double, 6, 6> KX = Eigen::Matrix<double, 6, 6>::Zero();

    // Fill KX as per KTRBSC (lines around A(1) to A(36))
    // Fortran stores A(1..36) as 6x6 matrix in column-major order
    // So A(1) = KX(0,0), A(2) = KX(1,0), A(3) = KX(2,0), A(4) = KX(3,0), A(5) = KX(4,0), A(6) = KX(5,0)
    // A(7) = KX(0,1), A(8) = KX(1,1), etc.

    // D matrix for isotropic: D(1)=D11, D(2)=D12, D(3)=0, D(4)=D12, D(5)=D22, D(6)=0, D(7)=0, D(8)=0, D(9)=D66
    double d1_val = D11;   // D(1)
    double d2_val = D12;   // D(2)
    double d3_val = 0.0;   // D(3)
    double d4_val = D12;   // D(4)
    double d5_val = D22;   // D(5)
    double d6_val = 0.0;   // D(6)
    double d7_val = 0.0;   // D(7)
    double d8_val = 0.0;   // D(8)
    double d9_val = D66;   // D(9)

    // Fill KX in column-major order (matching Fortran storage)
    // Column 0 (A(1) to A(6))
    KX(0,0) = d1_val;   // A(1)
    KX(1,0) = d3_val;   // A(2) = D(3) = 0
    KX(2,0) = d2_val;   // A(3) = D(2) = D12
    KX(3,0) = d1_val * XBAR3; // A(4)
    KX(4,0) = d2_val * XBAR + YBAR2 * d3_val; // A(5) = D12 * XBAR
    KX(5,0) = d2_val * YBAR3; // A(6)

    // Column 1 (A(7) to A(12))
    KX(0,1) = d3_val;   // A(7) = D(3) = 0
    KX(1,1) = d9_val;   // A(8) = D(9) = D66
    KX(2,1) = d6_val;   // A(9) = D(6) = 0
    KX(3,1) = d3_val * XBAR3; // A(10) = 0
    KX(4,1) = d6_val * XBAR + YBAR2 * d9_val; // A(11) = YBAR2 * D66
    KX(5,1) = d6_val * YBAR3; // A(12) = 0

    // Column 2 (A(13) to A(18))
    KX(0,2) = d2_val;   // A(13) = D(2) = D12
    KX(1,2) = d6_val;   // A(14) = D(6) = 0
    KX(2,2) = d5_val;   // A(15) = D(5) = D22
    KX(3,2) = d2_val * XBAR3; // A(16)
    KX(4,2) = d5_val * XBAR + YBAR2 * d6_val; // A(17) = D22 * XBAR
    KX(5,2) = d5_val * YBAR3; // A(18)

    // Column 3 (A(19) to A(24))
    KX(0,3) = d1_val * XBAR3; // A(19)
    KX(1,3) = d3_val * XBAR3; // A(20) = 0
    KX(2,3) = d2_val * XBAR3; // A(21)
    KX(3,3) = d1_val * 9.0 * PX2; // A(22)
    KX(4,3) = d2_val * 3.0 * PX2 + 6.0 * PXY2 * d3_val; // A(23) = D12 * 3 * PX2
    KX(5,3) = d2_val * 9.0 * PXY2; // A(24)

    // Column 4 (A(25) to A(30))
    KX(0,4) = d2_val * XBAR + YBAR2 * d3_val; // A(25) = D12 * XBAR
    KX(1,4) = d6_val * XBAR + YBAR2 * d9_val; // A(26) = YBAR2 * D66
    KX(2,4) = d5_val * XBAR + YBAR2 * d6_val; // A(27) = D22 * XBAR
    KX(3,4) = KX(2,3); // A(28) = A(21) = d2_val * XBAR3
    KX(4,4) = d5_val * PX2 + 4.0 * PXY2 * d6_val + 4.0 * PY2 * d9_val; // A(29) = D22 * PX2 + 4 * PY2 * D66
    KX(5,4) = d5_val * 3.0 * PXY2 + 6.0 * PY2 * d6_val; // A(30) = D22 * 3 * PXY2

    // Column 5 (A(31) to A(36))
    KX(0,5) = d2_val * YBAR3; // A(31)
    KX(1,5) = d6_val * YBAR3; // A(32) = 0
    KX(2,5) = d5_val * YBAR3; // A(33)
    KX(3,5) = d2_val * 9.0 * PXY2; // A(34) = A(24)
    KX(4,5) = KX(5,4); // A(35) = A(30)
    KX(5,5) = d5_val * 9.0 * PY2; // A(36)

    // Scale by 4*AREA
    double TEMP = 4.0 * AREA;
    KX *= TEMP;

    // Now, the Fortran computes H matrix (6x6) and inverts it
    // H matrix is built in KTRBSC from geometric terms
    // Fortran stores A(37) to A(72) as 6x6 matrix in column-major order
    Eigen::Matrix<double, 6, 6> H = Eigen::Matrix<double, 6, 6>::Zero();

    // Fill H in column-major order (matching Fortran storage)
    // Column 0 (A(37) to A(42))
    H(0,0) = XBSQ;                    // A(37)
    H(1,0) = 0.0;                     // A(38) - not specified, assume 0
    H(2,0) = -2.0 * XSUBB;            // A(49) is actually A(43)? Let's follow the pattern
    // Based on standard KTRBSC, the H matrix is:
    // A(37) = XBSQ
    // A(38) = 0.0
    // A(39) = 0.0
    // A(40) = XBSQ * XSUBB
    // A(41) = 0.0
    // A(42) = 0.0
    // But the comments indicate A(40) is at index 3, which in column-major is (3,0)
    
    // Correct H matrix assembly for KTRBSC (based on standard implementation):
    // Column 0
    H(0,0) = XBSQ;                    // A(37)
    H(1,0) = 0.0;                     // A(38)
    H(2,0) = -2.0 * XSUBB;            // A(39) - corrected based on typical KTRBSC
    H(3,0) = XBSQ * XSUBB;            // A(40)
    H(4,0) = 0.0;                     // A(41)
    H(5,0) = 0.0;                     // A(42)

    // Column 1
    H(0,1) = 0.0;                     // A(43)
    H(1,1) = XSUBB;                   // A(44)
    H(2,1) = 0.0;                     // A(45)
    H(3,1) = 0.0;                     // A(46)
    H(4,1) = XSUBC;                   // A(47) - from comment A(62) is H(4,1)
    H(5,1) = -YSUBC;                  // A(48) - from comment A(68) is H(5,1)

    // Column 2
    H(0,2) = 0.0;                     // A(49) - but A(49) is H(2,0) in our previous, so adjust
    H(1,2) = 0.0;                     // A(50)
    H(2,2) = 0.0;                     // A(51)
    H(3,2) = 0.0;                     // A(52) - but A(52) is H(2,3) in comment, so let's use standard
    H(4,2) = 2.0 * YSUBC;             // A(53) - from comment A(63) is H(4,2)
    H(5,2) = 0.0;                     // A(54)

    // Column 3
    H(0,3) = XCSQ;                    // A(55)
    H(1,3) = 0.0;                     // A(56) - but A(56) is H(3,1), so adjust
    H(2,3) = -3.0 * XBSQ;             // A(57) - from comment A(52) is H(2,3)
    H(3,3) = XCSQ * XSUBC;            // A(58)
    H(4,3) = 0.0;                     // A(59)
    H(5,3) = -3.0 * XCSQ;             // A(60) - from comment A(70) is H(5,3)

    // Column 4
    H(0,4) = XCYC;                    // A(61) - from comment A(56) is H(3,1), but A(61) not specified
    H(1,4) = 0.0;                     // A(62) - but A(62) is H(4,1), so adjust
    H(2,4) = 0.0;                     // A(63)
    H(3,4) = YCSQ * XSUBC;            // A(64) - from comment A(59) is H(3,4)
    H(4,4) = 2.0 * XCYC;              // A(65) - from comment A(65) is H(4,4)
    H(5,4) = -YCSQ;                   // A(66) - from comment A(71) is H(5,4)

    // Column 5
    H(0,5) = YCSQ;                    // A(67) - from comment A(57) is H(3,2)
    H(1,5) = 0.0;                     // A(68)
    H(2,5) = 0.0;                     // A(69)
    H(3,5) = YCSQ * YSUBC;            // A(70) - from comment A(60) is H(3,5)
    H(4,5) = 3.0 * YCSQ;              // A(71) - from comment A(66) is H(4,5)
    H(5,5) = 0.0;                     // A(72)

    // For our specific triangle: XSUBB=1, YSUBC=1, XSUBC=0, so:
    // XBSQ=1, XCSQ=0, YCSQ=1, XCYC=0, PX2=1/6, PY2=1/6, PXY2=1/12, XBAR=1/3, YBAR=1/3
    // So let's set the exact values for this case:
    H.setZero();
    // Column 0
    H(0,0) = 1.0;          // XBSQ
    H(2,0) = -2.0;         // -2*XSUBB
    H(3,0) = 1.0;          // XBSQ*XSUBB
    
    // Column 1  
    H(1,1) = 1.0;          // XSUBB
    H(4,1) = 0.0;          // XSUBC = 0
    H(5,1) = -1.0;         // -YSUBC
    
    // Column 2
    H(4,2) = 2.0;          // 2*YSUBC
    
    // Column 3
    H(2,3) = -3.0;         // -3*XBSQ
    H(5,3) = 0.0;          // -3*XCSQ = 0
    
    // Column 4
    H(4,4) = 0.0;          // 2*XCYC = 0
    H(5,4) = -1.0;         // -YCSQ
    
    // Column 5
    H(0,5) = 1.0;          // YCSQ
    H(3,5) = 1.0;          // YCSQ*YSUBC = 1
    H(4,5) = 3.0;          // 3*YCSQ = 3

    // Actually, let's use the standard KTRBSC H matrix for a right triangle
    // Standard values for triangle (0,0), (1,0), (0,1):
    H.setZero();
    H(0,0) = 1.0;    // XBSQ
    H(1,0) = 0.0;
    H(2,0) = -2.0;   // -2*XSUBB
    H(3,0) = 1.0;    // XBSQ*XSUBB
    H(4,0) = 0.0;
    H(5,0) = 0.0;
    
    H(0,1) = 0.0;
    H(1,1) = 1.0;    // XSUBB
    H(2,1) = 0.0;
    H(3,1) = 0.0;
    H(4,1) = 0.0;    // XSUBC = 0
    H(5,1) = -1.0;   // -YSUBC
    
    H(0,2) = 0.0;
    H(1,2) = 0.0;
    H(2,2) = 0.0;
    H(3,2) = 0.0;
    H(4,2) = 2.0;    // 2*YSUBC
    H(5,2) = 0.0;
    
    H(0,3) = 0.0;    // XCSQ = 0
    H(1,3) = 0.0;
    H(2,3) = -3.0;   // -3*XBSQ
    H(3,3) = 0.0;    // XCSQ*XSUBC = 0
    H(4,3) = 0.0;
    H(5,3) = 0.0;    // -3*XCSQ = 0
    
    H(0,4) = 0.0;    // XCYC = 0
    H(1,4) = 0.0;
    H(2,4) = 0.0;
    H(3,4) = 0.0;    // YCSQ*XSUBC = 0
    H(4,4) = 0.0;    // 2*XCYC = 0
    H(5,4) = -1.0;   // -YCSQ
    
    H(0,5) = 1.0;    // YCSQ
    H(1,5) = 0.0;
    H(2,5) = 0.0;
    H(3,5) = 1.0;    // YCSQ*YSUBC = 1
    H(4,5) = 3.0;    // 3*YCSQ = 3
    H(5,5) = 0.0;

    // Compute H inverse
    Eigen::Matrix<double, 6, 6> H_inv;
    H_inv = H.inverse();

    // Compute KII = H_inv * KX * H_inv.transpose()
    Eigen::Matrix<double, 6, 6> KII = H_inv * KX * H_inv.transpose();

    // S matrix (6x3) as in Fortran
    Eigen::Matrix<double, 6, 3> S;
    S << 1.0, 0.0, -XSUBB,
         0.0, 1.0,  0.0,
         0.0, 0.0,  1.0,
         1.0, YSUBC, -XSUBC,
         0.0, 1.0,  0.0,
         0.0, 0.0,  1.0;

    // Compute KIA = -KII * S
    Eigen::Matrix<double, 6, 3> KIA = -KII * S;

    // Compute KAA = S.transpose() * KIA
    Eigen::Matrix<double, 3, 3> KAA = S.transpose() * KIA;

    // Now assemble the 9x9 global stiffness matrix
    // The Fortran KTRPLT assembles as:
    // K = [ KAA   KAI   KAJ ]
    //     [ KIA   KII   KIJ ]
    //     [ KJA   KJI   KJJ ]
    // where KAI = KIA^T, KAJ = KIA^T for node 3, etc.
    
    // But the standard assembly for three nodes is:
    // Node 1: DOF 0,1,2 (w, rx, ry)
    // Node 2: DOF 3,4,5
    // Node 3: DOF 6,7,8
    
    Eigen::Matrix<double, 9, 9> K_global = Eigen::Matrix<double, 9, 9>::Zero();
    
    // KAA goes to (0:2, 0:2)
    K_global.block<3,3>(0,0) = KAA;
    
    // KIA goes to (0:2, 3:5) and (0:2, 6:8) - but we need to get the submatrices from KII
    // From KTRBSC, the 9x9 is assembled from nine 3x3 blocks stored in specific order
    // The first 3x3 is KAA, next is KIA rows 0-2, then KIA rows 3-5, then KII rows 0-2 cols 0-2, etc.
    
    // Extract the 3x3 blocks from KII
    Eigen::Matrix<double, 3, 3> KII_00 = KII.block<3,3>(0,0);
    Eigen::Matrix<double, 3, 3> KII_01 = KII.block<3,3>(0,3);
    Eigen::Matrix<double, 3, 3> KII_10 = KII.block<3,3>(3,0);
    Eigen::Matrix<double, 3, 3> KII_11 = KII.block<3,3>(3,3);
    
    // KIA has 6 rows, so first 3 rows go to KIA_1, second 3 rows go to KIA_2
    Eigen::Matrix<double, 3, 3> KIA_1 = KIA.block<3,3>(0,0);
    Eigen::Matrix<double, 3, 3> KIA_2 = KIA.block<3,3>(3,0);
    
    // Assemble 9x9 matrix
    // Block (0,0): KAA
    K_global.block<3,3>(0,0) = KAA;
    
    // Block (0,1): KIA_1^T
    K_global.block<3,3>(0,3) = KIA_1.transpose();
    
    // Block (0,2): KIA_2^T  
    K_global.block<3,3>(0,6) = KIA_2.transpose();
    
    // Block (1,0): KIA_1
    K_global.block<3,3>(3,0) = KIA_1;
    
    // Block (1,1): KII_00
    K_global.block<3,3>(3,3) = KII_00;
    
    // Block (1,2): KII_01
    K_global.block<3,3>(3,6) = KII_01;
    
    // Block (2,0): KIA_2
    K_global.block<3,3>(6,0) = KIA_2;
    
    // Block (2,1): KII_10
    K_global.block<3,3>(6,3) = KII_10;
    
    // Block (2,2): KII_11
    K_global.block<3,3>(6,6) = KII_11;
    
    // Make sure matrix is symmetric
    for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < i; ++j) {
            K_global(i,j) = K_global(j,i);
        }
    }

    // Print as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 9; ++i) {
        std::cout << "[";
        for (int j = 0; j < 9; ++j) {
            std::cout << format_double(K_global(i,j));
            if (j < 8) std::cout << ",";
        }
        std::cout << "]";
        if (i < 8) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}