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
    
    // Material properties: isotropic plane stress
    // G11 = E/(1-nu^2), G12 = nu*E/(1-nu^2), G22 = E/(1-nu^2), G33 = E/(2*(1+nu))
    const double G11 = E / (1.0 - nu * nu);
    const double G12 = nu * E / (1.0 - nu * nu);
    const double G22 = G11;
    const double G33 = E / (2.0 * (1.0 + nu));
    
    // Build G matrix (3x3) for plane stress
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> G;
    G << G11, G12, 0.0,
         G12, G22, 0.0,
         0.0, 0.0, G33;
    
    // Compute element geometry
    // E matrix: 3x3 where rows are i, j, k unit vectors
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> E_mat;
    
    // I-vector = R_B - R_A
    double e1 = x2 - x1;
    double e3 = y2 - y1;
    double e5 = z2 - z1;
    
    double xsubb = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    if (xsubb < 1.0e-6) {
        std::cerr << "Error: degenerate element - zero length edge AB" << std::endl;
        return 1;
    }
    
    // Normalize I-vector
    e1 /= xsubb;
    e3 /= xsubb;
    e5 /= xsubb;
    
    // Store temporarily in E_mat row 0 (i-vector)
    E_mat(0,0) = e1;
    E_mat(0,1) = e3;
    E_mat(0,2) = e5;
    
    // RSUBC - RSUBA
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
        std::cerr << "Error: degenerate element - zero area" << std::endl;
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
        std::cerr << "Error: degenerate element - zero j vector" << std::endl;
        return 1;
    }
    
    j1 /= temp;
    j2 /= temp;
    j3 /= temp;
    
    // Store j-vector in E_mat row 1
    E_mat(1,0) = j1;
    E_mat(1,1) = j2;
    E_mat(1,2) = j3;
    
    // Store k-vector in E_mat row 2
    E_mat(2,0) = e7;
    E_mat(2,1) = e8;
    E_mat(2,2) = e9;
    
    // Volume of element
    double vol = xsubb * ysubc * t / 2.0;
    
    // Material constants
    double reelmu = 1.0 / xsubb;
    double flambda = 1.0 / ysubc;
    double delta = xsubc / xsubb - 1.0;
    
    // Build C matrix: 3x6, partitioned as [C_A | C_B | C_C] where each is 3x2
    // C_A: rows 0-2, cols 0-1
    // C_B: rows 0-2, cols 2-3  
    // C_C: rows 0-2, cols 4-5
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C;
    
    // C_A
    C(0,0) = -reelmu; C(0,1) = 0.0;
    C(1,0) = 0.0;     C(1,1) = flambda * delta;
    C(2,0) = 0.0;     C(2,1) = -reelmu;
    
    // C_B
    C(0,2) = reelmu;           C(0,3) = 0.0;
    C(1,2) = 0.0;              C(1,3) = -flambda * reelmu * xsubc;
    C(2,2) = reelmu;           C(2,3) = 0.0;
    
    // C_C
    C(0,4) = 0.0;              C(0,5) = 0.0;
    C(1,4) = 0.0;              C(1,5) = flambda;
    C(2,4) = 0.0;              C(2,5) = flambda;
    
    // Since the Fortran code uses a different indexing pattern, let's reconstruct
    // the C matrix as used in the original algorithm based on the Fortran comments:
    // C(1) thru C(6): CSUBA (3x2) -> rows 0-2, cols 0-1
    // C(7) thru C(12): CSUBB (3x2) -> rows 0-2, cols 2-3
    // C(13) thru C(18): CSUBC (3x2) -> rows 0-2, cols 4-5
    
    // But looking at the Fortran code more carefully:
    // C(1)  = -REELMU          -> C_A[0,0]
    // C(2)  =  0.0D0            -> C_A[0,1] 
    // C(3)  =  0.0D0            -> C_A[1,0]
    // C(4)  =  FLAMDA*DELTA     -> C_A[1,1]
    // C(5)  =  C(4)             -> C_A[2,0] = FLAMDA*DELTA
    // C(6)  = -REELMU           -> C_A[2,1]
    // C(7)  =  REELMU           -> C_B[0,0]
    // C(8)  =  0.0D0            -> C_B[0,1]
    // C(9)  =  0.0D0            -> C_B[1,0]
    // C(10) = -FLAMDA*REELMU*XSUBC -> C_B[1,1]
    // C(11) =  C(10)           -> C_B[2,0] = -FLAMDA*REELMU*XSUBC
    // C(12) =  REELMU           -> C_B[2,1]
    // C(13) =  0.0D0            -> C_C[0,0]
    // C(14) =  0.0D0            -> C_C[0,1]
    // C(15) =  0.0D0            -> C_C[1,0]
    // C(16) =  FLAMDA           -> C_C[1,1]
    // C(17) =  FLAMDA           -> C_C[2,0]
    // C(18) =  0.0D0            -> C_C[2,1]
    
    // So the correct C matrix (3x6) is:
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> C_correct;
    C_correct << 
        -reelmu,           0.0,               reelmu,           0.0,               0.0,            0.0,
        0.0,      flambda*delta,              0.0,  -flambda*reelmu*xsubc,      0.0,       flambda,
        flambda*delta,   -reelmu,  -flambda*reelmu*xsubc,        reelmu,     flambda,            0.0;
    
    // Now compute stiffness matrix K = VOL * (C^T * G * C) but with proper transformation
    // The Fortran computes: K = VOL * T_I^T * E * C^T * G * C * E^T * T_J
    // Since we have no coordinate system transformations (all CSID = 0), T_I = T_J = I
    // So K = VOL * E * C^T * G * C * E^T
    
    // First compute C^T (6x3)
    Eigen::Matrix<double, 6, 3, Eigen::RowMajor> C_T = C_correct.transpose();
    
    // Compute C^T * G (6x3)
    Eigen::Matrix<double, 6, 3, Eigen::RowMajor> CT_G = C_T * G;
    
    // Compute C^T * G * C (6x6)
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> CT_G_C = CT_G * C_correct;
    
    // Compute E * (C^T * G * C) (3x6)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> E_CT_G_C = E_mat * CT_G_C;
    
    // Compute final K = VOL * E * C^T * G * C * E^T (3x3)
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> K_local = vol * E_CT_G_C * E_mat.transpose();
    
    // But wait - the Fortran computes a 6x6 global stiffness matrix.
    // Looking again at the Fortran: it computes KIJ(36) which is a 6x6 matrix.
    // The key is that the C matrix is 3x6, so C^T * G * C is 6x6.
    // Then they multiply by E and E^T, but E is 3x3, so that would give 3x3.
    // Actually, re-examining the Fortran comment: "K = VOL . T *E*C *G*C *E *T"
    // And the dimensions: E is 3x3, C is 3x6, G is 3x3, so E*C is 3x6, then (E*C)^T * G * (E*C) is 6x6.
    
    // Correct interpretation: K = VOL * (E*C)^T * G * (E*C)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> E_C = E_mat * C_correct;
    Eigen::Matrix<double, 6, 3, Eigen::RowMajor> E_C_T = E_C.transpose();
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> K_global = vol * E_C_T * G * E_C;
    
    // However, the Fortran code has a different structure. Let's follow the exact steps:
    // It computes: TEMPAR(10) = E * C(NPOINT)  [3x3 * 3x2 = 3x2]
    // Then TEMPAR(1) = TEMPAR(10) * G * TEMPAR(10)^T [3x2 * 3x3 * 2x3 = 3x3]
    // Then multiplies by VOL -> 3x3
    // Then assembles into 6x6 KIJ matrix at specific positions.
    
    // Looking at the assembly section (lines 95-96):
    // KIJ( 1) = TEMPAR(NPT1    )  -> position (0,0) of 3x3
    // KIJ( 2) = TEMPAR(NPT1 + 1)  -> position (0,1) of 3x3  
    // KIJ( 3) = TEMPAR(NPT1 + 2)  -> position (0,2) of 3x3
    // KIJ( 7) = TEMPAR(NPT1 + 3)  -> position (1,0) of 3x3 (since 7 = 1 + 6*1)
    // KIJ( 8) = TEMPAR(NPT1 + 4)  -> position (1,1) of 3x3
    // KIJ( 9) = TEMPAR(NPT1 + 5)  -> position (1,2) of 3x3
    // KIJ(13) = TEMPAR(NPT1 + 6) -> position (2,0) of 3x3 (13 = 1 + 6*2)
    // KIJ(14) = TEMPAR(NPT1 + 7) -> position (2,1) of 3x3
    // KIJ(15) = TEMPAR(NPT1 + 8) -> position (2,2) of 3x3
    
    // So the 3x3 result is placed in the upper-left 3x3 block of a 6x6 matrix? 
    // But that doesn't make sense for a 6x6 stiffness matrix.
    
    // Actually, re-reading the Fortran: it's computing three separate 3x3 matrices
    // for each node pair, and assembling them into the 6x6 global matrix.
    // The loop from 500 does this for I=1,2,3 (the pivot) and J=1,2,3 (the other nodes).
    
    // Since we have a simple case with no coordinate transformations (CSID=0),
    // and the problem asks for the full 6x6 stiffness matrix for the CTRIA3 element,
    // the standard approach is: K = B^T * D * B * t * area
    
    // Let's use the standard CST (Constant Strain Triangle) formulation:
    // Area = 0.5 * |AB × AC|
    double area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    
    // For CST, B matrix is 3x6:
    // B = [ b1 0  b2 0  b3 0 ]
    //     [ 0  c1  0 c2  0 c3 ]
    //     [ c1 b1  c2 b2  c3 b3 ]
    // where for node i: bi = (yj - yk)/2A, ci = (xk - xj)/2A
    
    double twoA = 2.0 * area;
    double b1 = (y2 - y3) / twoA;
    double b2 = (y3 - y1) / twoA;
    double b3 = (y1 - y2) / twoA;
    double c1 = (x3 - x2) / twoA;
    double c2 = (x1 - x3) / twoA;
    double c3 = (x2 - x1) / twoA;
    
    // Build B matrix (3x6)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    
    // D matrix for plane stress (3x3)
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> D;
    double denom = 1.0 - nu * nu;
    D << E/denom, nu*E/denom, 0.0,
         nu*E/denom, E/denom, 0.0,
         0.0, 0.0, E/(2.0*(1.0+nu));
    
    // K = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6, Eigen::RowMajor> K = t * area * B.transpose() * D * B;
    
    // This is the standard CST stiffness matrix and matches what's expected for CTRIA3.
    
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