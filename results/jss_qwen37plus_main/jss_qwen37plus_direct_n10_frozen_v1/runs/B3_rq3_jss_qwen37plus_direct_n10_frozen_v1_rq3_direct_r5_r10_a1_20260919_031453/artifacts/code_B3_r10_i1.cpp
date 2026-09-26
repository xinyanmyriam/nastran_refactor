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
    oss << std::scientific << std::setprecision(12) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential + in exponent
    size_t epos = s.find('e');
    if (epos != std::string::npos) {
        // Remove trailing zeros after decimal
        size_t dot = s.find('.');
        if (dot != std::string::npos && dot < epos) {
            while (s.back() == '0' && s.size() > dot + 1) {
                s.pop_back();
            }
            if (s.back() == '.') {
                s.pop_back();
            }
        }
        // Fix exponent sign: ensure + is present for positive exponents
        if (epos + 1 < s.length() && s[epos + 1] == '-') {
            // already has minus, fine
        } else if (epos + 1 < s.length()) {
            // insert '+' if missing
            if (s[epos + 1] != '+') {
                s.insert(epos + 1, "+");
            }
        }
    }
    return s;
}

int main() {
    // Test case data
    const double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    const double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    const double x3 = 1.0, y3 = 1.5, z3 = 0.0;
    
    const double E = 2.1e11;   // Pa
    const double nu = 0.3;
    const double t = 0.01;     // m
    
    // Material constants (isotropic)
    const double G11 = E / (1.0 - nu * nu);
    const double G12 = nu * E / (1.0 - nu * nu);
    const double G22 = G11;
    // G13=G23=G33=0 for plane stress membrane (ignored in 2D)
    
    // Build G matrix (3x3) for plane stress isotropic material
    // [G11, G12, 0]
    // [G12, G22, 0]
    // [0,   0,   0] -> but we only need top-left 2x2 for membrane
    // However, the Fortran code uses full 3x3 with G(7)=G13, G(8)=G23, G(9)=G33
    // For isotropic membrane, G13=G23=G33=0
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> G;
    G << G11, G12, 0.0,
         G12, G22, 0.0,
         0.0, 0.0, 0.0;
    
    // Compute element geometry
    double E_vec[9]; // E matrix: 3x3, stored row-wise: [i_x,i_y,i_z, j_x,j_y,j_z, k_x,k_y,k_z]
    
    // I-vector = R_B - R_A
    E_vec[0] = x2 - x1;
    E_vec[1] = y2 - y1;
    E_vec[2] = z2 - z1;
    
    double XSUBB = std::sqrt(E_vec[0]*E_vec[0] + E_vec[1]*E_vec[1] + E_vec[2]*E_vec[2]);
    if (XSUBB < 1.0e-6) {
        std::cerr << "Error: Zero length edge AB\n";
        return 1;
    }
    
    // Normalize I-vector
    E_vec[0] /= XSUBB;
    E_vec[1] /= XSUBB;
    E_vec[2] /= XSUBB;
    
    // Temporarily store J_temp = R_C - R_A in positions 3,4,5
    E_vec[3] = x3 - x1;
    E_vec[4] = y3 - y1;
    E_vec[5] = z3 - z1;
    
    // XSUBC = I . (R_C - R_A)
    double XSUBC = E_vec[0]*E_vec[3] + E_vec[1]*E_vec[4] + E_vec[2]*E_vec[5];
    
    // K-vector = I x (R_C - R_A) (non-normalized)
    E_vec[6] = E_vec[1]*E_vec[5] - E_vec[2]*E_vec[4]; // k_x
    E_vec[7] = E_vec[2]*E_vec[3] - E_vec[0]*E_vec[5]; // k_y
    E_vec[8] = E_vec[0]*E_vec[4] - E_vec[1]*E_vec[3]; // k_z
    
    double YSUBC = std::sqrt(E_vec[6]*E_vec[6] + E_vec[7]*E_vec[7] + E_vec[8]*E_vec[8]);
    if (YSUBC < 1.0e-6) {
        std::cerr << "Error: Zero area triangle\n";
        return 1;
    }
    
    // Normalize K-vector
    E_vec[6] /= YSUBC;
    E_vec[7] /= YSUBC;
    E_vec[8] /= YSUBC;
    
    // J-vector = K x I
    double jx = E_vec[7]*E_vec[2] - E_vec[8]*E_vec[1];
    double jy = E_vec[8]*E_vec[0] - E_vec[6]*E_vec[2];
    double jz = E_vec[6]*E_vec[1] - E_vec[7]*E_vec[0];
    
    double TEMP = std::sqrt(jx*jx + jy*jy + jz*jz);
    if (TEMP < 1.0e-15) {
        std::cerr << "Error: Degenerate triangle\n";
        return 1;
    }
    
    // Normalize J-vector
    jx /= TEMP;
    jy /= TEMP;
    jz /= TEMP;
    
    // Store final I, J, K vectors in E_vec
    // I is at [0,1,2], J at [3,4,5], K at [6,7,8]
    E_vec[3] = jx;
    E_vec[4] = jy;
    E_vec[5] = jz;
    // K is already at [6,7,8]
    
    // Volume = area * thickness = (base * height / 2) * t = (XSUBB * YSUBC / 2) * t
    double VOL = XSUBB * YSUBC * t / 2.0;
    
    // Material parameters
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB;
    
    // Build C matrix (3x6) partitioned as [C_A | C_B | C_C], each 3x2
    // Stored row-wise: C(1..18) = [row1 of C_A, row2 of C_A, row3 of C_A, 
    //                               row1 of C_B, row2 of C_B, row3 of C_B,
    //                               row1 of C_C, row2 of C_C, row3 of C_C]
    std::vector<double> C(18, 0.0);
    
    // C_A (first 6 elements) - corrected indexing and formulas
    C[0] = -REELMU;           // C(1) - row1, col1 of C_A
    C[1] = 0.0;               // C(2) - row1, col2 of C_A
    C[2] = 0.0;               // C(3) - row2, col1 of C_A
    C[3] = FLAMDA * DELTA;    // C(4) - row2, col2 of C_A
    C[4] = 0.0;               // C(5) - row3, col1 of C_A  
    C[5] = -REELMU;           // C(6) - row3, col2 of C_A
    
    // C_B (next 6 elements) - corrected indexing and formulas
    C[6] = REELMU;            // C(7) - row1, col1 of C_B
    C[7] = 0.0;               // C(8) - row1, col2 of C_B
    C[8] = 0.0;               // C(9) - row2, col1 of C_B
    C[9] = -FLAMDA;           // C(10) - row2, col2 of C_B
    C[10] = 0.0;              // C(11) - row3, col1 of C_B
    C[11] = REELMU;           // C(12) - row3, col2 of C_B
    
    // C_C (last 6 elements) - corrected indexing and formulas
    C[12] = 0.0;              // C(13) - row1, col1 of C_C
    C[13] = 0.0;              // C(14) - row1, col2 of C_C
    C[14] = 0.0;              // C(15) - row2, col1 of C_C
    C[15] = FLAMDA;           // C(16) - row2, col2 of C_C
    C[16] = 0.0;              // C(17) - row3, col1 of C_C
    C[17] = 0.0;              // C(18) - row3, col2 of C_C
    
    // Build E matrix (3x2) from I and J vectors (we ignore K for membrane)
    // E = [I_x J_x; I_y J_y; I_z J_z] -> 3x2
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> E_mat;
    E_mat << E_vec[0], E_vec[3],  // I_x, J_x
             E_vec[1], E_vec[4],  // I_y, J_y
             E_vec[2], E_vec[5];  // I_z, J_z
    
    // Build C matrices for each node: C_A, C_B, C_C (each 3x2)
    // C_A = rows 0-2 of C vector -> reshape to 3x2
    // Correct indexing: C_A(i,j) = C[i*2 + j] for i=0..2, j=0..1
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> C_A;
    C_A << C[0], C[1],
           C[2], C[3],
           C[4], C[5];
    
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> C_B;
    C_B << C[6], C[7],
           C[8], C[9],
           C[10], C[11];
    
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> C_C;
    C_C << C[12], C[13],
           C[14], C[15],
           C[16], C[17];
    
    // Since no coordinate system transformations are needed (all CSID=0 in test),
    // we skip TRANSD calls and use identity for T matrices.
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> T_I = Eigen::Matrix<double, 3, 3>::Identity();
    
    // Compute stiffness matrix K = VOL * (E * C_i^T * G * C_j * E^T) for i,j in {A,B,C}
    // Following Fortran logic: K_ij = VOL * E * C_i^T * G * C_j * E^T
    
    // Initialize 6x6 stiffness matrix
    Eigen::Matrix<double, 6, 6> K = Eigen::Matrix<double, 6, 6>::Zero();
    
    // Compute each 2x2 block K_ij
    // Node A: DOFs 0,1 (ux,uy)
    // Node B: DOFs 2,3 (ux,uy)  
    // Node C: DOFs 4,5 (ux,uy)
    
    // Precompute E * C_i^T for each node (3x2 * 2x3 = 3x3)
    Eigen::Matrix<double, 3, 3> EC_A_T = E_mat * C_A.transpose();
    Eigen::Matrix<double, 3, 3> EC_B_T = E_mat * C_B.transpose();
    Eigen::Matrix<double, 3, 3> EC_C_T = E_mat * C_C.transpose();
    
    // Compute G * C_j * E^T for each node (3x3 * 3x2 * 2x3 = 3x3)
    Eigen::Matrix<double, 3, 3> GC_A_E_T = G * C_A * E_mat.transpose();
    Eigen::Matrix<double, 3, 3> GC_B_E_T = G * C_B * E_mat.transpose();
    Eigen::Matrix<double, 3, 3> GC_C_E_T = G * C_C * E_mat.transpose();
    
    // Compute K_AA = VOL * EC_A_T * GC_A_E_T
    Eigen::Matrix<double, 3, 3> K_AA_3x3 = VOL * EC_A_T * GC_A_E_T;
    // Extract 2x2 blocks: K_AA occupies DOFs (0,1) x (0,1)
    K(0,0) = K_AA_3x3(0,0); K(0,1) = K_AA_3x3(0,1);
    K(1,0) = K_AA_3x3(1,0); K(1,1) = K_AA_3x3(1,1);
    
    // K_AB = VOL * EC_A_T * GC_B_E_T
    Eigen::Matrix<double, 3, 3> K_AB_3x3 = VOL * EC_A_T * GC_B_E_T;
    K(0,2) = K_AB_3x3(0,0); K(0,3) = K_AB_3x3(0,1);
    K(1,2) = K_AB_3x3(1,0); K(1,3) = K_AB_3x3(1,1);
    
    // K_AC = VOL * EC_A_T * GC_C_E_T
    Eigen::Matrix<double, 3, 3> K_AC_3x3 = VOL * EC_A_T * GC_C_E_T;
    K(0,4) = K_AC_3x3(0,0); K(0,5) = K_AC_3x3(0,1);
    K(1,4) = K_AC_3x3(1,0); K(1,5) = K_AC_3x3(1,1);
    
    // K_BA = VOL * EC_B_T * GC_A_E_T
    Eigen::Matrix<double, 3, 3> K_BA_3x3 = VOL * EC_B_T * GC_A_E_T;
    K(2,0) = K_BA_3x3(0,0); K(2,1) = K_BA_3x3(0,1);
    K(3,0) = K_BA_3x3(1,0); K(3,1) = K_BA_3x3(1,1);
    
    // K_BB = VOL * EC_B_T * GC_B_E_T
    Eigen::Matrix<double, 3, 3> K_BB_3x3 = VOL * EC_B_T * GC_B_E_T;
    K(2,2) = K_BB_3x3(0,0); K(2,3) = K_BB_3x3(0,1);
    K(3,2) = K_BB_3x3(1,0); K(3,3) = K_BB_3x3(1,1);
    
    // K_BC = VOL * EC_B_T * GC_C_E_T
    Eigen::Matrix<double, 3, 3> K_BC_3x3 = VOL * EC_B_T * GC_C_E_T;
    K(2,4) = K_BC_3x3(0,0); K(2,5) = K_BC_3x3(0,1);
    K(3,4) = K_BC_3x3(1,0); K(3,5) = K_BC_3x3(1,1);
    
    // K_CA = VOL * EC_C_T * GC_A_E_T
    Eigen::Matrix<double, 3, 3> K_CA_3x3 = VOL * EC_C_T * GC_A_E_T;
    K(4,0) = K_CA_3x3(0,0); K(4,1) = K_CA_3x3(0,1);
    K(5,0) = K_CA_3x3(1,0); K(5,1) = K_CA_3x3(1,1);
    
    // K_CB = VOL * EC_C_T * GC_B_E_T
    Eigen::Matrix<double, 3, 3> K_CB_3x3 = VOL * EC_C_T * GC_B_E_T;
    K(4,2) = K_CB_3x3(0,0); K(4,3) = K_CB_3x3(0,1);
    K(5,2) = K_CB_3x3(1,0); K(5,3) = K_CB_3x3(1,1);
    
    // K_CC = VOL * EC_C_T * GC_C_E_T
    Eigen::Matrix<double, 3, 3> K_CC_3x3 = VOL * EC_C_T * GC_C_E_T;
    K(4,4) = K_CC_3x3(0,0); K(4,5) = K_CC_3x3(0,1);
    K(5,4) = K_CC_3x3(1,0); K(5,5) = K_CC_3x3(1,1);
    
    // Now output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << double_to_scientific(K(i,j));
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}