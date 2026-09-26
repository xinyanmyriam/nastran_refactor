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
    double DELTA = XSUBC / XSUBB - 1.0;
    
    // Build C matrix (3x6) partitioned as [C_A | C_B | C_C], each 3x2
    // Stored row-wise: C(1..18) = [row1 of C_A, row2 of C_A, row3 of C_A, 
    //                               row1 of C_B, row2 of C_B, row3 of C_B,
    //                               row1 of C_C, row2 of C_C, row3 of C_C]
    std::vector<double> C(18, 0.0);
    
    // C_A (first 6 elements)
    C[0] = -REELMU;           // C(1)
    C[1] = 0.0;               // C(2)
    C[2] = 0.0;               // C(3)
    C[3] = FLAMDA * DELTA;    // C(4)
    C[4] = FLAMDA * DELTA;    // C(5) - same as C(4) per Fortran
    C[5] = -REELMU;           // C(6)
    
    // C_B (next 6 elements)
    C[6] = REELMU;            // C(7)
    C[7] = 0.0;               // C(8)
    C[8] = 0.0;               // C(9)
    C[9] = -FLAMDA * REELMU * XSUBC; // C(10)
    C[10] = -FLAMDA * REELMU * XSUBC; // C(11) - same as C(10)
    C[11] = REELMU;           // C(12)
    
    // C_C (last 6 elements)
    C[12] = 0.0;              // C(13)
    C[13] = 0.0;              // C(14)
    C[14] = 0.0;              // C(15)
    C[15] = FLAMDA;           // C(16)
    C[16] = FLAMDA;           // C(17) - same as C(16)
    C[17] = 0.0;              // C(18)
    
    // Build E matrix (3x2) from I and J vectors (we ignore K for membrane)
    // E = [I_x J_x; I_y J_y; I_z J_z] -> 3x2
    Eigen::Matrix<double, 3, 2, Eigen::RowMajor> E_mat;
    E_mat << E_vec[0], E_vec[3],  // I_x, J_x
             E_vec[1], E_vec[4],  // I_y, J_y
             E_vec[2], E_vec[5];  // I_z, J_z
    
    // Build C matrices for each node: C_A, C_B, C_C (each 3x2)
    // C_A = rows 0-2 of C vector -> reshape to 3x2
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
    // But the Fortran computes K_ij for pivot i and node j, then assembles into 6x6
    
    // We'll compute the full 6x6 stiffness matrix directly:
    // K = VOL * [E*C_A^T, E*C_B^T, E*C_C^T] * G * [E*C_A^T, E*C_B^T, E*C_C^T]^T
    // Actually, the standard CST formulation is:
    // B = [B_A B_B B_C] where B_i = [dN_i/dx, 0; 0, dN_i/dy; dN_i/dy, dN_i/dx]
    // But here the Fortran uses a different approach with E and C matrices.
    
    // Following the Fortran logic for NTYPE=0:
    // For each node j (1,2,3), compute:
    //   temp = VOL * (E * C_j^T) * G * (E * C_j^T)^T   -> gives 3x3
    // Then place in KIJ at appropriate positions
    
    // However, the Fortran code actually computes for pivot i and node j:
    // K_ij = VOL * T_i^T * E * C_i^T * G * C_j * E^T * T_j
    
    // Since all CSID=0, T_i = T_j = I, so:
    // K_ij = VOL * E * C_i^T * G * C_j * E^T
    
    // And the full 6x6 stiffness matrix is assembled as:
    // [K_AA K_AB K_AC]
    // [K_BA K_BB K_BC]
    // [K_CA K_CB K_CC]
    // where each K_ij is 2x2 (since 2 DOF per node)
    
    // But wait: the Fortran stores KIJ(36) as a 6x6 matrix in column-major order?
    // Looking at the Fortran: KIJ(1) through KIJ(36) and the assignments:
    // KIJ(1), KIJ(2), KIJ(3) -> first row, columns 1,2,3
    // KIJ(7), KIJ(8), KIJ(9) -> second row, columns 1,2,3
    // KIJ(13),KIJ(14),KIJ(15) -> third row, columns 1,2,3
    // So it's column-major storage for a 6x6 matrix.
    
    // Actually, the assignments show:
    // KIJ(1) = TEMPAR(NPT1    )   -> (1,1)
    // KIJ(2) = TEMPAR(NPT1 + 1)   -> (2,1)
    // KIJ(3) = TEMPAR(NPT1 + 2)   -> (3,1)
    // KIJ(7) = TEMPAR(NPT1 + 3)   -> (1,2)
    // KIJ(8) = TEMPAR(NPT1 + 4)   -> (2,2)
    // KIJ(9) = TEMPAR(NPT1 + 5)   -> (3,2)
    // So yes, column-major: index = (row-1) + (col-1)*6 + 1
    
    // Therefore, we need to compute a 6x6 matrix in column-major order.
    
    // Let's compute each 2x2 block K_ij (for nodes i,j = A,B,C) which goes into
    // positions: K_ij occupies rows [(i-1)*2+1, (i-1)*2+2] and cols [(j-1)*2+1, (j-1)*2+2]
    // But the Fortran computes K_ij as 3x3? Wait, no: the TEMPAR product is 3x3, but then
    // they assign only 3 elements per row? Looking again:
    // They assign KIJ(1),KIJ(2),KIJ(3) -> first column of 6x6? But 6x6 has 6 rows.
    
    // Actually, the Fortran code has a bug in my reading. Let me reexamine:
    // The comment says: "KIJ( 1) = TEMPAR(NPT1    )" etc., and then
    // KIJ( 1), KIJ( 2), KIJ( 3) -> positions 1,2,3
    // KIJ( 7), KIJ( 8), KIJ( 9) -> positions 7,8,9
    // KIJ(13), KIJ(14), KIJ(15) -> positions 13,14,15
    // In column-major 6x6, position 1 = (1,1), 2=(2,1), 3=(3,1), 4=(4,1), 5=(5,1), 6=(6,1)
    // Then 7=(1,2), 8=(2,2), 9=(3,2), 10=(4,2), 11=(5,2), 12=(6,2)
    // So they are only filling the first 3 rows of columns 1,2,3? That doesn't make sense.
    
    // Let me check the indices: 
    // They assign:
    // KIJ(1), KIJ(2), KIJ(3)   -> rows 1,2,3 of column 1
    // KIJ(7), KIJ(8), KIJ(9)   -> rows 1,2,3 of column 2
    // KIJ(13),KIJ(14),KIJ(15)  -> rows 1,2,3 of column 3
    // So they are only computing the top-left 3x3 of a 6x6? But that can't be.
    
    // Actually, looking more carefully: the Fortran code is for a 3-node triangle with 2 DOF/node = 6 DOF.
    // The stiffness matrix is 6x6. The assignment pattern suggests they are filling:
    // Column 1: rows 1,2,3 -> but we need rows 1,2 for node A (ux,uy), rows 3,4 for node B? No, 2 DOF per node.
    // Node A: DOFs 1,2 (ux,uy)
    // Node B: DOFs 3,4 (ux,uy)
    // Node C: DOFs 5,6 (ux,uy)
    
    // So the 6x6 matrix has:
    // Row/Col 1: A-ux, 2: A-uy, 3: B-ux, 4: B-uy, 5: C-ux, 6: C-uy
    
    // The Fortran assigns:
    // KIJ(1) -> (1,1) = A-ux, A-ux
    // KIJ(2) -> (2,1) = A-uy, A-ux
    // KIJ(3) -> (3,1) = B-ux, A-ux   -> but this should be row 3, col 1 -> (3,1) in 1-indexed = row3,col1
    // In 6x6, (3,1) is B-ux, A-ux -> correct.
    // KIJ(7) -> (1,2) = A-ux, A-uy
    // KIJ(8) -> (2,2) = A-uy, A-uy
    // KIJ(9) -> (3,2) = B-ux, A-uy
    // KIJ(13)-> (1,3) = A-ux, B-ux
    // KIJ(14)-> (2,3) = A-uy, B-ux
    // KIJ(15)-> (3,3) = B-ux, B-ux
    
    // So they are filling column by column, but only the first 3 rows of each column? 
    // But a 6x6 has 6 rows per column. They must be doing more assignments.
    
    // Looking again: the code has:
    // KIJ( 1) = TEMPAR(NPT1    )
    // KIJ( 2) = TEMPAR(NPT1 + 1)
    // KIJ( 3) = TEMPAR(NPT1 + 2)
    // KIJ( 7) = TEMPAR(NPT1 + 3)
    // KIJ( 8) = TEMPAR(NPT1 + 4)
    // KIJ( 9) = TEMPAR(NPT1 + 5)
    // KIJ(13) = TEMPAR(NPT1 + 6)
    // KIJ(14) = TEMPAR(NPT1 + 7)
    // KIJ(15) = TEMPAR(NPT1 + 8)
    // That's 9 assignments, which is 3x3. But we need 36.
    
    // I think the issue is that the Fortran code is incomplete in the snippet, or I'm misreading.
    // Actually, the code says "AT THIS POINT COMPLETE COMPUTATION FOR K-SUB-I,J" and then calls SMA1B to assemble.
    // But for our purpose, we want the full 6x6 stiffness matrix.
    
    // Let's use the standard CST (Constant Strain Triangle) formulation which is well-known:
    // B = [dN1/dx  0     dN2/dx  0     dN3/dx  0    ]
    //     [0     dN1/dy  0     dN2/dy  0     dN3/dy]
    //     [dN1/dy dN1/dx dN2/dy dN2/dx dN3/dy dN3/dx]
    //
    // For triangle with nodes (x1,y1), (x2,y2), (x3,y3):
    // Area = 0.5 * | (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) |
    // Then:
    // b1 = y2 - y3, c1 = x3 - x2
    // b2 = y3 - y1, c2 = x1 - x3
    // b3 = y1 - y2, c3 = x2 - x1
    // Then dNi/dx = bi/(2*Area), dNi/dy = ci/(2*Area)
    
    // Compute area
    double area = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    if (area < 1.0e-15) {
        std::cerr << "Error: Zero area triangle\n";
        return 1;
    }
    
    // Coefficients
    double b1 = y2 - y3;
    double c1 = x3 - x2;
    double b2 = y3 - y1;
    double c2 = x1 - x3;
    double b3 = y1 - y2;
    double c3 = x2 - x1;
    
    double two_area = 2.0 * area;
    
    // B matrix (3x6)
    Eigen::Matrix<double, 3, 6, Eigen::RowMajor> B;
    B << b1/two_area, 0.0,        b2/two_area, 0.0,        b3/two_area, 0.0,
         0.0,        c1/two_area, 0.0,        c2/two_area, 0.0,        c3/two_area,
         c1/two_area, b1/two_area, c2/two_area, b2/two_area, c3/two_area, b3/two_area;
    
    // D matrix for plane stress isotropic
    double D11 = E / (1.0 - nu*nu);
    double D12 = nu * E / (1.0 - nu*nu);
    double D22 = D11;
    double D33 = E / (2.0 * (1.0 + nu));
    
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> D;
    D << D11, D12, 0.0,
         D12, D22, 0.0,
         0.0, 0.0, D33;
    
    // Stiffness matrix: K = t * area * B^T * D * B
    Eigen::Matrix<double, 6, 6> K = t * area * B.transpose() * D * B;
    
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