#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Helper function to compute determinant and inverse of 4x4 matrix
// Returns determinant, and stores inverse in inv_out (column-major order as Fortran)
double invert_4x4(const Eigen::Matrix<double, 4, 4, Eigen::ColMajor>& A, Eigen::Matrix<double, 4, 4, Eigen::ColMajor>& inv_out) {
    // Use Eigen's LU decomposition for robust inversion
    Eigen::FullPivLU<Eigen::Matrix<double, 4, 4, Eigen::ColMajor>> lu(A);
    double det = lu.determinant();
    
    if (std::abs(det) < 1e-15) {
        // Singular matrix
        inv_out.setZero();
        return 0.0;
    }
    
    inv_out = lu.inverse();
    return det;
}

int main() {
    // Test case: Tetrahedron nodes
    // N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    std::vector<std::vector<double>> coords = {
        {0.0, 0.0, 0.0}, // N1
        {1.0, 0.0, 0.0}, // N2
        {0.0, 1.0, 0.0}, // N3
        {0.0, 0.0, 1.0}  // N4
    };

    // Material properties
    const double E = 200e9;      // Pa
    const double nu = 0.3;

    // Simulate ECPT array (100 elements, double precision)
    // We only need indices 1-22 (1-based) -> 0-21 (0-based)
    std::vector<double> ECPT(100, 0.0);

    // Element ID = 1
    ECPT[0] = 1.0;
    // Material ID = 1
    ECPT[1] = 1.0;
    // Grid point IDs: 1,2,3,4
    ECPT[2] = 1.0;
    ECPT[3] = 2.0;
    ECPT[4] = 3.0;
    ECPT[5] = 4.0;
    // Coord sys ID for each node (all 0)
    ECPT[6] = 0.0; // N1 coord sys ID
    // Coordinates for N1
    ECPT[7] = coords[0][0]; // X1
    ECPT[8] = coords[0][1]; // Y1
    ECPT[9] = coords[0][2]; // Z1
    ECPT[10] = 0.0; // N2 coord sys ID
    ECPT[11] = coords[1][0]; // X2
    ECPT[12] = coords[1][1]; // Y2
    ECPT[13] = coords[1][2]; // Z2
    ECPT[14] = 0.0; // N3 coord sys ID
    ECPT[15] = coords[2][0]; // X3
    ECPT[16] = coords[2][1]; // Y3
    ECPT[17] = coords[2][2]; // Z3
    ECPT[18] = 0.0; // N4 coord sys ID
    ECPT[19] = coords[3][0]; // X4
    ECPT[20] = coords[3][1]; // Y4
    ECPT[21] = coords[3][2]; // Z4
    // Element temperature = 0.0
    ECPT[22] = 0.0;

    // Simulate COMMON blocks
    // SMA1DP: C(72), G(36), H(16), TEMP(12), T(9), CT(18), GCT(18), KIJ(36), HDETER, TEMP1, NGPT, DIREC, KOUNT, TVOL
    std::vector<double> C(72, 0.0);
    std::vector<double> G(36, 0.0);
    std::vector<double> H(16, 0.0);
    std::vector<double> TEMP(12, 0.0);
    std::vector<double> T(9, 0.0);
    std::vector<double> CT(18, 0.0);
    std::vector<double> GCT(18, 0.0);
    std::vector<double> KIJ(36, 0.0);
    double HDETER = 0.0;
    double TEMP1 = 0.0;
    int NGPT = 0;
    int DIREC = 0;
    int KOUNT = 0;
    double TVOL = 0.0;

    // MATOUT: E, GG, NU, RHO, ALPHA, TSUB0, GSUBE, SIGT, SIGC, SIGS
    double E_mat = 0.0;
    double GG = 0.0;
    double NU = 0.0;
    double RHO = 0.0;
    double ALPHA = 0.0;
    double TSUB0 = 0.0;
    double GSUBE = 0.0;
    double SIGT = 0.0;
    double SIGC = 0.0;
    double SIGS = 0.0;

    // SMA1CL: IOPT4, K4GGSW, NPVT, ISKP(17), NOGOO
    int IOPT4 = 0;
    int K4GGSW = 0;
    int NPVT = 0; // pivot point — we'll use node 1 (grid ID 1) as pivot
    std::vector<int> ISKP(17, 0);
    int NOGOO = 0;

    // SMA1IO: DUM1(10), IFKGG, DUM2(1), IF4GG, DUM3(23)
    int IFKGG = 0;
    int IF4GG = 0;

    // SYSTEM: SYSBUF, OUT, NOGO
    int OUT = 0;
    bool NOGO = false;

    // BLANK: SKIP(16), VOLUME, SURFAC
    double VOLUME = 1.0; // assume unit volume scaling
    double SURFAC = 0.0;

    // MATIN: MATID, INFLAG, ELTEMP
    int MATID = 0;
    int INFLAG = 0;
    double ELTEMP = 0.0;

    // HMTOUT: MATBUF(7)
    std::vector<double> MATBUF(7, 0.0);

    // Simulate the KTETRA logic step-by-step

    // Step 1: Fill H matrix (4x4)
    // H = [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    H[0] = 1.0;
    H[1] = ECPT[7];  // x1
    H[2] = ECPT[8];  // y1
    H[3] = ECPT[9];  // z1
    H[4] = 1.0;
    H[5] = ECPT[11]; // x2
    H[6] = ECPT[12]; // y2
    H[7] = ECPT[13]; // z2
    H[8] = 1.0;
    H[9] = ECPT[15]; // x3
    H[10] = ECPT[16]; // y3
    H[11] = ECPT[17]; // z3
    H[12] = 1.0;
    H[13] = ECPT[19]; // x4
    H[14] = ECPT[20]; // y4
    H[15] = ECPT[21]; // z4

    // Convert H to Eigen 4x4 matrix (column-major as Fortran)
    // Fortran stores H(1) to H(16) in column-major order for a 4x4 matrix.
    // So H(1)=H11, H(2)=H21, H(3)=H31, H(4)=H41, H(5)=H12, etc.
    Eigen::Matrix<double, 4, 4, Eigen::ColMajor> H_fortran;
    for (int j = 0; j < 4; ++j) { // column index
        for (int i = 0; i < 4; ++i) { // row index
            H_fortran(i, j) = H[j*4 + i];
        }
    }

    // Invert H and get determinant
    Eigen::Matrix<double, 4, 4, Eigen::ColMajor> H_inv;
    HDETER = invert_4x4(H_fortran, H_inv);

    // Check determinant sign for orientation
    if (HDETER == 0.0) {
        std::cerr << "Error: Degenerate tetrahedron\n";
        return 1;
    }

    // Set DIREC based on sign of determinant
    DIREC = (HDETER > 0.0) ? 1 : -1;

    // Take absolute value for volume computation
    HDETER = std::abs(HDETER);

    // Material data: E, nu
    E_mat = E;
    NU = nu;

    // Compute temp1 = (1+nu)*(1-2*nu)
    TEMP1 = (1.0 + NU) * (1.0 - 2.0 * NU);
    if (std::abs(TEMP1) < 1e-12) {
        std::cerr << "Error: Illegal nu value\n";
        return 1;
    }

    // Fill G matrix (6x6 stress-strain matrix for isotropic material)
    // G is stored in column-major order as a 36-element vector
    // Indices: 1-6,7-12,...,31-36 → 0-5,6-11,...,30-35
    const double factor = E_mat / TEMP1;
    G[0]  = factor * (1.0 - NU);   // G11
    G[6]  = G[0];                  // G22
    G[12] = G[0];                  // G33
    G[1]  = factor * NU;           // G12
    G[2]  = G[1];                  // G13
    G[5]  = G[1];                  // G21
    G[7]  = G[1];                  // G23
    G[10] = G[1];                  // G31
    G[11] = G[1];                  // G32
    // Shear terms: G44 = G55 = G66 = G/2 = E/(2*(1+nu))
    GG = E_mat / (2.0 * (1.0 + NU));
    G[21] = GG; // G44 (index 21 = 4th diagonal -> 3*6 + 3 = 21)
    G[28] = GG; // G55 (index 28 = 5th diagonal -> 4*6 + 4 = 28)
    G[35] = GG; // G66 (index 35 = 6th diagonal -> 5*6 + 5 = 35)

    // Fill C matrices (4 matrices of size 6x3 each, stored consecutively in C(72))
    // But we don't actually use C for stiffness — we build B directly.

    // Compute the gradients from H_inv:
    // Since H = [1 x y z; ...], then H_inv * [1;x;y;z] = [N1;N2;N3;N4]
    // So the shape function gradients are in rows 1-3 of H_inv (0-indexed rows 1,2,3).
    Eigen::Matrix<double, 3, 4> dN_dx_dy_dz;
    dN_dx_dy_dz.row(0) = H_inv.row(1).transpose(); // dN_i/dx
    dN_dx_dy_dz.row(1) = H_inv.row(2).transpose(); // dN_i/dy
    dN_dx_dy_dz.row(2) = H_inv.row(3).transpose(); // dN_i/dz

    // Now build the full 6x12 B matrix (for 4 nodes * 3 DOF)
    Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();

    for (int i = 0; i < 4; ++i) {
        double dNdx = dN_dx_dy_dz(0, i);
        double dNdy = dN_dx_dy_dz(1, i);
        double dNdz = dN_dx_dy_dz(2, i);

        // u_x row: dN_i/dx at col 3*i
        B(0, 3*i + 0) = dNdx;
        // u_y row: dN_i/dy at col 3*i+1
        B(1, 3*i + 1) = dNdy;
        // u_z row: dN_i/dz at col 3*i+2
        B(2, 3*i + 2) = dNdz;
        // gamma_xy = du_y/dx + du_x/dy
        B(3, 3*i + 0) = dNdy; // du_x/dy
        B(3, 3*i + 1) = dNdx; // du_y/dx
        // gamma_xz = du_z/dx + du_x/dz
        B(4, 3*i + 0) = dNdz; // du_x/dz
        B(4, 3*i + 2) = dNdx; // du_z/dx
        // gamma_yz = du_z/dy + du_y/dz
        B(5, 3*i + 1) = dNdz; // du_y/dz
        B(5, 3*i + 2) = dNdy; // du_z/dy
    }

    // Volume factor: V = |det(H)| / 6
    double volume = HDETER / 6.0;

    // Construct 6x6 material matrix D (same as G in NASTRAN)
    Eigen::Matrix<double, 6, 6> D = Eigen::Matrix<double, 6, 6>::Zero();
    D(0,0) = G[0];  D(0,1) = G[1];  D(0,2) = G[2];
    D(1,0) = G[5];  D(1,1) = G[6];  D(1,2) = G[7];
    D(2,0) = G[10]; D(2,1) = G[11]; D(2,2) = G[12];
    D(3,3) = G[21];
    D(4,4) = G[28];
    D(5,5) = G[35];

    // Symmetrize D (it's symmetric)
    D(1,0) = D(0,1);
    D(2,0) = D(0,2);
    D(2,1) = D(1,2);

    // Compute element stiffness matrix: K = B^T * D * B * volume
    Eigen::Matrix<double, 12, 12> K = volume * B.transpose() * D * B;

    // Output stiffness matrix in plain text, row by row, space separated, scientific notation
    std::cout << "stiffness_matrix:\n";
    std::cout << std::setprecision(12) << std::scientific;
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            std::cout << K(i, j);
            if (j < 11) std::cout << " ";
        }
        std::cout << "\n";
    }

    return 0;
}