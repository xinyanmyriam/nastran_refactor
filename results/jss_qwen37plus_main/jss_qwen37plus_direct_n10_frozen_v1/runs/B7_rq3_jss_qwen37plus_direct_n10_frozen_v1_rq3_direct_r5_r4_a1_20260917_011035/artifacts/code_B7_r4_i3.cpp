#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// Helper function to compute determinant and inverse of 4x4 matrix
// Returns determinant, and stores inverse in inv_out (column-major order as Fortran expects)
double invert_4x4(const Eigen::Matrix<double, 4, 4>& A, Eigen::Matrix<double, 4, 4>& inv_out) {
    // Use Eigen's LU decomposition for robust inversion
    Eigen::FullPivLU<Eigen::Matrix4d> lu(A);
    double det = lu.determinant();
    
    if (std::abs(det) < 1e-15) {
        // Singular matrix
        inv_out.setZero();
        return 0.0;
    }
    
    inv_out = lu.inverse();
    return det;
}

// Helper: multiply 3x6 matrix C by 6x6 matrix G -> result is 3x6
Eigen::Matrix<double, 3, 6> multiply_CG(const Eigen::Matrix<double, 6, 3>& C, const Eigen::Matrix<double, 6, 6>& G) {
    // C is 6x3, G is 6x6 => C^T * G would be 3x6, but Fortran code uses C as 6x3 and multiplies C * G? 
    // Looking at Fortran: GMMATD(C(NPOINT),6,3,1, G(1),6,6,0, GCT(1)) with TRAN=1 on first arg means transpose C
    // So it computes C^T * G -> 3x6
    return C.transpose() * G;
}

// Helper: multiply 3x6 matrix GCT by 6x3 matrix CT -> result is 3x3
Eigen::Matrix3d multiply_GCT_CT(const Eigen::Matrix<double, 3, 6>& GCT, const Eigen::Matrix<double, 6, 3>& CT) {
    // GCT is 3x6, CT is 6x3 => GCT * CT is 3x3
    return GCT * CT;
}

// Helper: multiply 3x6 matrix GCT by 6x3 matrix C -> result is 3x3
Eigen::Matrix3d multiply_GCT_C(const Eigen::Matrix<double, 3, 6>& GCT, const Eigen::Matrix<double, 6, 3>& C) {
    return GCT * C;
}

int main() {
    // Test case: Tetrahedron N1=(0,0,0), N2=(1,0,0), N3=(0,1,0), N4=(0,0,1)
    // E = 200e9, nu = 0.3
    const double E = 200e9;
    const double nu = 0.3;

    // ECPT array simulation (Fortran indexing starts at 1, we use 0-based)
    // ECPT(1) = element ID = 1
    // ECPT(2) = material ID = 1
    // ECPT(3) = grid point 1 = 1
    // ECPT(4) = grid point 2 = 2
    // ECPT(5) = grid point 3 = 3
    // ECPT(6) = grid point 4 = 4
    // ECPT(7) = coord sys ID for grid 1 = 0 (global)
    // ECPT(8-10) = x1,y1,z1 = 0,0,0
    // ECPT(11) = coord sys ID for grid 2 = 0
    // ECPT(12-14) = x2,y2,z2 = 1,0,0
    // ECPT(15) = coord sys ID for grid 3 = 0
    // ECPT(16-18) = x3,y3,z3 = 0,1,0
    // ECPT(19) = coord sys ID for grid 4 = 0
    // ECPT(20-22) = x4,y4,z4 = 0,0,1
    // ECPT(23) = element temperature = 0.0

    std::vector<double> ECPT(100, 0.0);
    ECPT[0] = 1.0;   // element ID
    ECPT[1] = 1.0;   // material ID
    ECPT[2] = 1.0;   // grid point 1
    ECPT[3] = 2.0;   // grid point 2
    ECPT[4] = 3.0;   // grid point 3
    ECPT[5] = 4.0;   // grid point 4
    ECPT[6] = 0.0;   // coord sys ID grid 1
    ECPT[7] = 0.0;   // x1
    ECPT[8] = 0.0;   // y1
    ECPT[9] = 0.0;   // z1
    ECPT[10] = 0.0;  // coord sys ID grid 2
    ECPT[11] = 1.0;  // x2
    ECPT[12] = 0.0;  // y2
    ECPT[13] = 0.0;  // z2
    ECPT[14] = 0.0;  // coord sys ID grid 3
    ECPT[15] = 0.0;  // x3
    ECPT[16] = 1.0;  // y3
    ECPT[17] = 0.0;  // z3
    ECPT[18] = 0.0;  // coord sys ID grid 4
    ECPT[19] = 0.0;  // x4
    ECPT[20] = 0.0;  // y4
    ECPT[21] = 1.0;  // z4
    ECPT[22] = 0.0;  // element temperature

    // Simulate common blocks
    double VOLUME = 1.0;  // volume scaling factor
    double SURFAC = 0.0; // surface area scaling factor
    int NGPT = 4;         // number of integration points for tetra
    double HDETER = 0.0;  // determinant of H matrix
    int DIREC = 0;        // orientation flag
    int KOUNT = 0;        // counter
    double TVOL = 0.0;    // total volume accumulator
    int NPVT = 0;         // pivot point (not used in basic case)
    int IOPT = 0;         // option flag
    int JTYPE = 0;        // element type (0 for tetra)

    // Build H matrix (4x4) as described in Fortran
    // H = [1 x1 y1 z1; 1 x2 y2 z2; 1 x3 y3 z3; 1 x4 y4 z4]
    Eigen::Matrix<double, 4, 4> H;
    H << 1.0, ECPT[7], ECPT[8], ECPT[9],
         1.0, ECPT[11], ECPT[12], ECPT[13],
         1.0, ECPT[15], ECPT[16], ECPT[17],
         1.0, ECPT[19], ECPT[20], ECPT[21];

    // Invert H and get determinant
    Eigen::Matrix<double, 4, 4> H_inv;
    HDETER = invert_4x4(H, H_inv);

    // Check for bad geometry
    if (std::abs(HDETER) < 1e-15) {
        std::cerr << "Error: Bad or reverse geometry\n";
        return 1;
    }

    // Set direction flag based on determinant sign
    if (DIREC == 0) {
        DIREC = (HDETER > 0.0) ? 1 : -1;
    }

    // Material properties
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    if (std::abs(TEMP1) < 1e-6) {
        std::cerr << "Error: Illegal value of nu\n";
        return 1;
    }

    // Build 6x6 G matrix (stress-strain matrix for isotropic material)
    // G = [C11 C12 C12 0   0   0;
    //      C12 C11 C12 0   0   0;
    //      C12 C12 C11 0   0   0;
    //      0   0   0   C44 0   0;
    //      0   0   0   0   C44 0;
    //      0   0   0   0   0   C44]
    // where C11 = E*(1-nu)/((1+nu)*(1-2*nu)), C12 = E*nu/((1+nu)*(1-2*nu)), C44 = E/(2*(1+nu))
    double C11 = E * (1.0 - nu) / TEMP1;
    double C12 = E * nu / TEMP1;
    double C44 = E / (2.0 * (1.0 + nu));

    Eigen::Matrix<double, 6, 6> G = Eigen::Matrix<double, 6, 6>::Zero();
    G(0,0) = C11; G(0,1) = C12; G(0,2) = C12;
    G(1,0) = C12; G(1,1) = C11; G(1,2) = C12;
    G(2,0) = C12; G(2,1) = C12; G(2,2) = C11;
    G(3,3) = C44;
    G(4,4) = C44;
    G(5,5) = C44;

    // Build C matrices (4 matrices of size 6x3 each)
    // From Fortran: for i=1..4, C(i) = H(i+4), C(i+4) = H(i+8), C(i+8) = H(i+12)
    // Then set other entries as per the pattern
    std::vector<Eigen::Matrix<double, 6, 3>> C_matrices(4);
    for (int i = 0; i < 4; ++i) {
        // Get row i of H (0-indexed), which corresponds to H(i+1) in Fortran
        // H rows: row0 = [1,x1,y1,z1], row1 = [1,x2,y2,z2], etc.
        double x = H(i,1);
        double y = H(i,2);
        double z = H(i,3);

        // Fill C matrix for this node
        // Pattern from Fortran: C(J+1)=H(I+4), C(J+5)=H(I+8), C(J+9)=H(I+12),
        // C(J+11)=H(I+12), C(J+12)=H(I+8), C(J+13)=H(I+12), C(J+15)=H(I+4), C(J+16)=H(I+8), C(J+17)=H(I+4)
        // But looking at indices, it's building derivative matrices
        // Actually, for tetrahedral element, the B-matrix derivatives are constant:
        // dN1/dx = h11, dN1/dy = h12, dN1/dz = h13, etc.
        // The standard approach is: B = [dN1/dx 0 0 dN2/dx 0 0 dN3/dx 0 0 dN4/dx 0 0;
        //                              0 dN1/dy 0 0 dN2/dy 0 0 dN3/dy 0 0 dN4/dy 0;
        //                              0 0 dN1/dz 0 0 dN2/dz 0 0 dN3/dz 0 0 dN4/dz;
        //                              dN1/dy dN1/dx 0 dN2/dy dN2/dx 0 dN3/dy dN3/dx 0 dN4/dy dN4/dx 0;
        //                              0 dN1/dz dN1/dy 0 dN2/dz dN2/dy 0 dN3/dz dN3/dy 0 dN4/dz dN4/dy;
        //                              dN1/dz 0 dN1/dx dN2/dz 0 dN2/dx dN3/dz 0 dN3/dx dN4/dz 0 dN4/dx]
        // Where dNi/dx, dNi/dy, dNi/dz come from H_inv (since N = H_inv * [1;x;y;z])
        
        // The gradient of shape functions is given by the last 3 columns of H_inv (excluding first column)
        // Since N = H_inv * [1; x; y; z], then dN/dx = H_inv(:,1), dN/dy = H_inv(:,2), dN/dz = H_inv(:,3)
        Eigen::Vector4d dNdx = H_inv.col(1);
        Eigen::Vector4d dNdy = H_inv.col(2);
        Eigen::Vector4d dNdz = H_inv.col(3);

        // Build 6x3 C matrix for node i (which corresponds to derivatives for shape function i)
        // But Fortran builds 4 separate 6x3 matrices, one for each node
        // From the pattern in Fortran, it appears they are building:
        // C_i = [dNi/dx, 0, 0;
        //        0, dNi/dy, 0;
        //        0, 0, dNi/dz;
        //        dNi/dy, dNi/dx, 0;
        //        0, dNi/dz, dNi/dy;
        //        dNi/dz, 0, dNi/dx]
        C_matrices[i] << dNdx(i), 0.0, 0.0,
                           0.0, dNdy(i), 0.0,
                           0.0, 0.0, dNdz(i),
                           dNdy(i), dNdx(i), 0.0,
                           0.0, dNdz(i), dNdy(i),
                           dNdz(i), 0.0, dNdx(i);
    }

    // Compute HDETER scaling factor
    // From Fortran: divide by 6.0 for tetra
    double scale_factor = std::abs(HDETER) / 6.0;

    // Initialize global stiffness matrix (12x12)
    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();

    // Process each of the 4 nodes
    for (int i = 0; i < 4; ++i) {
        // Get C matrix for this node (6x3)
        const Eigen::Matrix<double, 6, 3>& C_i = C_matrices[i];

        // Compute GCT = C_i^T * G  (3x6)
        Eigen::Matrix<double, 3, 6> GCT = C_i.transpose() * G;

        // Scale by determinant factor
        GCT *= scale_factor;

        // Now compute contribution to stiffness: K_i = GCT * C_i (3x3)
        Eigen::Matrix3d K_i = GCT * C_i;

        // Map 3x3 local stiffness to global 12x12 stiffness
        // Node i has DOFs at positions [3*i, 3*i+1, 3*i+2] for ux, uy, uz
        int base = 3 * i;
        K.block<3,3>(base, base) += K_i;
    }

    // Output as plain formatted matrix (matching original JSON structure intent)
    // Print: {"stiffness_matrix": [[...], [...], ...]}
    std::cout << "{\"stiffness_matrix\": [" << std::endl;
    for (int i = 0; i < 12; ++i) {
        std::cout << "  [";
        for (int j = 0; j < 12; ++j) {
            std::cout << std::setprecision(12) << std::scientific << K(i,j);
            if (j < 11) std::cout << ", ";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
        std::cout << std::endl;
    }
    std::cout << "]}" << std::endl;

    return 0;
}