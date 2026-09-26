#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0 * (1.0 + nu));

    // H matrix: column 1 = [1,X1,Y1,Z1], column 2 = [1,X2,Y2,Z2], etc.
    // (Fortran stores H(1..16) column-major as [1,X1,Y1,Z1, 1,X2,Y2,Z2, ...])
    Matrix4d H;
    H << 1.0, 1.0, 1.0, 1.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;

    double HDETER = std::abs(H.determinant()) / 6.0;
    Matrix4d Hinv = H.inverse();

    // Isotropic 6x6 material stress-strain matrix G
    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double g1 = E * (1.0 - nu) / TEMP1;
    double g2 = E * nu / TEMP1;
    G(0,0) = g1; G(1,1) = g1; G(2,2) = g1;
    G(0,1) = g2; G(0,2) = g2;
    G(1,0) = g2; G(1,2) = g2;
    G(2,0) = g2; G(2,1) = g2;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    // Standard strain-displacement B matrix (6x3) per node.
    // Strain order: [exx, eyy, ezz, gxy, gyz, gzx]; DOF order: (ux, uy, uz).
    // Shape-function derivatives come from columns 2,3,4 of Hinv.
    Matrix<double,6,3> B[4];
    for (int i = 0; i < 4; i++) {
        B[i].setZero();
        double dx = Hinv(i, 1);  // dN/dx
        double dy = Hinv(i, 2);  // dN/dy
        double dz = Hinv(i, 3);  // dN/dz
        B[i](0,0) = dx;
        B[i](1,1) = dy;
        B[i](2,2) = dz;
        B[i](3,0) = dy; B[i](3,1) = dx;
        B[i](4,1) = dz; B[i](4,2) = dy;
        B[i](5,0) = dz; B[i](5,2) = dx;
    }

    // Assemble full 12x12 stiffness: K_ij = HDETER * B_i^T * G * B_j
    Matrix<double,12,12> K = Matrix<double,12,12>::Zero();
    for (int i = 0; i < 4; i++) {
        Matrix<double,3,6> GCT = B[i].transpose() * G * HDETER;
        for (int j = 0; j < 4; j++) {
            K.block<3,3>(3*i, 3*j) = GCT * B[j];
        }
    }

    // Output as JSON
    std::cout << "{\"stiffness_matrix\":[";
    for (int r = 0; r < 12; r++) {
        std::cout << "[";
        for (int c = 0; c < 12; c++) {
            std::cout << std::scientific << std::setprecision(6) << K(r,c);
            if (c < 11) std::cout << ",";
        }
        std::cout << "]";
        if (r < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}