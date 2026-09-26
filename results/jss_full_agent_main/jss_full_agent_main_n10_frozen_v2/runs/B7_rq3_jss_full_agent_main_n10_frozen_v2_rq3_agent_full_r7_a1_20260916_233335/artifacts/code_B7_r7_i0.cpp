#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0 * (1.0 + nu));

    // H is 4x4, column-major in Fortran:
    // H(1..4)=col1=[1,x1,y1,z1], H(5..8)=col2=[1,x2,y2,z2], etc.
    Matrix4d H;
    H << 1.0, 1.0, 1.0, 1.0,
         0.0, 1.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0,
         0.0, 0.0, 0.0, 1.0;

    double hdeter = H.determinant();
    Matrix4d Hinv = H.inverse();
    double HDETER = std::fabs(hdeter) / 6.0;

    // Isotropic material stress-strain matrix G (6x6)
    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double g1 = E * (1.0 - nu) / TEMP1;
    double g2 = E * nu / TEMP1;
    G(0,0) = g1; G(1,0) = g2; G(2,0) = g2;
    G(0,1) = g2; G(1,1) = g1; G(2,1) = g2;
    G(0,2) = g2; G(1,2) = g2; G(2,2) = g1;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    // Build the strain-displacement matrices C_i (6x3) for each node.
    // Shape function derivatives are rows of Hinv:
    //   dN_i/dx = Hinv(i,1), dN_i/dy = Hinv(i,2), dN_i/dz = Hinv(i,3)
    MatrixXd Cfull = MatrixXd::Zero(6, 12);
    for (int I = 1; I <= 4; I++) {
        double dNdx = Hinv(I-1, 1);
        double dNdy = Hinv(I-1, 2);
        double dNdz = Hinv(I-1, 3);
        int c0 = 3*(I-1);
        Cfull(0, c0+0) = dNdx;
        Cfull(1, c0+1) = dNdy;
        Cfull(2, c0+2) = dNdz;
        Cfull(3, c0+0) = dNdy; Cfull(3, c0+1) = dNdx;
        Cfull(4, c0+1) = dNdz; Cfull(4, c0+2) = dNdy;
        Cfull(5, c0+0) = dNdz; Cfull(5, c0+2) = dNdx;
    }

    // Assemble 12x12 stiffness: K_ij = C_i^T * G * C_j * V
    MatrixXd K = MatrixXd::Zero(12, 12);
    for (int i = 0; i < 4; i++) {
        MatrixXd Ci = Cfull.block<6,3>(0, 3*i);
        for (int j = 0; j < 4; j++) {
            MatrixXd Cj = Cfull.block<6,3>(0, 3*j);
            Matrix3d Kij = Ci.transpose() * G * Cj * HDETER;
            K.block<3,3>(3*i, 3*j) = Kij;
        }
    }

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