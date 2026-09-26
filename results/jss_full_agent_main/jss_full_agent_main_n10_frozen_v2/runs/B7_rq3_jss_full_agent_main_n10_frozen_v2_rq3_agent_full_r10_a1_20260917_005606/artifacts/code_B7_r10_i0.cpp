#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0 * (1.0 + nu));

    Matrix4d H;
    H << 1.0, 0.0, 0.0, 0.0,
         1.0, 1.0, 0.0, 0.0,
         1.0, 0.0, 1.0, 0.0,
         1.0, 0.0, 0.0, 1.0;

    double HDETER = H.determinant();
    Matrix4d Hinv = H.inverse();

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double g1 = E * (1.0 - nu) / TEMP1;
    double g2 = E * nu / TEMP1;
    G(0,0) = g1; G(1,1) = g1; G(2,2) = g1;
    G(0,1) = g2; G(0,2) = g2;
    G(1,0) = g2; G(1,2) = g2;
    G(2,0) = g2; G(2,1) = g2;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    double Hf[17];
    for (int c = 0; c < 4; c++)
        for (int r = 0; r < 4; r++)
            Hf[c*4 + r + 1] = Hinv(r, c);

    Matrix<double,6,3> C[4];
    for (int i = 0; i < 4; i++) {
        C[i] = Matrix<double,6,3>::Zero();
        int I = i + 1;
        C[i](0,0) = Hf[I+4];
        C[i](4,0) = Hf[I+8];
        C[i](2,1) = Hf[I+12];
        C[i](4,1) = Hf[I+12];
        C[i](5,1) = Hf[I+8];
        C[i](0,2) = Hf[I+12];
        C[i](2,2) = Hf[I+4];
        C[i](3,2) = Hf[I+8];
        C[i](4,2) = Hf[I+4];
    }

    HDETER = HDETER / 6.0;

    Matrix<double,12,12> K = Matrix<double,12,12>::Zero();
    for (int I = 0; I < 4; I++) {
        for (int J = 0; J < 4; J++) {
            Matrix<double,3,3> kij = C[I].transpose() * G * C[J] * HDETER;
            K.block<3,3>(3*I, 3*J) = kij;
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