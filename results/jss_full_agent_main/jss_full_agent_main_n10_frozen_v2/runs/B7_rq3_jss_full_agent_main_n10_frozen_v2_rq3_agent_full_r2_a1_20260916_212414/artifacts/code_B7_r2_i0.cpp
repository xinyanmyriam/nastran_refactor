#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double x[4] = {0.0, 1.0, 0.0, 0.0};
    double y[4] = {0.0, 0.0, 1.0, 0.0};
    double z[4] = {0.0, 0.0, 0.0, 1.0};

    double E  = 200e9;
    double nu = 0.3;
    double GG = E / (2.0 * (1.0 + nu));

    Matrix4d H;
    for (int i = 0; i < 4; ++i) {
        H(i, 0) = 1.0;
        H(i, 1) = x[i];
        H(i, 2) = y[i];
        H(i, 3) = z[i];
    }

    double HDETER = std::fabs(H.determinant());
    Matrix4d Hinv = H.inverse();

    Matrix<double, 6, 6> G = Matrix<double, 6, 6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double g1 = E * (1.0 - nu) / TEMP1;
    double g2 = E * nu / TEMP1;

    G(0, 0) = g1; G(1, 1) = g1; G(2, 2) = g1;
    G(0, 1) = g2; G(1, 0) = g2;
    G(0, 2) = g2; G(2, 0) = g2;
    G(1, 2) = g2; G(2, 1) = g2;
    G(3, 3) = GG; G(4, 4) = GG; G(5, 5) = GG;

    // C matrices use Hinv columns 2,3,4 (0-indexed 1,2,3)
    Matrix<double, 6, 3> Cmat[4];
    for (int I = 0; I < 4; ++I) {
        Cmat[I] = Matrix<double, 6, 3>::Zero();
        double hI4  = Hinv(I, 1);  // H(I+4)  -> col 2
        double hI8  = Hinv(I, 2);  // H(I+8)  -> col 3
        double hI12 = Hinv(I, 3);  // H(I+12) -> col 4

        Cmat[I](0, 0) = hI4;
        Cmat[I](4, 0) = hI8;
        Cmat[I](2, 1) = hI12;
        Cmat[I](4, 1) = hI12;
        Cmat[I](5, 1) = hI8;
        Cmat[I](0, 2) = hI12;
        Cmat[I](2, 2) = hI4;
        Cmat[I](3, 2) = hI8;
        Cmat[I](4, 2) = hI4;
    }

    HDETER = HDETER / 6.0;

    Matrix<double, 12, 12> K = Matrix<double, 12, 12>::Zero();
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Matrix3d Kij = Cmat[i].transpose() * G * Cmat[j] * HDETER;
            K.block<3,3>(3*i, 3*j) = Kij;
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << std::scientific << std::setprecision(6) << K(i, j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}