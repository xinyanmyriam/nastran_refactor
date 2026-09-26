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

    double hdeter = H.determinant();
    Matrix4d Hinv = H.inverse();
    hdeter = std::abs(hdeter);

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double temp1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double c1 = E * (1.0 - nu) / temp1;
    double c2 = E * nu / temp1;

    G(0,0) = c1; G(1,1) = c1; G(2,2) = c1;
    G(0,1) = c2; G(0,2) = c2;
    G(1,0) = c2; G(1,2) = c2;
    G(2,0) = c2; G(2,1) = c2;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    Matrix<double,6,3> C[4];
    for (int i = 0; i < 4; ++i) {
        C[i].setZero();
        double dNdx = Hinv(i, 1);
        double dNdy = Hinv(i, 2);
        double dNdz = Hinv(i, 3);

        C[i](0, 0) = dNdx;
        C[i](4, 0) = dNdz;
        C[i](2, 1) = dNdy;
        C[i](4, 1) = dNdz;
        C[i](5, 1) = dNdy;
        C[i](0, 2) = dNdx;
        C[i](2, 2) = dNdy;
        C[i](3, 2) = dNdz;
        C[i](4, 2) = dNdx;
    }

    double hdet = hdeter / 6.0;

    Matrix<double,12,12> K = Matrix<double,12,12>::Zero();

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            Matrix3d KIJ = C[i].transpose() * G * C[j] * hdet;
            K.block<3,3>(3*i, 3*j) = KIJ;
        }
    }

    std::cout << std::setprecision(10);
    std::cout << "{\"stiffness_matrix\":[";
    for (int r = 0; r < 12; ++r) {
        std::cout << "[";
        for (int c = 0; c < 12; ++c) {
            std::cout << std::scientific << K(r, c);
            if (c < 11) std::cout << ",";
        }
        std::cout << "]";
        if (r < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}