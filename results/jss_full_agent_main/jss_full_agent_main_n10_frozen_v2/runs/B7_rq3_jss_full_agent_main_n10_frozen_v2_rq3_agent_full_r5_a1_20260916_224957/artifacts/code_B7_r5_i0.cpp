#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double x[4] = {0.0, 1.0, 0.0, 0.0};
    double y[4] = {0.0, 0.0, 1.0, 0.0};
    double z[4] = {0.0, 0.0, 0.0, 1.0};

    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0 * (1.0 + nu));

    Matrix4d H;
    for (int i = 0; i < 4; i++) {
        H(i, 0) = 1.0;
        H(i, 1) = x[i];
        H(i, 2) = y[i];
        H(i, 3) = z[i];
    }

    Matrix4d Hinv = H.inverse();
    double HDETER = H.determinant();

    Matrix<double, 6, 6> G = Matrix<double, 6, 6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double g1 = E * (1.0 - nu) / TEMP1;
    double g2 = E * nu / TEMP1;
    G(0,0) = g1; G(1,1) = g1; G(2,2) = g1;
    G(0,1) = g2; G(0,2) = g2;
    G(1,0) = g2; G(1,2) = g2;
    G(2,0) = g2; G(2,1) = g2;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    Matrix<double, 6, 3> C[4];
    for (int i = 0; i < 4; i++) {
        C[i].setZero();
        double h2 = Hinv(i, 1);
        double h3 = Hinv(i, 2);
        double h4 = Hinv(i, 3);
        C[i](0, 0) = h2;
        C[i](4, 0) = h3;
        C[i](2, 1) = h4;
        C[i](4, 1) = h4;
        C[i](5, 1) = h3;
        C[i](0, 2) = h4;
        C[i](2, 2) = h2;
        C[i](3, 2) = h3;
        C[i](4, 2) = h2;
    }

    HDETER = std::fabs(HDETER) / 6.0;

    MatrixXd K = MatrixXd::Zero(12, 12);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Matrix3d Kij = C[i].transpose() * G * C[j] * HDETER;
            K.block<3,3>(3*i, 3*j) = Kij;
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int r = 0; r < 12; r++) {
        std::cout << "[";
        for (int c = 0; c < 12; c++) {
            std::cout << std::scientific << std::setprecision(6) << K(r, c);
            if (c < 11) std::cout << ",";
        }
        std::cout << "]";
        if (r < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}