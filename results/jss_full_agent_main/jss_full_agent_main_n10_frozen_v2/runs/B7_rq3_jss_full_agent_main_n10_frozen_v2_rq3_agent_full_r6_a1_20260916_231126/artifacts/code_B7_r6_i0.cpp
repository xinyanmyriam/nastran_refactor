#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0 * (1.0 + nu));

    double bx[4] = {-1.0, 1.0, 0.0, 0.0};
    double by[4] = {-1.0, 0.0, 1.0, 0.0};
    double bz[4] = {-1.0, 0.0, 0.0, 1.0};

    Matrix<double, 6, 6> D = Matrix<double, 6, 6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double d1 = E * (1.0 - nu) / TEMP1;
    double d2 = E * nu / TEMP1;
    D(0,0) = d1; D(1,1) = d1; D(2,2) = d1;
    D(0,1) = d2; D(0,2) = d2;
    D(1,0) = d2; D(1,2) = d2;
    D(2,0) = d2; D(2,1) = d2;
    D(3,3) = GG; D(4,4) = GG; D(5,5) = GG;

    Matrix<double, 6, 3> B[4];
    for (int i = 0; i < 4; i++) {
        B[i].setZero();
        B[i](0,0) = bx[i];
        B[i](1,1) = by[i];
        B[i](2,2) = bz[i];
        B[i](3,0) = by[i]; B[i](3,1) = bx[i];
        B[i](4,1) = bz[i]; B[i](4,2) = by[i];
        B[i](5,0) = bz[i]; B[i](5,2) = bx[i];
    }

    double V = 1.0 / 6.0;

    MatrixXd K = MatrixXd::Zero(12, 12);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            Matrix3d Kij = V * B[i].transpose() * D * B[j];
            K.block<3,3>(i*3, j*3) = Kij;
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        std::cout << "[";
        for (int j = 0; j < 12; j++) {
            std::cout << std::scientific << std::setprecision(6) << K(i, j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}