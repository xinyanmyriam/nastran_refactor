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

    Matrix4d Hinv = H.inverse();
    double HDETER = std::abs(H.determinant());

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0 + nu) * (1.0 - 2.0 * nu);
    double g1 = E * (1.0 - nu) / TEMP1;
    double g2 = E * nu / TEMP1;
    G(0,0) = g1; G(1,1) = g1; G(2,2) = g1;
    G(0,1) = g2; G(0,2) = g2;
    G(1,0) = g2; G(1,2) = g2;
    G(2,0) = g2; G(2,1) = g2;
    G(3,3) = GG; G(4,4) = GG; G(5,5) = GG;

    HDETER = HDETER / 6.0;

    MatrixXd C[4];
    for (int I = 1; I <= 4; I++) {
        double a = Hinv(I-1, 1);
        double b = Hinv(I-1, 2);
        double c = Hinv(I-1, 3);
        MatrixXd Ci = MatrixXd::Zero(6, 3);
        Ci(0,0) = a; Ci(0,2) = c;
        Ci(2,1) = c; Ci(2,2) = a;
        Ci(3,1) = c; Ci(3,2) = b;
        Ci(4,0) = b; Ci(4,1) = c; Ci(4,2) = a;
        Ci(5,1) = b;
        C[I-1] = Ci;
    }

    MatrixXd K = MatrixXd::Zero(12, 12);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            MatrixXd Kij = C[i].transpose() * G * C[j] * HDETER;
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