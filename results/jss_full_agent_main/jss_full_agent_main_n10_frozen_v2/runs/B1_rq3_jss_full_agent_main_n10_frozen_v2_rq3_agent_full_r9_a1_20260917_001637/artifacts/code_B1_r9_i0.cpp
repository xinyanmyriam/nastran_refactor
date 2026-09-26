#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    Vector3d nodeA(0.0, 0.0, 0.0);
    Vector3d nodeB(2.0, 0.0, 0.0);
    double E = 200e9;
    double A = 0.01;
    double G = 76.923e9;
    double J = 5e-6;

    Vector3d d = nodeA - nodeB;
    double XL = d.norm();
    Vector3d XN = d / XL;
    Matrix3d N = XN * XN.transpose();

    double DSCL = A * E / XL;
    double DSCR = J * G / XL;

    MatrixXd K = MatrixXd::Zero(12, 12);

    K.block<3,3>(0,0) = DSCL * N;
    K.block<3,3>(3,3) = DSCR * N;
    K.block<3,3>(6,6) = DSCL * N;
    K.block<3,3>(9,9) = DSCR * N;

    K.block<3,3>(0,6) = -DSCL * N;
    K.block<3,3>(6,0) = -DSCL * N;
    K.block<3,3>(3,9) = -DSCR * N;
    K.block<3,3>(9,3) = -DSCR * N;

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        std::cout << "[";
        for (int j = 0; j < 12; j++) {
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}