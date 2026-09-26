#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    Eigen::Vector3d A(0.0, 0.0, 0.0);
    Eigen::Vector3d B(2.0, 0.0, 0.0);
    double E = 200e9;
    double Area = 0.01;
    double G = 76.923e9;
    double J = 5e-6;

    Eigen::Vector3d d = A - B;
    double XL = d.norm();
    Eigen::Vector3d XN = d / XL;

    double DSCL = Area * E / XL;
    double DSCR = J * G / XL;

    Eigen::Matrix3d N = XN * XN.transpose();

    Eigen::Matrix<double, 12, 12> K = Eigen::Matrix<double, 12, 12>::Zero();
    Eigen::Matrix3d Kt = DSCL * N;
    Eigen::Matrix3d Kr = DSCR * N;

    K.block<3,3>(0,0) = Kt;
    K.block<3,3>(3,3) = Kr;
    K.block<3,3>(6,6) = Kt;
    K.block<3,3>(9,9) = Kr;
    K.block<3,3>(0,6) = -Kt;
    K.block<3,3>(6,0) = -Kt;
    K.block<3,3>(3,9) = -Kr;
    K.block<3,3>(9,3) = -Kr;

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}