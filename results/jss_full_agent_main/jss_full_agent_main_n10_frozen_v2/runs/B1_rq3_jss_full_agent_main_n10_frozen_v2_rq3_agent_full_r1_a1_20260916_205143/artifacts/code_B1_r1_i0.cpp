#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    double E = 200e9, A = 0.01, G = 76.923e9, J = 5e-6;
    Eigen::Vector3d pA(0,0,0), pB(2,0,0);
    
    Eigen::Vector3d d = pA - pB;
    double XL = d.norm();
    Eigen::Vector3d XN = d / XL;
    
    double DSCL = A * E / XL;
    double DSCR = J * G / XL;
    
    Eigen::Matrix3d N = XN * XN.transpose();
    
    Eigen::Matrix<double,6,6> Kblk = Eigen::Matrix<double,6,6>::Zero();
    Kblk.block<3,3>(0,0) = DSCL * N;
    // no torsion
    
    Eigen::Matrix<double,12,12> K = Eigen::Matrix<double,12,12>::Zero();
    K.block<6,6>(0,0) =  Kblk;
    K.block<6,6>(6,6) =  Kblk;
    K.block<6,6>(0,6) = -Kblk;
    K.block<6,6>(6,0) = -Kblk;
    
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            double v = K(i,j);
            if (v == 0.0) {
                std::cout << "0.000000e+00";
            } else {
                std::cout << std::scientific << std::setprecision(6) << v;
            }
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}