#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    double E = 200e9, A = 0.01, G = 76.923e9, J = 5e-6;
    Eigen::Vector3d pA(0,0,0), pB(2,0,0);
    Eigen::Vector3d d = pA - pB;
    double XL = d.norm();
    Eigen::Vector3d XN = d / XL;
    double DSCL = A*E/XL, DSCR = J*G/XL;
    Eigen::Matrix3d N = XN * XN.transpose();
    Eigen::MatrixXd K6 = Eigen::MatrixXd::Zero(6,6);
    K6.block<3,3>(0,0) = DSCL*N;
    K6.block<3,3>(3,0) = DSCR*N;
    Eigen::MatrixXd K = Eigen::MatrixXd::Zero(12,12);
    K.block<6,6>(0,0)=K6; K.block<6,6>(0,6)=-K6;
    K.block<6,6>(6,0)=-K6; K.block<6,6>(6,6)=K6;
    std::cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<12;i++) for(int j=0;j<12;j++){
        std::cout << std::scientific << std::setprecision(6) << K(i,j);
        if(!(i==11&&j==11)) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}