#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    // Test case parameters
    double E = 200e9, A = 0.01, G = 76.923e9, J = 5e-6;
    Eigen::Vector3d pA(0,0,0), pB(2,0,0);
    
    // Direction vector and length
    Eigen::Vector3d d = pA - pB;
    double XL = d.norm();
    Eigen::Vector3d XN = d / XL;
    
    // Axial and torsional stiffness scalars
    double DSCL = A * E / XL;
    double DSCR = J * G / XL;
    
    // N = XN * XN^T (outer product)
    Eigen::Matrix3d N = XN * XN.transpose();
    
    // Build 6x6 stiffness matrix
    Eigen::Matrix<double,6,6> KE = Eigen::Matrix<double,6,6>::Zero();
    KE.block<3,3>(0,0) = DSCL * N;  // axial block
    KE.block<3,3>(3,0) = DSCR * N;  // torsional block
    
    // Output as JSON
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << KE(i,j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}