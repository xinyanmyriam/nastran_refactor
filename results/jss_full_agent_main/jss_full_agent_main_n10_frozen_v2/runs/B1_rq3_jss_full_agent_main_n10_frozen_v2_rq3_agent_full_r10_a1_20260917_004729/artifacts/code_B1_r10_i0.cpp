#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

int main() {
    // Test case: Node A (0,0,0), Node B (2,0,0)
    Eigen::Vector3d A(0.0, 0.0, 0.0);
    Eigen::Vector3d B(2.0, 0.0, 0.0);
    double E = 200e9;      // Young's modulus (Pa)
    double Area = 0.01;    // Cross-sectional area (m^2)
    double G = 76.923e9;   // Shear modulus (Pa)
    double J = 5e-6;       // Polar moment of inertia (m^4)

    // Direction vector from A to B (matches Fortran: X = D(1)-D(4))
    Eigen::Vector3d d = A - B;
    double XL = d.norm();

    // Normalized direction vector
    Eigen::Vector3d XN = d / XL;

    // Axial and torsional stiffness coefficients
    double DSCL = Area * E / XL;
    double DSCR = J * G / XL;

    // N matrix = outer product of direction vector
    Eigen::Matrix3d N = XN * XN.transpose();

    // 6x6 stiffness matrix (3 DOF per node: ux,uy,uz,rx,ry,rz)
    Eigen::MatrixXd KE = Eigen::MatrixXd::Zero(6, 6);

    // Axial block (translational DOFs)
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            KE(i, j) = DSCL * N(i, j);

    // Torsional block (rotational DOFs)
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            KE(3 + i, 3 + j) = DSCR * N(i, j);

    // Output as JSON
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << KE(i, j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}