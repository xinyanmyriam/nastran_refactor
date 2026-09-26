#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    // Triangular membrane (CST) test case
    // Nodes: (0,0), (2,0), (1,1.5) m
    double x1=0.0, y1=0.0, x2=2.0, y2=0.0, x3=1.0, y3=1.5;
    double E_mod = 2.1e11;  // Pa
    double nu = 0.3;
    double t = 0.01;        // m

    // Triangle area
    double A = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

    // Strain-displacement (B) matrix coefficients
    double b1 = y2-y3, b2 = y3-y1, b3 = y1-y2;
    double c1 = x3-x2, c2 = x1-x3, c3 = x2-x1;

    MatrixXd B(3,6);
    B << b1,0,b2,0,b3,0,
         0,c1,0,c2,0,c3,
         c1,b1,c2,b2,c3,b3;
    B /= (2.0*A);

    // Isotropic plane-stress constitutive (D) matrix
    double f = E_mod/(1.0-nu*nu);
    Matrix3d D;
    D << f, nu*f, 0.0,
         nu*f, f, 0.0,
         0.0, 0.0, E_mod/(2.0*(1.0+nu));

    // Element stiffness matrix: K = t * A * B^T * D * B
    MatrixXd K = t * A * B.transpose() * D * B;

    // Print 6x6 stiffness matrix as JSON
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int r = 0; r < 6; r++) {
        std::cout << "[";
        for (int c = 0; c < 6; c++) {
            std::cout << K(r,c);
            if (c < 5) std::cout << ",";
        }
        std::cout << "]";
        if (r < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}