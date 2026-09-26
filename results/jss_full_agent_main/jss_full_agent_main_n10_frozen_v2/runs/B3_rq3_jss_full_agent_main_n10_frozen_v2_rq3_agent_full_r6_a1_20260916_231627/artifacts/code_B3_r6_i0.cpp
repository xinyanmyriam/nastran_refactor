#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double x1 = 0.0, y1 = 0.0;
    double x2 = 2.0, y2 = 0.0;
    double x3 = 1.0, y3 = 1.5;
    double E_mod = 2.1e11;
    double nu = 0.3;
    double t = 0.01;

    double A = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));

    MatrixXd B(3,6);
    B << (y2-y3), 0, (y3-y1), 0, (y1-y2), 0,
         0, (x3-x2), 0, (x1-x3), 0, (x2-x1),
         (x3-x2), (y2-y3), (x1-x3), (y3-y1), (x2-x1), (y1-y2);
    B /= (2.0*A);

    double c = E_mod / (1.0 - nu*nu);
    Matrix3d D;
    D << c, c*nu, 0,
         c*nu, c, 0,
         0, 0, c*(1.0-nu)/2.0;

    MatrixXd K = t * A * B.transpose() * D * B;

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; i++) {
        std::cout << "[";
        for (int j = 0; j < 6; j++) {
            std::cout << K(i,j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}