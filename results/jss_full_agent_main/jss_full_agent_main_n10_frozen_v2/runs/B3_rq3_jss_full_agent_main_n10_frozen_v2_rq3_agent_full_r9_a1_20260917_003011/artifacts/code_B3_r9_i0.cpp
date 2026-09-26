#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    double x3 = 1.0, y3 = 1.5, z3 = 0.0;
    double E_mod = 2.1e11;
    double nu = 0.3;
    double t = 0.01;

    double A = 0.5 * std::abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1));
    
    Matrix<double,3,6> B;
    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
    B << b1, 0, b2, 0, b3, 0,
         0, c1, 0, c2, 0, c3,
         c1, b1, c2, b2, c3, b3;
    B /= (2.0*A);
    
    Matrix3d D;
    double fac = E_mod / (1.0 - nu*nu);
    D << fac, fac*nu, 0,
         fac*nu, fac, 0,
         0, 0, fac*(1.0-nu)/2.0;
    
    Matrix<double,6,6> Kstd = t * A * B.transpose() * D * B;
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; i++) {
        std::cout << "[";
        for (int j = 0; j < 6; j++) {
            std::cout << Kstd(i,j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}