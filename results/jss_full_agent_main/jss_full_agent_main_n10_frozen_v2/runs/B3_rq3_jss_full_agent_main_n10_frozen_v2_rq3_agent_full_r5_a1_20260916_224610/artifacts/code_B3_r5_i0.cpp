#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    // Test case: triangular membrane (CST), 3 nodes, 2 DOF/node
    double x1 = 0.0, y1 = 0.0, z1 = 0.0;
    double x2 = 2.0, y2 = 0.0, z2 = 0.0;
    double x3 = 1.0, y3 = 1.5, z3 = 0.0;
    double E_mod = 2.1e11;
    double nu = 0.3;
    double t = 0.01;

    // Build local element coordinate system (I, J, K vectors)
    Vector3d Ivec, Jvec, Kvec;
    Ivec << x2 - x1, y2 - y1, z2 - z1;
    double XSUBB = Ivec.norm();
    Ivec /= XSUBB;

    Vector3d RC;
    RC << x3 - x1, y3 - y1, z3 - z1;
    double XSUBC = Ivec.dot(RC);

    Kvec = Ivec.cross(RC);
    double YSUBC = Kvec.norm();
    Kvec /= YSUBC;

    Jvec = Kvec.cross(Ivec);
    double TEMP = Jvec.norm();
    if (TEMP != 0.0) Jvec /= TEMP;

    double area = XSUBB * YSUBC / 2.0;

    // Local coordinates: node1 (0,0), node2 (XSUBB,0), node3 (XSUBC,YSUBC)
    double x1l = 0.0, y1l = 0.0;
    double x2l = XSUBB, y2l = 0.0;
    double x3l = XSUBC, y3l = YSUBC;

    double b1 = y2l - y3l;
    double b2 = y3l - y1l;
    double b3 = y1l - y2l;
    double c1 = x3l - x2l;
    double c2 = x1l - x3l;
    double c3 = x2l - x1l;

    // Strain-displacement matrix B (3x6)
    Matrix<double, 3, 6> B;
    B << b1, 0, b2, 0, b3, 0,
         0, c1, 0, c2, 0, c3,
         c1, b1, c2, b2, c3, b3;
    B /= (2.0 * area);

    // Isotropic plane-stress material matrix D (3x3)
    Matrix3d D;
    double fac = E_mod / (1.0 - nu*nu);
    D << fac, fac*nu, 0.0,
         fac*nu, fac, 0.0,
         0.0, 0.0, fac*(1.0-nu)/2.0;

    // Stiffness matrix K = t * A * B^T * D * B
    Matrix<double, 6, 6> K = t * area * B.transpose() * D * B;

    // Output as JSON
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << K(i,j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}