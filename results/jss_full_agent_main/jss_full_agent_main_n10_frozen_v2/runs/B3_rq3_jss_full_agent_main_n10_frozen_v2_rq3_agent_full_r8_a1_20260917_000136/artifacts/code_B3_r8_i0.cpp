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

    Vector3d ivec(x2-x1, y2-y1, z2-z1);
    double XSUBB = ivec.norm();
    ivec /= XSUBB;

    Vector3d rsubc(x3-x1, y3-y1, z3-z1);
    double XSUBC = ivec.dot(rsubc);

    Vector3d kvec = ivec.cross(rsubc);
    double YSUBC = kvec.norm();
    kvec /= YSUBC;

    Vector3d jvec = kvec.cross(ivec);
    jvec.normalize();

    double VOL = XSUBB * YSUBC * t / 2.0;
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB - 1.0;

    MatrixXd C(3, 6);
    C << -REELMU, 0.0, REELMU, 0.0, 0.0, 0.0,
          0.0, FLAMDA*DELTA, 0.0, -FLAMDA*REELMU*XSUBC, 0.0, FLAMDA,
          FLAMDA*DELTA, -REELMU, -FLAMDA*REELMU*XSUBC, REELMU, FLAMDA, 0.0;

    double G11 = E_mod / (1.0 - nu*nu);
    double G12 = nu * G11;
    double G33 = E_mod / (2.0 * (1.0 + nu));

    Matrix3d G;
    G << G11, G12, 0.0,
         G12, G11, 0.0,
         0.0, 0.0, G33;

    MatrixXd K = VOL * C.transpose() * G * C;

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