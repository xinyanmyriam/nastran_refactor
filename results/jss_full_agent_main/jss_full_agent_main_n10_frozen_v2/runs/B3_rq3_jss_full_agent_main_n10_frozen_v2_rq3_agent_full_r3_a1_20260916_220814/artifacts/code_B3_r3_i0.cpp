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

    Vector3d ivec, jvec, kvec;
    ivec(0) = x2 - x1;
    ivec(1) = y2 - y1;
    ivec(2) = z2 - z1;
    double XSUBB = ivec.norm();
    ivec /= XSUBB;

    Vector3d rc;
    rc(0) = x3 - x1;
    rc(1) = y3 - y1;
    rc(2) = z3 - z1;
    double XSUBC = ivec.dot(rc);

    kvec = ivec.cross(rc);
    double YSUBC = kvec.norm();
    kvec /= YSUBC;

    jvec = kvec.cross(ivec);
    double TEMP = jvec.norm();
    if (TEMP != 0.0) jvec /= TEMP;

    double VOL = XSUBB * YSUBC * t / 2.0;
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB - 1.0;

    MatrixXd C(3, 6);
    C(0,0) = -REELMU;  C(0,1) = 0.0;
    C(1,0) = 0.0;      C(1,1) = FLAMDA*DELTA;
    C(2,0) = FLAMDA*DELTA; C(2,1) = -REELMU;
    C(0,2) = REELMU;   C(0,3) = 0.0;
    C(1,2) = 0.0;      C(1,3) = -FLAMDA*REELMU*XSUBC;
    C(2,2) = -FLAMDA*REELMU*XSUBC; C(2,3) = REELMU;
    C(0,4) = 0.0;      C(0,5) = 0.0;
    C(1,4) = 0.0;      C(1,5) = FLAMDA;
    C(2,4) = FLAMDA;   C(2,5) = 0.0;

    MatrixXd Emat(3, 2);
    Emat.col(0) = ivec;
    Emat.col(1) = jvec;

    double G11 = E_mod / (1.0 - nu*nu);
    double G12 = nu * E_mod / (1.0 - nu*nu);
    double G22 = G11;
    double G33 = E_mod / (2.0 * (1.0 + nu));

    Matrix3d G;
    G << G11, G12, 0.0,
         G12, G22, 0.0,
         0.0, 0.0, G33;

    // K_IJ = VOL * E * C_I^T * G * C_J * E^T  (3x3), take 2x2 upper-left
    MatrixXd K = MatrixXd::Zero(6, 6);
    for (int I = 0; I < 3; ++I) {
        MatrixXd CI = C.block(0, 2*I, 3, 2);
        MatrixXd common = Emat * CI.transpose() * G * VOL;  // 3x3
        for (int J = 0; J < 3; ++J) {
            MatrixXd CJ = C.block(0, 2*J, 3, 2);
            MatrixXd block = common * CJ * Emat.transpose();  // 3x3
            K.block(2*I, 2*J, 2, 2) = block.topLeftCorner(2,2);
        }
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 6; ++i) {
        std::cout << "[";
        for (int j = 0; j < 6; ++j) {
            std::cout << K(i, j);
            if (j < 5) std::cout << ",";
        }
        std::cout << "]";
        if (i < 5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}