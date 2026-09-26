#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double x1=0.0, y1=0.0, z1=0.0;
    double x2=2.0, y2=0.0, z2=0.0;
    double x3=1.0, y3=1.5, z3=0.0;
    double E_mod = 2.1e11;
    double nu = 0.3;
    double t = 0.01;

    double e1 = x2 - x1;
    double e3 = y2 - y1;
    double e5 = z2 - z1;
    double XSUBB = std::sqrt(e1*e1 + e3*e3 + e5*e5);
    e1 /= XSUBB; e3 /= XSUBB; e5 /= XSUBB;

    double e2 = x3 - x1;
    double e4 = y3 - y1;
    double e6 = z3 - z1;

    double XSUBC = e1*e2 + e3*e4 + e5*e6;

    double e7 = e3*e6 - e5*e4;
    double e8 = e5*e2 - e1*e6;
    double e9 = e1*e4 - e3*e2;
    double YSUBC = std::sqrt(e7*e7 + e8*e8 + e9*e9);
    e7 /= YSUBC; e8 /= YSUBC; e9 /= YSUBC;

    e2 = e5*e8 - e3*e9;
    e4 = e1*e9 - e5*e7;
    e6 = e3*e7 - e1*e8;
    double TEMP = std::sqrt(e2*e2 + e4*e4 + e6*e6);
    e2 /= TEMP; e4 /= TEMP; e6 /= TEMP;

    MatrixXd Emat(3,2);
    Emat(0,0)=e1; Emat(1,0)=e3; Emat(2,0)=e5;
    Emat(0,1)=e2; Emat(1,1)=e4; Emat(2,1)=e6;

    double VOL = XSUBB*YSUBC*t/2.0;
    double REELMU = 1.0/XSUBB;
    double FLAMDA = 1.0/YSUBC;
    double DELTA = XSUBC/XSUBB - 1.0;

    MatrixXd C(3,6);
    C(0,0) = -REELMU;  C(0,1) = 0.0;
    C(1,0) = 0.0;      C(1,1) = FLAMDA*DELTA;
    C(2,0) = FLAMDA*DELTA; C(2,1) = -REELMU;

    C(0,2) = REELMU;   C(0,3) = 0.0;
    C(1,2) = 0.0;      C(1,3) = -FLAMDA*REELMU*XSUBC;
    C(2,2) = -FLAMDA*REELMU*XSUBC; C(2,3) = REELMU;

    C(0,4) = 0.0;      C(0,5) = 0.0;
    C(1,4) = 0.0;      C(1,5) = FLAMDA;
    C(2,4) = FLAMDA;   C(2,5) = 0.0;

    double G11 = E_mod/(1.0 - nu*nu);
    double G12 = nu*E_mod/(1.0 - nu*nu);
    double G33 = E_mod/(2.0*(1.0+nu));
    MatrixXd G(3,3);
    G << G11, G12, 0.0,
         G12, G11, 0.0,
         0.0, 0.0, G33;

    MatrixXd K = VOL * C.transpose() * G * C;

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<6;i++){
        std::cout << "[";
        for (int j=0;j<6;j++){
            std::cout << K(i,j);
            if (j<5) std::cout << ",";
        }
        std::cout << "]";
        if (i<5) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}