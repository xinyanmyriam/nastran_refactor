#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double E = 200e9;
    double G = 76.923e9;
    double A = 0.01;
    double I1 = 8.333e-6;
    double I2 = 8.333e-6;
    double FJ = 1.667e-5;
    double K1 = 0.0;
    double K2 = 0.0;
    double I12 = 0.0;

    Vector3d gpa(0, 0, 0);
    Vector3d gpb(2, 0, 0);
    Vector3d smallv(0, 0, 1);

    Vector3d veci = gpa - gpb;
    veci = -veci;
    double FL = veci.norm();
    veci /= FL;

    Vector3d smalv0 = smallv;
    double fl2 = smalv0.norm();
    smalv0 /= fl2;

    Vector3d veck = veci.cross(smalv0);
    double fll = veck.norm();
    veck /= fll;

    Vector3d vecj = veck.cross(veci);
    fll = vecj.norm();
    vecj /= fll;

    double L = FL;
    double LSQ = L * L;
    double LCUBE = LSQ * L;

    double EI1 = E * I1;
    double EI2 = E * I2;

    double R1 = 12.0 * EI1 / LCUBE;
    double R2 = 12.0 * EI2 / LCUBE;

    double SK1 = 0.25 * R1 * LSQ + EI1 / L;
    double SK2 = 0.25 * R2 * LSQ + EI2 / L;
    double SK3 = 0.25 * R1 * LSQ - EI1 / L;
    double SK4 = 0.25 * R2 * LSQ - EI2 / L;

    double AEL = A * E / L;
    double LR1 = L * R1 / 2.0;
    double LR2 = L * R2 / 2.0;
    double GJL = G * FJ / L;

    MatrixXd KE = MatrixXd::Zero(12, 12);
    auto setKE = [&](int idx, double val) {
        int i = (idx - 1) % 12;
        int j = (idx - 1) / 12;
        KE(i, j) = val;
    };

    setKE(1, AEL); setKE(7, -AEL);
    setKE(14, R1); setKE(18, LR1); setKE(20, -R1); setKE(24, LR1);
    setKE(27, R2); setKE(29, -LR2); setKE(33, -R2); setKE(35, -LR2);
    setKE(40, GJL); setKE(46, -GJL);
    setKE(51, -LR2); setKE(53, SK2); setKE(57, LR2); setKE(59, SK4);
    setKE(62, LR1); setKE(66, SK1); setKE(68, -LR1); setKE(72, SK3);
    setKE(73, -AEL); setKE(79, AEL);
    setKE(86, -R1); setKE(90, -LR1); setKE(92, R1); setKE(96, -LR1);
    setKE(99, -R2); setKE(101, LR2); setKE(105, R2); setKE(107, LR2);
    setKE(112, -GJL); setKE(118, GJL);
    setKE(123, -LR2); setKE(125, SK4); setKE(129, LR2); setKE(131, SK2);
    setKE(134, LR1); setKE(138, SK3); setKE(140, -LR1); setKE(144, SK1);

    Matrix3d AT;
    AT(0,0)=veci(0); AT(1,0)=veci(1); AT(2,0)=veci(2);
    AT(0,1)=vecj(0); AT(1,1)=vecj(1); AT(2,1)=vecj(2);
    AT(0,2)=veck(0); AT(1,2)=veck(1); AT(2,2)=veck(2);

    Matrix3d Gmat = AT;
    Matrix<double,6,6> W = Matrix<double,6,6>::Zero();
    W.block<3,3>(0,0) = Gmat;
    W.block<3,3>(3,3) = Gmat;

    Matrix<double,12,12> T = Matrix<double,12,12>::Zero();
    T.block<6,6>(0,0) = W;
    T.block<6,6>(6,6) = W;

    Matrix<double,12,12> K_global = T.transpose() * KE * T;

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        std::cout << "[";
        for (int j = 0; j < 12; j++) {
            std::cout << std::scientific << std::setprecision(6) << K_global(i, j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}