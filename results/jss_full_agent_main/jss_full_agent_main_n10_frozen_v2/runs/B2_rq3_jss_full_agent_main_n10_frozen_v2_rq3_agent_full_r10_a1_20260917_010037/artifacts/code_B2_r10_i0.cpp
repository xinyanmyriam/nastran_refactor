#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace Eigen;

int main() {
    double E = 200e9;
    double G = 76.923e9;
    double A = 0.01;
    double Iy = 8.333e-6;
    double Iz = 8.333e-6;
    double J = 1.667e-5;
    double L = 2.0;

    double I1 = Iy;
    double I2 = Iz;
    double FJ = J;
    double K1 = 0.0, K2 = 0.0, I12 = 0.0;

    double LSQ = L*L;
    double LCUBE = LSQ*L;

    double EI1 = E*I1;
    double EI2 = E*I2;

    double R1, R2;
    if (K1 == 0.0 || I12 != 0.0) {
        R1 = 12.0*EI1/LCUBE;
    } else {
        double GAK1 = G*A*K1;
        R1 = (12.0*EI1*GAK1)/(GAK1*LCUBE + 12.0*L*EI1);
    }
    if (K2 == 0.0 || I12 != 0.0) {
        R2 = 12.0*EI2/LCUBE;
    } else {
        double GAK2 = G*A*K2;
        R2 = (12.0*EI2*GAK2)/(GAK2*LCUBE + 12.0*L*EI2);
    }

    double SK1 = 0.25*R1*LSQ + EI1/L;
    double SK2 = 0.25*R2*LSQ + EI2/L;
    double SK3 = 0.25*R1*LSQ - EI1/L;
    double SK4 = 0.25*R2*LSQ - EI2/L;

    double AEL = A*E/L;
    double LR1 = L*R1/2.0;
    double LR2 = L*R2/2.0;
    double GJL = G*FJ/L;

    MatrixXd KE = MatrixXd::Zero(12,12);

    auto setKE = [&](int idx, double val) {
        int i = (idx-1) % 12;
        int j = (idx-1) / 12;
        KE(i,j) = val;
    };

    setKE(  1,  AEL);
    setKE(  7, -AEL);
    setKE( 14,  R1);
    setKE( 18,  LR1);
    setKE( 20, -R1);
    setKE( 24,  LR1);
    setKE( 27,  R2);
    setKE( 29, -LR2);
    setKE( 33, -R2);
    setKE( 35, -LR2);
    setKE( 40,  GJL);
    setKE( 46, -GJL);
    setKE( 51, -LR2);
    setKE( 53,  SK2);
    setKE( 57,  LR2);
    setKE( 59,  SK4);
    setKE( 62,  LR1);
    setKE( 66,  SK1);
    setKE( 68, -LR1);
    setKE( 72,  SK3);
    setKE( 73, -AEL);
    setKE( 79,  AEL);
    setKE( 86, -R1);
    setKE( 90, -LR1);
    setKE( 92,  R1);
    setKE( 96, -LR1);
    setKE( 99, -R2);
    setKE(101,  LR2);
    setKE(105,  R2);
    setKE(107,  LR2);
    setKE(112, -GJL);
    setKE(118,  GJL);
    setKE(123, -LR2);
    setKE(125,  SK4);
    setKE(129,  LR2);
    setKE(131,  SK2);
    setKE(134,  LR1);
    setKE(138,  SK3);
    setKE(140, -LR1);
    setKE(144,  SK1);

    if (I12 != 0.0) {
        double BETA = 12.0*E*I12/LCUBE;
        double LB = L*BETA/2.0;
        double L2B3 = LSQ*BETA/3.0;
        double L2B6 = LSQ*BETA/6.0;
        setKE( 15,  BETA);
        setKE( 17, -LB);
        setKE( 21, -BETA);
        setKE( 23, -LB);
        setKE( 26,  BETA);
        setKE( 30,  LB);
        setKE( 32, -BETA);
        setKE( 36,  LB);
        setKE( 50, -LB);
        setKE( 54, -L2B3);
        setKE( 56,  LB);
        setKE( 60, -L2B6);
        setKE( 63,  LB);
        setKE( 65, -L2B3);
        setKE( 69, -LB);
        setKE( 71, -L2B6);
        setKE( 87, -BETA);
        setKE( 89,  LB);
        setKE( 93,  BETA);
        setKE( 95,  LB);
        setKE( 98, -BETA);
        setKE(102, -LB);
        setKE(104,  BETA);
        setKE(108, -LB);
        setKE(122, -LB);
        setKE(126, -L2B6);
        setKE(128,  LB);
        setKE(132, -L2B3);
        setKE(135,  LB);
        setKE(137, -L2B6);
        setKE(141, -LB);
        setKE(143, -L2B3);
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        std::cout << "[";
        for (int j = 0; j < 12; j++) {
            std::cout << std::scientific << std::setprecision(6) << KE(i,j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}