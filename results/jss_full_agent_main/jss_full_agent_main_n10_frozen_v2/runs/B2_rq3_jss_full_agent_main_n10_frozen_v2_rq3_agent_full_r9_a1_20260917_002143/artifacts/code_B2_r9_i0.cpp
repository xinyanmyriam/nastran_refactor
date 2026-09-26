#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace Eigen;

int main() {
    double E  = 200e9;
    double G  = 76.923e9;
    double A  = 0.01;
    double I1 = 8.333e-6;
    double I2 = 8.333e-6;
    double FJ = 1.667e-5;
    double K1 = 0.0;
    double K2 = 0.0;
    double I12 = 0.0;

    Vector3d PA(0.0, 0.0, 0.0);
    Vector3d PB(2.0, 0.0, 0.0);
    Vector3d SMALLV(0.0, 0.0, 1.0);

    Vector3d VECI = PA - PB;
    VECI = -VECI;
    double FL = VECI.norm();
    VECI /= FL;

    Vector3d SMALV0 = SMALLV;
    double fl2 = SMALV0.norm();
    if (fl2 > 0.0) SMALV0 /= fl2;

    Vector3d VECK = VECI.cross(SMALV0);
    double FLL = VECK.norm();
    VECK /= FLL;

    Vector3d VECJ = VECK.cross(VECI);
    FLL = VECJ.norm();
    VECJ /= FLL;

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

    std::vector<double> KE(145, 0.0);

    KE[1]  =  AEL;
    KE[7]  = -AEL;
    KE[14] =  R1;
    KE[18] =  LR1;
    KE[20] = -R1;
    KE[24] =  LR1;
    KE[27] =  R2;
    KE[29] = -LR2;
    KE[33] = -R2;
    KE[35] = -LR2;
    KE[40] =  GJL;
    KE[46] = -GJL;
    KE[51] = -LR2;
    KE[53] =  SK2;
    KE[57] =  LR2;
    KE[59] =  SK4;
    KE[62] =  LR1;
    KE[66] =  SK1;
    KE[68] = -LR1;
    KE[72] =  SK3;
    KE[73] = -AEL;
    KE[79] =  AEL;
    KE[86] = -R1;
    KE[90] = -LR1;
    KE[92] =  R1;
    KE[96] = -LR1;
    KE[99] = -R2;
    KE[101] =  LR2;
    KE[105] =  R2;
    KE[107] =  LR2;
    KE[112] = -GJL;
    KE[118] =  GJL;
    KE[123] = -LR2;
    KE[125] =  SK4;
    KE[129] =  LR2;
    KE[131] =  SK2;
    KE[134] =  LR1;
    KE[138] =  SK3;
    KE[140] = -LR1;
    KE[144] =  SK1;

    MatrixXd K(12, 12);
    for (int col = 0; col < 12; ++col)
        for (int row = 0; row < 12; ++row)
            K(row, col) = KE[col * 12 + row + 1];

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; ++i) {
        std::cout << "[";
        for (int j = 0; j < 12; ++j) {
            std::cout << K(i, j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}