#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

int main() {
    double E = 200e9, G = 76.923e9, A = 0.01;
    double Iy = 8.333e-6, Iz = 8.333e-6, J = 1.667e-5;
    double L = 2.0;
    double K1 = 0.0, K2 = 0.0, I12 = 0.0;

    double LSQ = L*L;
    double LCUBE = LSQ*L;

    double EI1 = E*Iy;
    double EI2 = E*Iz;

    double R1, R2;
    if (K1 == 0.0 || I12 != 0.0) R1 = 12.0*EI1/LCUBE;
    else {
        double GAK1 = G*A*K1;
        R1 = (12.0*EI1*GAK1)/(GAK1*LCUBE + 12.0*L*EI1);
    }
    if (K2 == 0.0 || I12 != 0.0) R2 = 12.0*EI2/LCUBE;
    else {
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
    double GJL = G*J/L;

    Eigen::MatrixXd KE = Eigen::MatrixXd::Zero(12,12);
    auto set = [&](int idx, double v){
        int col = (idx-1)/12;
        int row = (idx-1)%12;
        KE(row,col) = v;
    };
    set(1, AEL); set(7,-AEL);
    set(14,R1); set(18,LR1); set(20,-R1); set(24,LR1);
    set(27,R2); set(29,-LR2); set(33,-R2); set(35,-LR2);
    set(40,GJL); set(46,-GJL);
    set(51,-LR2); set(53,SK2); set(57,LR2); set(59,SK4);
    set(62,LR1); set(66,SK1); set(68,-LR1); set(72,SK3);
    set(73,-AEL); set(79,AEL);
    set(86,-R1); set(90,-LR1); set(92,R1); set(96,-LR1);
    set(99,-R2); set(101,LR2); set(105,R2); set(107,LR2);
    set(112,-GJL); set(118,GJL);
    set(123,-LR2); set(125,SK4); set(129,LR2); set(131,SK2);
    set(134,LR1); set(138,SK3); set(140,-LR1); set(144,SK1);

    std::cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<12;i++){
        std::cout << "[";
        for (int j=0;j<12;j++){
            std::cout << std::scientific << std::setprecision(6) << KE(i,j);
            if (j<11) std::cout << ",";
        }
        std::cout << "]";
        if (i<11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}