#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    const double PI = 3.14159265358979323846;
    const double DCR = 0.017453292;

    Vector3d GPA(1.0, 0.0, 0.0);
    Vector3d GPB(0.0, 0.0, 1.0);
    Vector3d SMALLV(0.0, 1.0, 0.0);
    double R = 1.0;
    double BETAR = 90.0;
    double E = 200e9;
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));

    double ro = 0.05;
    double t_wall = 0.005;
    double ri = ro - t_wall;
    double A = PI * (ro*ro - ri*ri);
    double I1 = PI/4.0 * (std::pow(ro,4) - std::pow(ri,4));
    double I2 = I1;
    double FJ = 2.0 * I1;

    double K1 = 1.0, K2 = 1.0;

    Vector3d VECJ = SMALLV.normalized();
    Vector3d SMALV0 = GPB - GPA;
    SMALV0.normalize();
    Vector3d VECK = SMALV0.cross(VECJ);
    VECK.normalize();
    Vector3d VECI = VECJ.cross(VECK);
    VECI.normalize();

    double T = BETAR * DCR;
    double RA  = R / (A * E);
    double RV1 = R / (2.0 * K1 * G * A);
    double RV2 = K1 / K2 * RV1;
    double RT  = R / (G * FJ * 2.0);
    double RB0 = R / (E * I2 * 2.0);
    double RB1 = R / (E * I1);
    double R2  = R * R;

    double ST  = std::sin(BETAR * DCR);
    double CT  = std::cos(BETAR * DCR);
    double S2T = std::sin(2.0 * BETAR * DCR);
    double C2T = std::cos(2.0 * BETAR * DCR);

    Matrix<double,6,6> F = Matrix<double,6,6>::Zero();

    F(0,0) += 0.25*RA*(2.0*T + S2T);
    F(1,1) += 0.25*RA*(2.0*T - S2T);
    F(0,1) += 0.50*RA*ST*ST;

    F(0,0) += 0.5*RV1*(2.0*T - S2T);
    F(1,1) += 0.5*RV1*(2.0*T + S2T);
    F(2,2) += 2.0*RV2*T;
    F(0,1) -= RV1*ST*ST;

    F(2,2) += 0.5*RT*R2*(6.0*T + S2T - 8.0*ST);
    F(3,3) += 0.5*RT*(2.0*T + S2T);
    F(4,4) += 0.5*RT*(2.0*T - S2T);
    F(2,3) += RT*R*(ST - T*CT);
    F(2,4) += RT*R*(2.0 - 2.0*CT - T*ST);
    F(3,4) += 0.5*RT*(1.0 - C2T);

    F(0,0) += 0.25*RB1*R2*(2.0*T*(2.0 + C2T) - 3.0*S2T);
    F(1,1) += 0.25*RB1*R2*(2.0*T*(2.0 - C2T) + 3.0*S2T - 8.0*ST);
    F(2,2) += 0.50*RB0*R2*(2.0*T - S2T);
    F(3,3) += 0.50*RB0*(2.0*T - S2T);
    F(4,4) += 0.50*RB0*(2.0*T + S2T);
    F(5,5) += RB1*T;
    F(0,1) += 0.25*RB1*R2*(1.0 + 3.0*C2T + 2.0*T*S2T - 4.0*CT);
    F(0,5) -= RB1*R*(ST - T*CT);
    F(1,5) += RB1*R*(T*ST + CT - 1.0);
    F(2,3) += RB0*R*(ST - T*CT);
    F(2,4) -= RB0*R*T*ST;
    F(3,4) -= 0.50*RB0*(1.0 - C2T);

    Matrix<double,6,6> DF = Matrix<double,6,6>::Zero();
    for (int I = 0; I < 6; I++)
        for (int K = I; K < 6; K++) { DF(K,I)=F(I,K); DF(I,K)=DF(K,I); }

    Matrix<double,6,6> DFinv = DF.inverse();

    Matrix<double,6,6> H = Matrix<double,6,6>::Zero();
    for (int K = 0; K < 6; K++) H(K,K) = -1.0;
    H(3,2) = -(R*(1.0 - CT));
    H(4,2) = (R*ST);
    H(5,0) = -H(3,2);
    H(5,1) = -H(4,2);

    Matrix<double,12,12> S = Matrix<double,12,12>::Zero();
    for (int K = 0; K < 6; K++)
        for (int I = K; I < 6; I++) S(K+6, I+6) = DFinv(K, I);
    for (int IR = 0; IR < 6; IR++)
        for (int IC = 0; IC < 6; IC++) {
            double sum = 0.0;
            for (int IN = 0; IN < 6; IN++) sum += H(IR,IN) * DFinv(IN,IC);
            S(IR, IC+6) = sum;
        }
    for (int IR = 0; IR < 6; IR++)
        for (int IC = IR; IC < 6; IC++) {
            double sum = 0.0;
            for (int IN = 0; IN < 6; IN++) sum += S(IR, IN+6) * H(IC, IN);
            S(IR, IC) = sum;
        }
    for (int I = 0; I < 12; I++)
        for (int K = I; K < 12; K++) S(K, I) = S(I, K);

    Matrix3d Amat;
    Amat.row(0) = VECI.transpose();
    Amat.row(1) = VECJ.transpose();
    Amat.row(2) = VECK.transpose();
    Matrix<double,12,12> Tmat = Matrix<double,12,12>::Zero();
    for (int b = 0; b < 4; b++) Tmat.block<3,3>(3*b, 3*b) = Amat;
    Matrix<double,12,12> Kg = Tmat.transpose() * S * Tmat;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        std::cout << "[";
        for (int j = 0; j < 12; j++) {
            std::cout << Kg(i,j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}