#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static const double DCR = 0.017453292;

int main() {
    double E = 200e9;
    double nu = 0.3;
    double G = E / (2.0 * (1.0 + nu));
    
    double ro = 0.05;
    double t_wall = 0.005;
    double ri = ro - t_wall;
    
    double A = M_PI * (ro*ro - ri*ri);
    double I1 = M_PI / 4.0 * (pow(ro,4) - pow(ri,4));
    double I2 = I1;
    double FJ = 2.0 * I1;
    
    Vector3d GPA(1.0, 0.0, 0.0);
    Vector3d GPB(0.0, 0.0, 1.0);
    Vector3d SMALLV(0.0, 1.0, 0.0);
    
    double R = 1.0;
    double BETAR = 90.0;
    
    double K1 = 1.0, K2 = 1.0, KX = 1.0, KY = 1.0, KZ = 1.0;
    
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
    
    double ST  = sin(T);
    double CT  = cos(T);
    double S2T = sin(2.0 * T);
    double C2T = cos(2.0 * T);
    
    MatrixXd F = MatrixXd::Zero(6,6);
    
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
    
    MatrixXd DF = MatrixXd::Zero(6,6);
    for (int i = 0; i < 6; i++)
        for (int k = i; k < 6; k++) {
            DF(k,i) = F(i,k);
            DF(i,k) = DF(k,i);
        }
    
    DF = DF.inverse();
    
    MatrixXd H = MatrixXd::Zero(6,6);
    for (int k = 0; k < 6; k++) H(k,k) = -1.0;
    H(3,2) = -(R*(1.0 - CT));
    H(4,2) = (R*ST);
    H(5,0) = -H(3,2);
    H(5,1) = -H(4,2);
    
    MatrixXd S = MatrixXd::Zero(12,12);
    
    for (int k = 0; k < 6; k++)
        for (int i = k; i < 6; i++)
            S(k+6, i+6) = DF(k, i);
    
    for (int ir = 0; ir < 6; ir++)
        for (int ic = 0; ic < 6; ic++) {
            double sum = 0.0;
            for (int in = 0; in < 6; in++)
                sum += H(ir,in) * DF(in,ic);
            S(ir, ic+6) = sum;
        }
    
    for (int ir = 0; ir < 6; ir++)
        for (int ic = ir; ic < 6; ic++) {
            double sum = 0.0;
            for (int in = 0; in < 6; in++)
                sum += S(ir, in+6) * H(ic, in);
            S(ir, ic) = sum;
        }
    
    for (int i = 0; i < 12; i++)
        for (int k = i; k < 12; k++)
            S(k,i) = S(i,k);
    
    Matrix3d AT;
    AT.row(0) = VECI.transpose();
    AT.row(1) = VECJ.transpose();
    AT.row(2) = VECK.transpose();
    
    Matrix3d Gmat = AT;
    
    MatrixXd W = MatrixXd::Zero(6,6);
    W.block<3,3>(0,0) = Gmat;
    W.block<3,3>(3,3) = Gmat;
    
    MatrixXd Wfull = MatrixXd::Zero(12,12);
    Wfull.block<6,6>(0,0) = W;
    Wfull.block<6,6>(6,6) = W;
    
    MatrixXd Kglobal = Wfull.transpose() * S * Wfull;
    
    cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        cout << "[";
        for (int j = 0; j < 12; j++) {
            cout << scientific << setprecision(6) << Kglobal(i,j);
            if (j < 11) cout << ",";
        }
        cout << "]";
        if (i < 11) cout << ",";
    }
    cout << "]}" << endl;
    
    return 0;
}