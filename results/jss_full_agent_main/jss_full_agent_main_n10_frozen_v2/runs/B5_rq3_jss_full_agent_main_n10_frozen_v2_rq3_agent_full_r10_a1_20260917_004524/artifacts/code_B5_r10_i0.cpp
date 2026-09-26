#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

double G11, G12, G13, G22, G23, G33;

void MAT_isotropic(double E, double nu) {
    double factor = E / (1.0 - nu*nu);
    G11 = factor; G12 = factor*nu; G13 = 0.0;
    G22 = factor; G23 = 0.0; G33 = factor*(1.0-nu)/2.0;
}

MatrixXd inverd(const MatrixXd& A, int n, bool& ising) {
    FullPivLU<MatrixXd> lu(A);
    if (!lu.isInvertible()) { ising = true; return MatrixXd::Zero(n,n); }
    ising = false;
    return lu.inverse();
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double EYE = t*t*t/12.0;
    MAT_isotropic(E, nu);
    double X1=0,Y1=0,Z1=0, X2=1,Y2=0,Z2=0, X3=0,Y3=1,Z3=0;
    MatrixXd E6 = MatrixXd::Zero(6,3);
    double e11 = X2-X1, e14 = Y2-Y1, e17 = Z2-Z1;
    double XSUBB = sqrt(e11*e11+e14*e14+e17*e17);
    e11/=XSUBB; e14/=XSUBB; e17/=XSUBB;
    double e2 = X3-X1, e5 = Y3-Y1, e8 = Z3-Z1;
    double XSUBC = e11*e2+e14*e5+e17*e8;
    double e1 = e14*e8-e5*e17, e4 = e2*e17-e11*e8, e7 = e11*e5-e2*e14;
    double YSUBC = sqrt(e1*e1+e4*e4+e7*e7);
    e1/=YSUBC; e4/=YSUBC; e7/=YSUBC;
    double e12 = e4*e17-e14*e7, e15 = e11*e7-e1*e17, e18 = e1*e14-e11*e4;
    double TEMP = sqrt(e12*e12+e15*e15+e18*e18);
    e12/=TEMP; e15/=TEMP; e18/=TEMP;
    E6(0,0)=e1; E6(3,0)=e4; E6(0,1)=e7;
    E6(4,1)=e11; E6(1,2)=e14; E6(4,2)=e17;
    E6(5,1)=e12; E6(2,2)=e15; E6(5,2)=e18;
    Matrix3d G;
    G << G11,G12,G13, G12,G22,G23, G13,G23,G33;
    Matrix3d D = G * EYE;
    double d1=D(0,0),d2=D(1,0),d3=D(2,0),d5=D(1,1),d6=D(2,1),d9=D(2,2);
    double AREA = XSUBB*YSUBC/2.0;
    double XBAR = (XSUBB+XSUBC)/3.0, YBAR = YSUBC/3.0;
    double XCSQ=XSUBC*XSUBC, YCSQ=YSUBC*YSUBC, XBSQ=XSUBB*XSUBB;
    double XCYC=XSUBC*YSUBC;
    double PX2=(XBSQ+XSUBB*XSUBC+XCSQ)/6.0, PY2=YCSQ/6.0;
    double PXY2=YSUBC*(XSUBB+2.0*XSUBC)/12.0;
    double XBAR3=3.0*XBAR, YBAR3=3.0*YBAR, YBAR2=2.0*YBAR;
    MatrixXd KX = MatrixXd::Zero(6,6);
    KX(0,0)=d1; KX(0,1)=d3; KX(0,2)=d2; KX(0,3)=d1*XBAR3; KX(0,4)=d2*XBAR+YBAR2*d3; KX(0,5)=d2*YBAR3;
    KX(1,0)=KX(0,1); KX(1,1)=d9; KX(1,2)=d6; KX(1,3)=d3*XBAR3; KX(1,4)=d6*XBAR+YBAR2*d9; KX(1,5)=d6*YBAR3;
    KX(2,0)=KX(0,2); KX(2,1)=KX(1,2); KX(2,2)=d5; KX(2,3)=d2*XBAR3; KX(2,4)=d5*XBAR+YBAR2*d6; KX(2,5)=d5*YBAR3;
    KX(3,0)=KX(0,3); KX(3,1)=KX(1,3); KX(3,2)=KX(2,3); KX(3,3)=d1*9.0*PX2; KX(3,4)=d2*3.0*PX2+6.0*PXY2*d3; KX(3,5)=d2*9.0*PXY2;
    KX(4,0)=KX(0,4); KX(4,1)=KX(1,4); KX(4,2)=KX(2,4); KX(4,3)=KX(3,4); KX(4,4)=d5*PX2+4.0*PXY2*d6+4.0*PY2*d9; KX(4,5)=d5*3.0*PXY2+6.0*PY2*d6;
    KX(5,0)=KX(0,5); KX(5,1)=KX(1,5); KX(5,2)=KX(2,5); KX(5,3)=KX(3,5); KX(5,4)=KX(4,5); KX(5,5)=d5*9.0*PY2;
    KX *= 4.0*AREA;
    MatrixXd HBAR = MatrixXd::Zero(6,6);
    HBAR(0,0)=XBSQ; HBAR(3,0)=XBSQ*XSUBB; HBAR(1,1)=XSUBB;
    HBAR(0,2)=-2.0*XSUBB; HBAR(3,2)=-3.0*XBSQ;
    HBAR(0,3)=XCSQ; HBAR(1,3)=XCYC; HBAR(2,3)=YCSQ;
    HBAR(3,3)=XCSQ*XSUBC; HBAR(4,3)=YCSQ*XSUBC; HBAR(5,3)=YCSQ*YSUBC;
    HBAR(1,4)=XSUBC; HBAR(2,4)=YSUBC*2.0;
    HBAR(1,5)=XCYC*2.0; HBAR(2,5)=YCSQ*3.0;
    HBAR(0,4)=-2.0*XSUBC; HBAR(0,5)=-YSUBC;
    HBAR(3,4)=-3.0*XCSQ; HBAR(3,5)=-YCSQ;
    bool ising;
    MatrixXd Hinv = inverd(HBAR, 6, ising);
    MatrixXd KQ = KX * Hinv;
    MatrixXd KII = Hinv.transpose() * KQ;
    MatrixXd S = MatrixXd::Zero(6,3);
    S(0,0)=1.0; S(2,0)=-XSUBB; S(4,0)=1.0;
    S(2,1)=1.0; S(3,1)=1.0; S(4,1)=YSUBC; S(5,1)=-XSUBC;
    S(1,2)=1.0; S(5,2)=1.0;
    MatrixXd KIA_before = KII * S;
    MatrixXd KAA = S.transpose() * KIA_before;
    MatrixXd KIA = -KIA_before;
    MatrixXd KU = MatrixXd::Zero(9,9);
    KU.block<3,3>(0,0) = KAA;
    KU.block<3,3>(3,0) = KIA.block<3,3>(0,0);
    KU.block<3,3>(6,0) = KIA.block<3,3>(3,0);
    KU.block<3,3>(0,3) = KIA.block<3,3>(0,0).transpose();
    KU.block<3,3>(0,6) = KIA.block<3,3>(3,0).transpose();
    KU.block<3,3>(3,3) = KII.block<3,3>(0,0);
    KU.block<3,3>(3,6) = KII.block<3,3>(0,3);
    KU.block<3,3>(6,3) = KII.block<3,3>(3,0);
    KU.block<3,3>(6,6) = KII.block<3,3>(3,3);
    cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 9; i++) {
        cout << "[";
        for (int j = 0; j < 9; j++) {
            cout << scientific << setprecision(6) << KU(i,j);
            if (j < 8) cout << ",";
        }
        cout << "]";
        if (i < 8) cout << ",";
    }
    cout << "]}" << endl;
    return 0;
}