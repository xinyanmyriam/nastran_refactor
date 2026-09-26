#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace Eigen;
using namespace std;

struct Arr1 {
    vector<double> v;
    Arr1(int n) : v(n+1, 0.0) {}
    double& operator()(int i) { return v[i]; }
    double operator()(int i) const { return v[i]; }
};

MatrixXd gmmatd(const MatrixXd& A, int transA, const MatrixXd& B, int transB) {
    MatrixXd At = transA ? A.transpose() : A;
    MatrixXd Bt = transB ? B.transpose() : B;
    return At * Bt;
}

MatrixXd inverd(const MatrixXd& M, int& ising) {
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = 2; return MatrixXd::Zero(M.rows(), M.cols()); }
    ising = 1;
    return lu.inverse();
}

struct MatOut { double G11,G12,G13,G22,G23,G33; };

bool ktrbsc(Arr1& A,
            double X1,double Y1,double Z1,
            double X2,double Y2,double Z2,
            double X3,double Y3,double Z3,
            double EYE, double ANGLE, double T2,
            const MatOut& mat,
            double& XSUBB_out, double& XSUBC_out, double& YSUBC_out) {
    Vector3d Ivec, Jvec, Kvec;
    Vector3d d2(X2-X1, Y2-Y1, Z2-Z1);
    double XSUBB = d2.norm();
    if (XSUBB <= 1.0e-6) return false;
    Ivec = d2 / XSUBB;
    Vector3d d1(X3-X1, Y3-Y1, Z3-Z1);
    double XSUBC = Ivec.dot(d1);
    Kvec = Ivec.cross(d1);
    double YSUBC = Kvec.norm();
    if (YSUBC <= 1.0e-6) return false;
    Kvec = Kvec / YSUBC;
    Jvec = Kvec.cross(Ivec);
    Jvec.normalize();
    XSUBB_out = XSUBB; XSUBC_out = XSUBC; YSUBC_out = YSUBC;

    Matrix3d G;
    G << mat.G11, mat.G12, mat.G13,
         mat.G12, mat.G22, mat.G23,
         mat.G13, mat.G23, mat.G33;
    Matrix3d D = G * EYE;

    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    double XCSQ = XSUBC*XSUBC;
    double YCSQ = YSUBC*YSUBC;
    double XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ + XSUBB*XSUBC + XCSQ)/6.0;
    double PY2 = YCSQ/6.0;
    double PXY2 = YSUBC*(XSUBB + 2.0*XSUBC)/12.0;
    double XBAR3 = 3.0*XBAR;
    double YBAR3 = 3.0*YBAR;
    double YBAR2 = 2.0*YBAR;

    double Dv[10];
    Dv[1]=D(0,0); Dv[2]=D(1,0); Dv[3]=D(2,0);
    Dv[4]=D(0,1); Dv[5]=D(1,1); Dv[6]=D(2,1);
    Dv[7]=D(0,2); Dv[8]=D(1,2); Dv[9]=D(2,2);

    MatrixXd KX = MatrixXd::Zero(6,6);
    auto sA = [&](int i, int j, double val){ KX(i-1,j-1)=val; };
    auto gA = [&](int i, int j){ return KX(i-1,j-1); };
    sA(1,1, Dv[1]);
    sA(2,1, Dv[3]);
    sA(3,1, Dv[2]);
    sA(4,1, Dv[1]*XBAR3);
    sA(5,1, Dv[2]*XBAR + YBAR2*Dv[3]);
    sA(6,1, Dv[2]*YBAR3);
    sA(1,2, gA(2,1));
    sA(2,2, Dv[9]);
    sA(3,2, Dv[6]);
    sA(4,2, Dv[3]*XBAR3);
    sA(5,2, Dv[6]*XBAR + YBAR2*Dv[9]);
    sA(6,2, Dv[6]*YBAR3);
    sA(1,3, gA(3,1));
    sA(2,3, gA(3,2));
    sA(3,3, Dv[5]);
    sA(4,3, Dv[2]*XBAR3);
    sA(5,3, Dv[5]*XBAR + YBAR2*Dv[6]);
    sA(6,3, Dv[5]*YBAR3);
    sA(1,4, gA(4,1));
    sA(2,4, gA(4,2));
    sA(3,4, gA(4,3));
    sA(4,4, Dv[1]*9.0*PX2);
    sA(5,4, Dv[2]*3.0*PX2 + 6.0*PXY2*Dv[3]);
    sA(6,4, Dv[2]*9.0*PXY2);
    sA(1,5, gA(5,1));
    sA(2,5, gA(5,2));
    sA(3,5, gA(5,3));
    sA(4,5, gA(5,4));
    sA(5,5, Dv[5]*PX2 + 4.0*PXY2*Dv[6] + 4.0*PY2*Dv[9]);
    sA(6,5, Dv[5]*3.0*PXY2 + 6.0*PY2*Dv[6]);
    sA(1,6, gA(6,1));
    sA(2,6, gA(6,2));
    sA(3,6, gA(6,3));
    sA(4,6, gA(6,4));
    sA(5,6, gA(6,5));
    sA(6,6, Dv[5]*9.0*PY2);
    double TEMP = 4.0*AREA;
    KX = KX * TEMP;

    MatrixXd HBAR = MatrixXd::Zero(6,6);
    auto sH = [&](int i, int j, double val){ HBAR(i-1,j-1)=val; };
    sH(1,1, XBSQ);
    sH(4,1, XBSQ*XSUBB);
    sH(2,2, XSUBB);
    sH(1,3, -2.0*XSUBB);
    sH(4,3, -3.0*XBSQ);
    sH(1,4, XCSQ);
    sH(2,4, XCYC);
    sH(3,4, YCSQ);
    sH(4,4, XCSQ*XSUBC);
    sH(5,4, YCSQ*XSUBC);
    sH(6,4, YCSQ*YSUBC);
    sH(2,5, XSUBC);
    sH(3,5, YSUBC*2.0);
    sH(5,5, XCYC*2.0);
    sH(6,5, YCSQ*3.0);
    sH(2,6, -2.0*XSUBC);
    sH(3,6, -YSUBC);
    sH(5,6, -3.0*XCSQ);
    sH(6,6, -YCSQ);

    int ising;
    MatrixXd Hinv = inverd(HBAR, ising);
    if (ising == 2) return false;

    MatrixXd KQ = KX * Hinv;
    MatrixXd KII = Hinv.transpose() * KQ;

    MatrixXd S = MatrixXd::Zero(6,3);
    S(0,0)=1.0; S(0,1)=0.0; S(0,2)=-XSUBB;
    S(1,0)=0.0; S(1,1)=1.0; S(1,2)=0.0;
    S(2,0)=0.0; S(2,1)=0.0; S(2,2)=1.0;
    S(3,0)=1.0; S(3,1)=YSUBC; S(3,2)=-XSUBC;
    S(4,0)=0.0; S(4,1)=1.0; S(4,2)=0.0;
    S(5,0)=0.0; S(5,1)=0.0; S(5,2)=1.0;

    MatrixXd KIA_pre = KII * S;
    MatrixXd KAA = S.transpose() * KIA_pre;
    MatrixXd KIA = -KIA_pre;

    for (int i=0;i<3;i++) for (int j=0;j<3;j++) A(1 + j*3 + i) = KAA(i,j);
    A(10)=KIA(0,0); A(11)=KIA(3,0); A(12)=KIA(0,1);
    A(13)=KIA(1,0); A(14)=KIA(4,0); A(15)=KIA(1,1);
    A(16)=KIA(2,0); A(17)=KIA(5,0); A(18)=KIA(2,1);
    A(19)=KIA(3,1); A(20)=KIA(0,2); A(21)=KIA(3,2);
    A(22)=KIA(4,1); A(23)=KIA(1,2); A(24)=KIA(4,2);
    A(25)=KIA(5,1); A(26)=KIA(2,2); A(27)=KIA(5,2);
    A(28)=KIA(0,0); A(29)=KIA(1,0); A(30)=KIA(2,0);
    A(31)=KIA(3,0); A(32)=KIA(4,0); A(33)=KIA(5,0);
    A(34)=KIA(0,1); A(35)=KIA(1,1); A(36)=KIA(2,1);
    A(37)=KII(0,0); A(38)=KII(1,0); A(39)=KII(2,0);
    A(40)=KII(3,0); A(41)=KII(4,0); A(42)=KII(5,0);
    A(43)=KII(0,1); A(44)=KII(1,1); A(45)=KII(2,1);
    A(46)=KII(3,1); A(47)=KII(4,1); A(48)=KII(5,1);
    A(49)=KII(0,2); A(50)=KII(1,2); A(51)=KII(2,2);
    A(52)=KII(3,2); A(53)=KII(4,2); A(54)=KII(5,2);
    for (int i=0;i<6;i++) for (int j=0;j<3;j++) A(55 + j*6 + i) = S(i,j);
    A(64)=KII(0,3); A(65)=KII(1,3); A(66)=KII(2,3);
    A(67)=KII(3,3); A(68)=KII(4,3); A(69)=KII(5,3);
    A(70)=KII(0,4); A(71)=KII(1,4); A(72)=KII(2,4);
    A(73)=KII(3,4); A(74)=KII(4,4); A(75)=KII(5,4);
    A(76)=KII(0,5); A(77)=KII(1,5); A(78)=KII(2,5);
    A(79)=KII(3,5); A(80)=KII(4,5); A(81)=KII(5,5);
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) A(145 + j*6 + i) = Hinv(i,j);
    return true;
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double EYE = t*t*t/12.0;
    MatOut mat;
    mat.G11 = E/(1.0-nu*nu);
    mat.G12 = nu*E/(1.0-nu*nu);
    mat.G13 = 0.0;
    mat.G22 = E/(1.0-nu*nu);
    mat.G23 = 0.0;
    mat.G33 = E/(2.0*(1.0+nu));

    double X1=0,Y1=0,Z1=0;
    double X2=1,Y2=0,Z2=0;
    double X3=0,Y3=1,Z3=0;

    Arr1 A(225);
    double XSUBB, XSUBC, YSUBC;
    bool ok = ktrbsc(A, X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3, EYE, 0.0, 0.0, mat, XSUBB, XSUBC, YSUBC);
    if (!ok) { cout << "{\"stiffness_matrix\":[]}" << endl; return 1; }

    MatrixXd K = MatrixXd::Zero(9,9);
    for (int I=0; I<3; I++) {
        for (int J=0; J<3; J++) {
            int base = 9*((I)*3+(J)) + 1;
            for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
                K(3*I+i, 3*J+j) = A(base + j*3 + i);
            }
        }
    }

    cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<9;i++) {
        cout << "[";
        for (int j=0;j<9;j++) {
            cout << scientific << setprecision(6) << K(i,j);
            if (j<8) cout << ",";
        }
        cout << "]";
        if (i<8) cout << ",";
    }
    cout << "]}" << endl;
    return 0;
}