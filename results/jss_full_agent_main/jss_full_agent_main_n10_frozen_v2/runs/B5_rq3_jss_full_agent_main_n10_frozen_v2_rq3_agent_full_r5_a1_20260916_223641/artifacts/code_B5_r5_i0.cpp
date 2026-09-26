#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct Common {
    double DEGRA = M_PI / 180.0;
    double G11=0,G12=0,G13=0,G22=0,G23=0,G33=0;
    double G2X211=0,G2X212=0,G2X222=0;
    double GSUBE=0;
    double E=0, nu=0, t=0, EYE=0;
};

static Common C;

MatrixXd inverd(const MatrixXd& M, int& ising) {
    int n = M.rows();
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = 2; return MatrixXd::Zero(n,n); }
    ising = 1;
    return lu.inverse();
}

void mat_routine() {
    double E = C.E, nu = C.nu;
    double f = E / (1.0 - nu*nu);
    C.G11 = f; C.G12 = f*nu; C.G13 = 0.0;
    C.G22 = f; C.G23 = 0.0; C.G33 = f*(1.0-nu)/2.0;
    C.G2X211 = 0; C.G2X212 = 0; C.G2X222 = 0;
}

struct KTRBSCResult {
    MatrixXd A;
    MatrixXd HINV;
    MatrixXd S;
};

KTRBSCResult ktrbsc(int iopt,
                    const Vector3d& p1, const Vector3d& p2, const Vector3d& p3,
                    double angle_deg) {
    KTRBSCResult res;

    Vector3d d2 = p2 - p1;
    Vector3d d1 = p3 - p1;
    double xsubb = d2.norm();
    Vector3d ivec = d2 / xsubb;
    Vector3d kvec = ivec.cross(d1);
    double ysubc = kvec.norm();
    kvec /= ysubc;
    Vector3d jvec = kvec.cross(ivec);
    jvec.normalize();
    double xsubc = ivec.dot(d1);

    mat_routine();
    Matrix3d G;
    G << C.G11, C.G12, C.G13,
         C.G12, C.G22, C.G23,
         C.G13, C.G23, C.G33;
    Matrix3d D = G * C.EYE;

    double area = xsubb * ysubc / 2.0;
    double xbar = (xsubb + xsubc) / 3.0;
    double ybar = ysubc / 3.0;
    double xcsq = xsubc*xsubc, ycsq = ysubc*ysubc, xbsq = xsubb*xsubb;
    double xcyc = xsubc*ysubc;
    double px2 = (xbsq + xsubb*xsubc + xcsq)/6.0;
    double py2 = ycsq/6.0;
    double pxy2 = ysubc*(xsubb + 2.0*xsubc)/12.0;
    double xbar3 = 3.0*xbar, ybar3 = 3.0*ybar, ybar2 = 2.0*ybar;

    double D1=D(0,0), D2=D(1,0), D3=D(2,0), D5=D(1,1), D6=D(2,1), D9=D(2,2);

    MatrixXd Kx = MatrixXd::Zero(6,6);
    auto sK = [&](int idx, double v){ int c=(idx-1)/6; int r=(idx-1)%6; Kx(r,c)=v; };
    auto gK = [&](int idx)->double{ int c=(idx-1)/6; int r=(idx-1)%6; return Kx(r,c); };
    sK(1, D1); sK(2, D3); sK(3, D2);
    sK(4, D1*xbar3); sK(5, D2*xbar + ybar2*D3); sK(6, D2*ybar3);
    sK(7, gK(2)); sK(8, D9); sK(9, D6);
    sK(10, D3*xbar3); sK(11, D6*xbar + ybar2*D9); sK(12, D6*ybar3);
    sK(13, gK(3)); sK(14, gK(9)); sK(15, D5);
    sK(16, D2*xbar3); sK(17, D5*xbar + ybar2*D6); sK(18, D5*ybar3);
    sK(19, gK(4)); sK(20, gK(10)); sK(21, gK(16));
    sK(22, D1*9.0*px2); sK(23, D2*3.0*px2 + 6.0*pxy2*D3); sK(24, D2*9.0*pxy2);
    sK(25, gK(5)); sK(26, gK(11)); sK(27, gK(17));
    sK(28, gK(23)); sK(29, D5*px2 + 4.0*pxy2*D6 + 4.0*py2*D9); sK(30, D5*3.0*pxy2 + 6.0*py2*D6);
    sK(31, gK(6)); sK(32, gK(12)); sK(33, gK(18));
    sK(34, gK(24)); sK(35, gK(30)); sK(36, D5*9.0*py2);
    double temp = 4.0*area;
    for (int i=1;i<=36;i++) sK(i, gK(i)*temp);

    MatrixXd Hbar = MatrixXd::Zero(6,6);
    auto sH = [&](int idx, double v){ int c=(idx-37)/6; int r=(idx-37)%6; Hbar(r,c)=v; };
    sH(37, xbsq); sH(40, xbsq*xsubb); sH(44, xsubb);
    sH(49, -2.0*xsubb); sH(52, -3.0*xbsq);
    sH(55, xcsq); sH(56, xcyc); sH(57, ycsq);
    sH(58, xcsq*xsubc); sH(59, ycsq*xsubc); sH(60, ycsq*ysubc);
    sH(62, xsubc); sH(63, ysubc*2.0);
    sH(65, xcyc*2.0); sH(66, ycsq*3.0);
    sH(67, -2.0*xsubc); sH(68, -ysubc);
    sH(70, -3.0*xcsq); sH(71, -ycsq);

    int ising;
    MatrixXd Hinv = inverd(Hbar, ising);

    MatrixXd KQHinv = Kx * Hinv;
    MatrixXd KII = Hinv.transpose() * KQHinv;

    MatrixXd S = MatrixXd::Zero(6,3);
    S(0,0)=1.0; S(0,2)=-xsubb;
    S(1,1)=1.0;
    S(2,2)=1.0;
    S(3,0)=1.0; S(3,1)=ysubc; S(3,2)=-xsubc;
    S(4,1)=1.0;
    S(5,2)=1.0;

    MatrixXd KIA = -(KII * S);
    MatrixXd KAA = S.transpose() * (KII * S);

    VectorXd A = VectorXd::Zero(145);
    auto sA = [&](int i, double v){ A(i)=v; };
    auto gA = [&](int i)->double{ return A(i); };

    for (int c=0;c<3;c++) for (int r=0;r<3;r++) sA(1+c*3+r, KAA(r,c));
    for (int c=0;c<3;c++) for (int r=0;r<6;r++) sA(46+c*6+r, KIA(r,c));
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) sA(109+c*6+r, KII(r,c));

    for (int i=28;i<=36;i++) sA(i, gA(i+18));
    sA(10,gA(46)); sA(11,gA(49)); sA(12,gA(52));
    sA(13,gA(47)); sA(14,gA(50)); sA(15,gA(53));
    sA(16,gA(48)); sA(17,gA(51)); sA(18,gA(54));
    sA(19,gA(55)); sA(20,gA(58)); sA(21,gA(61));
    sA(22,gA(56)); sA(23,gA(59)); sA(24,gA(62));
    sA(25,gA(57)); sA(26,gA(60)); sA(27,gA(63));
    sA(37,gA(109)); sA(38,gA(110)); sA(39,gA(111));
    sA(40,gA(115)); sA(41,gA(116)); sA(42,gA(117));
    sA(43,gA(121)); sA(44,gA(122)); sA(45,gA(123));
    sA(46,gA(112)); sA(47,gA(113)); sA(48,gA(114));
    sA(49,gA(118)); sA(50,gA(119)); sA(51,gA(120));
    sA(52,gA(124)); sA(53,gA(125)); sA(54,gA(126));
    sA(64,gA(127)); sA(65,gA(128)); sA(66,gA(129));
    sA(67,gA(133)); sA(68,gA(134)); sA(69,gA(135));
    sA(70,gA(139)); sA(71,gA(140)); sA(72,gA(141));
    sA(73,gA(130)); sA(74,gA(131)); sA(75,gA(132));
    sA(76,gA(136)); sA(77,gA(137)); sA(78,gA(138));
    sA(79,gA(142)); sA(80,gA(143)); sA(81,gA(144));

    MatrixXd KU = MatrixXd::Zero(9,9);
    for (int c=0;c<9;c++) for (int r=0;r<9;r++) KU(r,c) = gA(1+c*9+r);

    res.A = KU;
    res.HINV = Hinv;
    res.S = S;
    return res;
}

int main() {
    C.E = 200e9; C.nu = 0.3; C.t = 0.01;
    C.EYE = C.t*C.t*C.t/12.0;

    Vector3d pA(0,0,0), pB(1,0,0), pC(0,1,0);
    int M[9] = {1,2,4, 2,3,4, 3,1,4};

    Vector3d d2 = pB - pA, d1 = pC - pA;
    double xsubb = d2.norm();
    Vector3d ivec = d2/xsubb;
    Vector3d kvec = ivec.cross(d1);
    double ysubc = kvec.norm();
    kvec /= ysubc;
    Vector3d jvec = kvec.cross(ivec);
    jvec.normalize();
    double xsubc = ivec.dot(d1);

    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,1) = xsubb; R(0,2) = xsubc; R(1,2) = ysubc;
    R(0,3) = (xsubb + xsubc)/3.0; R(1,3) = ysubc/3.0;

    VectorXd KSUM = VectorXd::Zero(64);
    VectorXd G = VectorXd::Zero(37);
    auto sKSUM = [&](int i, double v){ KSUM(i)=v; };
    auto gKSUM = [&](int i)->double{ return KSUM(i); };
    auto sG = [&](int i, double v){ G(i)=v; };
    auto gG = [&](int i)->double{ return G(i); };

    for (int J=1; J<=3; J++) {
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0], SUBSCB = M[KM+1], SUBSCC = M[KM+2];

        Vector2d V, VV;
        for (int i=0;i<2;i++) {
            V(i) = R(i,SUBSCB-1) - R(i,SUBSCA-1);
            VV(i) = R(i,SUBSCC-1) - R(i,SUBSCA-1);
        }
        double XSUBB = V.norm();
        double U1 = V(0)/XSUBB, U2 = V(1)/XSUBB;
        double XSUBC = U1*VV(0) + U2*VV(1);
        double YSUBC = U1*VV(1) - U2*VV(0);

        Vector3d q1(R(0,SUBSCA-1), R(1,SUBSCA-1), 0);
        Vector3d q2(R(0,SUBSCB-1), R(1,SUBSCB-1), 0);
        Vector3d q3(R(0,SUBSCC-1), R(1,SUBSCC-1), 0);

        KTRBSCResult r = ktrbsc(2, q1, q2, q3, 0.0);

        Matrix3d Tmat;
        Tmat << 1.0, 0.0, 0.0,
                0.0, U1, -U2,
                0.0, U2, U1;

        VectorXd KUflat = VectorXd::Zero(82);
        for (int c=0;c<9;c++) for (int rr=0;rr<9;rr++) KUflat(1+c*9+rr) = r.A(rr,c);

        for (int I=1; I<=3; I++) {
            Matrix3d KUblk;
            int base = 27*I-8;
            for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KUblk(rr,c) = KUflat(base + c*3 + rr);
            Matrix3d TEMP9 = Tmat.transpose() * KUblk;
            Matrix3d PROD9 = TEMP9 * Tmat;
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] + 18;
            for (int K=1;K<=9;K++) {
                int NSUBC = NPOINT + K;
                sKSUM(NSUBC, gKSUM(NSUBC) + PROD9((K-1)%3, (K-1)/3));
            }
        }
        for (int K=1;K<=2;K++) {
            int NPOINT = KM + K;
            if (M[NPOINT-1] != 1) continue;
            Matrix3d KUblk;
            int base = 36*K-35;
            for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KUblk(rr,c) = KUflat(base + c*3 + rr);
            Matrix3d TEMP9 = Tmat.transpose() * KUblk;
            Matrix3d PROD9 = TEMP9 * Tmat;
            NPOINT = 9*1 - 9;
            for (int I=1;I<=9;I++) {
                int NSUBC = NPOINT + I;
                sKSUM(NSUBC, gKSUM(NSUBC) + PROD9((I-1)%3, (I-1)/3));
            }
            base = 18*K-8;
            for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KUblk(rr,c) = KUflat(base + c*3 + rr);
            TEMP9 = Tmat.transpose() * KUblk;
            PROD9 = TEMP9 * Tmat;
            NPOINT = KM + 3 - K;
            NPOINT = 9*M[NPOINT-1] - 9;
            for (int I=1;I<=9;I++) {
                int NSUBC = NPOINT + I;
                sKSUM(NSUBC, gKSUM(NSUBC) + PROD9((I-1)%3, (I-1)/3));
            }
        }

        double TEMP1 = XSUBB - XSUBC;
        double TEMP2 = YSUBC*YSUBC;
        double L1 = sqrt(XSUBC*XSUBC + TEMP2);
        double L2 = sqrt(TEMP1*TEMP1 + TEMP2);
        double S1 = XSUBC/L1, S2 = TEMP1/L2;
        double C1 = YSUBC/L1, C2 = YSUBC/L2;
        double X1 = XSUBC/2.0, Y1 = YSUBC/2.0;
        double X2 = (XSUBB + XSUBC)/2.0, Y2 = Y1;
        MatrixXd HQ(2,6);
        HQ(0,0) = -XSUBC*C1;
        HQ(0,1) = X1*S1 - Y1*C1;
        HQ(0,2) = 2.0*Y1*S1;
        HQ(0,3) = -3.0*X1*X1*C1;
        HQ(0,4) = Y1*(2.0*X1*S1 - Y1*C1);
        HQ(0,5) = 3.0*Y1*Y1*S1;
        HQ(1,0) = 2.0*X2*C2;
        HQ(1,1) = X2*S2 + Y2*C2;
        HQ(1,2) = 2.0*Y2*S2;
        HQ(1,3) = 3.0*X2*X2*C2;
        HQ(1,4) = Y2*(2.0*X2*S2 + Y2*C2);
        HQ(1,5) = 3.0*Y2*Y2*S2;

        MatrixXd PROD12 = HQ * r.HINV;
        MatrixXd HABC = -(PROD12 * r.S);
        HABC(0,0) += S1; HABC(0,1) += C1;
        HABC(1,0) += S2; HABC(1,1) -= C2;

        MatrixXd HPSI_A(2,3), HPSI_B(2,3), HPSI_C(2,3);
        HPSI_A(0,0)=HABC(0,0); HPSI_A(0,1)=HABC(0,1); HPSI_A(0,2)=HABC(0,2);
        HPSI_A(1,0)=HABC(1,0); HPSI_A(1,1)=HABC(1,1); HPSI_A(1,2)=HABC(1,2);
        HPSI_B(0,0)=PROD12(0,0); HPSI_B(0,1)=PROD12(0,1); HPSI_B(0,2)=PROD12(0,2);
        HPSI_B(1,0)=PROD12(1,0); HPSI_B(1,1)=PROD12(1,1); HPSI_B(1,2)=PROD12(1,2);
        HPSI_C(0,0)=PROD12(0,3); HPSI_C(0,1)=PROD12(0,4); HPSI_C(0,2)=PROD12(0,5);
        HPSI_C(1,0)=PROD12(1,3); HPSI_C(1,1)=PROD12(1,4); HPSI_C(1,2)=PROD12(1,5);

        MatrixXd HPSI[3] = {HPSI_A, HPSI_B, HPSI_C};

        for (int I=1;I<=3;I++) {
            MatrixXd TEMP9 = HPSI[I-1] * Tmat;
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] - 9;
            if (J==1) {
                for (int K=1;K<=6;K++) {
                    NPOINT++;
                    sG(NPOINT, gG(NPOINT) + TEMP9((K-1)%2, (K-1)/2));
                }
            } else if (J==2) {
                NPOINT += 3;
                for (int K=1;K<=6;K++) {
                    NPOINT++;
                    sG(NPOINT, gG(NPOINT) + TEMP9((K-1)%2, (K-1)/2));
                }
            } else {
                sG(NPOINT+7, gG(NPOINT+7) + TEMP9(0,0));
                sG(NPOINT+8, gG(NPOINT+8) + TEMP9(0,1));
                sG(NPOINT+9, gG(NPOINT+9) + TEMP9(0,2));
                sG(NPOINT+1, gG(NPOINT+1) + TEMP9(1,0));
                sG(NPOINT+2, gG(NPOINT+2) + TEMP9(1,1));
                sG(NPOINT+3, gG(NPOINT+3) + TEMP9(1,2));
            }
        }
    }

    MatrixXd E = MatrixXd::Zero(6,3);
    E(0,0)=kvec(0); E(0,1)=kvec(1); E(0,2)=kvec(2);
    E(4,0)=ivec(0); E(4,1)=ivec(1); E(4,2)=ivec(2);
    E(5,0)=jvec(0); E(5,1)=jvec(1); E(5,2)=jvec(2);
    MatrixXd TITE = E;

    Matrix3d G4;
    for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) G4(rr,c) = gG(28 + c*3 + rr);
    int ising;
    Matrix3d G4inv = inverd(G4, ising);

    Matrix3d Gpiv;
    for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) Gpiv(rr,c) = gG(1 + c*3 + rr);
    Matrix3d PROD9 = G4inv * Gpiv;

    Matrix3d KSUM55;
    for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KSUM55(rr,c) = gKSUM(55 + c*3 + rr);
    Matrix3d TEMP9 = PROD9.transpose() * KSUM55;

    MatrixXd KOUTs[3];
    for (int J=1;J<=3;J++) {
        Matrix3d KSUMJ4;
        int b = 9*J+19;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KSUMJ4(rr,c) = gKSUM(b + c*3 + rr);
        Matrix3d ARR9 = PROD9.transpose() * KSUMJ4.transpose();

        int NBEGIN = 9*J-9;
        for (int I=1;I<=9;I++) {
            int NPOINT = NBEGIN + I;
            sKSUM(NPOINT, gKSUM(NPOINT) - ARR9((I-1)%3, (I-1)/3));
        }

        Matrix3d GJ;
        b = 9*J-8;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) GJ(rr,c) = gG(b + c*3 + rr);
        ARR9 = G4inv * GJ;

        Matrix3d KSUMpiv4;
        b = 9*1+19;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KSUMpiv4(rr,c) = gKSUM(b + c*3 + rr);
        Matrix3d ARRAY9 = KSUMpiv4 * ARR9;

        for (int I=1;I<=9;I++) {
            int NPOINT = NBEGIN + I;
            sKSUM(NPOINT, gKSUM(NPOINT) - ARRAY9((I-1)%3, (I-1)/3));
        }

        ARRAY9 = TEMP9 * ARR9;
        for (int I=1;I<=9;I++) {
            int NPOINT = NBEGIN + I;
            sKSUM(NPOINT, gKSUM(NPOINT) + ARRAY9((I-1)%3, (I-1)/3));
        }

        MatrixXd TJTE = E;
        Matrix3d KSUMblk;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KSUMblk(rr,c) = gKSUM(NBEGIN+1 + c*3 + rr);
        MatrixXd TEMP18 = KSUMblk * TJTE.transpose();
        MatrixXd KOUT = TITE * TEMP18;
        KOUTs[J-1] = KOUT;
    }

    MatrixXd K9 = MatrixXd::Zero(9,9);
    for (int J=1;J<=3;J++) {
        MatrixXd KOUT = KOUTs[J-1];
        int pj = (J-1)*3;
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
            K9(i, j) += KOUT(i,j);
            K9(i, pj+j) += KOUT(i, 3+j);
            K9(pj+i, j) += KOUT(3+i, j);
            K9(pj+i, pj+j) += KOUT(3+i, 3+j);
        }
    }

    cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<9;i++) {
        cout << "[";
        for (int j=0;j<9;j++) {
            cout << scientific << setprecision(6) << K9(i,j);
            if (j<8) cout << ",";
        }
        cout << "]";
        if (i<8) cout << ",";
    }
    cout << "]}" << endl;
    return 0;
}