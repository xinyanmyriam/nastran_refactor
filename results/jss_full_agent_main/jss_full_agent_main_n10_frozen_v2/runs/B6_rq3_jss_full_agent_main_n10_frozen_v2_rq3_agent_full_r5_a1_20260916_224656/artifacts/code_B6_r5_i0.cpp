#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct MatOut {
    double G11, G12, G13, G22, G23, G33;
    double G2X211, G2X212, G2X222;
};

static MatOut mat_iso(double E, double nu, double t) {
    MatOut m;
    double I = t*t*t/12.0;
    double D0 = E*I/(1.0 - nu*nu);
    m.G11 = D0;
    m.G12 = D0*nu;
    m.G13 = 0.0;
    m.G22 = D0;
    m.G23 = 0.0;
    m.G33 = D0*(1.0-nu)/2.0;
    m.G2X211 = 0.0; m.G2X212 = 0.0; m.G2X222 = 0.0;
    return m;
}

static MatrixXd inverd(const MatrixXd& M, int& ising) {
    int n = M.rows();
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = 2; return MatrixXd::Zero(n,n); }
    ising = 0;
    return lu.inverse();
}

static void ktrbsc(int IOPT,
                   double x1,double y1,double z1,
                   double x2,double y2,double z2,
                   double x3,double y3,double z3,
                   double angle, double eye,
                   double E, double nu, double t,
                   MatrixXd& Aout, int& nogo) {
    nogo = 0;
    VectorXd Evec = VectorXd::Zero(18);
    double XSUBB, XSUBC, YSUBC;

    Evec(10) = x2 - x1;
    Evec(13) = y2 - y1;
    Evec(16) = z2 - z1;
    XSUBB = sqrt(Evec(10)*Evec(10) + Evec(13)*Evec(13) + Evec(16)*Evec(16));
    if (XSUBB <= 1.0e-6) { nogo = 1; return; }
    Evec(10) /= XSUBB; Evec(13) /= XSUBB; Evec(16) /= XSUBB;

    Evec(1) = x3 - x1;
    Evec(4) = y3 - y1;
    Evec(7) = z3 - z1;
    XSUBC = Evec(10)*Evec(1) + Evec(13)*Evec(4) + Evec(16)*Evec(7);

    Evec(0) = Evec(13)*Evec(7) - Evec(4)*Evec(16);
    Evec(3) = Evec(1)*Evec(16) - Evec(10)*Evec(7);
    Evec(6) = Evec(10)*Evec(4) - Evec(1)*Evec(13);
    YSUBC = sqrt(Evec(0)*Evec(0) + Evec(3)*Evec(3) + Evec(6)*Evec(6));
    if (YSUBC <= 1.0e-6) { nogo = 1; return; }
    Evec(0) /= YSUBC; Evec(3) /= YSUBC; Evec(6) /= YSUBC;

    Evec(11) = Evec(3)*Evec(16) - Evec(13)*Evec(6);
    Evec(14) = Evec(10)*Evec(6) - Evec(0)*Evec(16);
    Evec(17) = Evec(0)*Evec(13) - Evec(10)*Evec(3);
    double TEMP = sqrt(Evec(11)*Evec(11) + Evec(14)*Evec(14) + Evec(17)*Evec(17));
    Evec(11) /= TEMP; Evec(14) /= TEMP; Evec(17) /= TEMP;
    Evec(1) = 0; Evec(2) = 0; Evec(4) = 0; Evec(5) = 0; Evec(7) = 0; Evec(8) = 0;
    Evec(9) = 0; Evec(12) = 0; Evec(15) = 0;

    MatOut mo = mat_iso(E, nu, t);
    Matrix3d G;
    G << mo.G11, mo.G12, mo.G13,
         mo.G12, mo.G22, mo.G23,
         mo.G13, mo.G23, mo.G33;
    Matrix3d D = G * eye;

    double AREA = XSUBB*YSUBC/2.0;
    double XBAR = (XSUBB+XSUBC)/3.0;
    double YBAR = YSUBC/3.0;
    double XCSQ = XSUBC*XSUBC;
    double YCSQ = YSUBC*YSUBC;
    double XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ+XSUBB*XSUBC+XCSQ)/6.0;
    double PY2 = YCSQ/6.0;
    double PXY2 = YSUBC*(XSUBB+2.0*XSUBC)/12.0;
    double XBAR3 = 3.0*XBAR;
    double YBAR3 = 3.0*YBAR;
    double YBAR2 = 2.0*YBAR;

    VectorXd A = VectorXd::Zero(225);
    auto a = [&](int i)->double& { return A(i-1); };

    a(1) = D(0,0);
    a(2) = D(0,2);
    a(3) = D(0,1);
    a(4) = D(0,0)*XBAR3;
    a(5) = D(0,1)*XBAR + YBAR2*D(0,2);
    a(6) = D(0,1)*YBAR3;
    a(7) = a(2);
    a(8) = D(2,2);
    a(9) = D(1,2);
    a(10) = D(0,2)*XBAR3;
    a(11) = D(1,2)*XBAR + YBAR2*D(2,2);
    a(12) = D(1,2)*YBAR3;
    a(13) = a(3);
    a(14) = a(9);
    a(15) = D(1,1);
    a(16) = D(0,1)*XBAR3;
    a(17) = D(1,1)*XBAR + YBAR2*D(1,2);
    a(18) = D(1,1)*YBAR3;
    a(19) = a(4);
    a(20) = a(10);
    a(21) = a(16);
    a(22) = D(0,0)*9.0*PX2;
    a(23) = D(0,1)*3.0*PX2 + 6.0*PXY2*D(0,2);
    a(24) = D(0,1)*9.0*PXY2;
    a(25) = a(5);
    a(26) = a(11);
    a(27) = a(17);
    a(28) = a(23);
    a(29) = D(1,1)*PX2 + 4.0*PXY2*D(1,2) + 4.0*PY2*D(2,2);
    a(30) = D(1,1)*3.0*PXY2 + 6.0*PY2*D(1,2);
    a(31) = a(6);
    a(32) = a(12);
    a(33) = a(18);
    a(34) = a(24);
    a(35) = a(30);
    a(36) = D(1,1)*9.0*PY2;
    TEMP = 4.0*AREA;
    for (int i = 1; i <= 36; i++) a(i) *= TEMP;

    for (int i = 37; i <= 72; i++) a(i) = 0.0;
    a(37) = XBSQ;
    a(40) = XBSQ*XSUBB;
    a(44) = XSUBB;
    a(49) = -2.0*XSUBB;
    a(52) = -3.0*XBSQ;
    a(55) = XCSQ;
    a(56) = XCYC;
    a(57) = YCSQ;
    a(58) = XCSQ*XSUBC;
    a(59) = YCSQ*XSUBC;
    a(60) = YCSQ*YSUBC;
    a(62) = XSUBC;
    a(63) = YSUBC*2.0;
    a(65) = XCYC*2.0;
    a(66) = YCSQ*3.0;
    a(67) = -2.0*XSUBC;
    a(68) = -YSUBC;
    a(70) = -3.0*XCSQ;
    a(71) = -YCSQ;

    MatrixXd H = MatrixXd::Zero(6,6);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            H(i,j) = a(37 + i*6 + j);
    int ising;
    MatrixXd Hinv = inverd(H, ising);
    if (ising == 2) { nogo = 1; return; }
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            a(37 + i*6 + j) = Hinv(i,j);

    MatrixXd KX = MatrixXd::Zero(6,6);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            KX(i,j) = a(1 + i*6 + j);
    MatrixXd KQ = KX * Hinv;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            a(73 + i*6 + j) = KQ(i,j);

    MatrixXd KII = Hinv.transpose() * KQ;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            a(109 + i*6 + j) = KII(i,j);

    a(82) = 1.0; a(83) = 0.0; a(84) = -XSUBB;
    a(85) = 0.0; a(86) = 1.0; a(87) = 0.0;
    a(88) = 0.0; a(89) = 0.0; a(90) = 1.0;
    a(91) = 1.0; a(92) = YSUBC; a(93) = -XSUBC;
    a(94) = 0.0; a(95) = 1.0; a(96) = 0.0;
    a(97) = 0.0; a(98) = 0.0; a(99) = 1.0;

    MatrixXd S = MatrixXd::Zero(6,3);
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            S(i,j) = a(82 + i*3 + j);
    MatrixXd KIA = KII * S;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            a(46 + i*3 + j) = KIA(i,j);

    MatrixXd KAA = S.transpose() * KIA;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            a(1 + i*3 + j) = KAA(i,j);

    for (int i = 46; i <= 63; i++) a(i) = -a(i);

    for (int i = 28; i <= 36; i++) a(i) = a(i+18);
    a(10) = a(46); a(11) = a(49); a(12) = a(52);
    a(13) = a(47); a(14) = a(50); a(15) = a(53);
    a(16) = a(48); a(17) = a(51); a(18) = a(54);
    a(19) = a(55); a(20) = a(58); a(21) = a(61);
    a(22) = a(56); a(23) = a(59); a(24) = a(62);
    a(25) = a(57); a(26) = a(60); a(27) = a(63);
    a(37) = a(109); a(38) = a(110); a(39) = a(111);
    a(40) = a(115); a(41) = a(116); a(42) = a(117);
    a(43) = a(121); a(44) = a(122); a(45) = a(123);
    a(46) = a(112); a(47) = a(113); a(48) = a(114);
    a(49) = a(118); a(50) = a(119); a(51) = a(120);
    a(52) = a(124); a(53) = a(125); a(54) = a(126);
    a(64) = a(127); a(65) = a(128); a(66) = a(129);
    a(67) = a(133); a(68) = a(134); a(69) = a(135);
    a(70) = a(139); a(71) = a(140); a(72) = a(141);
    a(73) = a(130); a(74) = a(131); a(75) = a(132);
    a(76) = a(136); a(77) = a(137); a(78) = a(138);
    a(79) = a(142); a(80) = a(143); a(81) = a(144);

    Aout = A.head(81);
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double eye = t*t*t/12.0;

    double VQ[4][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0}};

    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};

    int NPIVOT = 1;
    int JNOT;
    if (NPIVOT - 2 <= 0) JNOT = NPIVOT + 2;
    else JNOT = NPIVOT - 2;

    double R[2][5];
    for (int i=0;i<2;i++) for(int j=0;j<5;j++) R[i][j]=0.0;

    double D1[3], D2[3], A1v[3];
    for (int i=0;i<3;i++) {
        D1[i] = VQ[2][i] - VQ[0][i];
        D2[i] = VQ[3][i] - VQ[1][i];
        A1v[i] = VQ[1][i] - VQ[0][i];
    }

    double KV[3], IV[3], JV[3];
    KV[0] = D1[1]*D2[2] - D2[1]*D1[2];
    KV[1] = D1[2]*D2[0] - D2[2]*D1[0];
    KV[2] = D1[0]*D2[1] - D2[0]*D1[1];
    double TEMP = sqrt(KV[0]*KV[0]+KV[1]*KV[1]+KV[2]*KV[2]);
    for (int i=0;i<3;i++) KV[i] /= TEMP;

    TEMP = (A1v[0]*KV[0]+A1v[1]*KV[1]+A1v[2]*KV[2])/2.0;
    for (int i=0;i<3;i++) IV[i] = A1v[i] - TEMP*KV[i];
    TEMP = sqrt(IV[0]*IV[0]+IV[1]*IV[1]+IV[2]*IV[2]);
    for (int i=0;i<3;i++) IV[i] /= TEMP;

    JV[0] = KV[1]*IV[2] - IV[1]*KV[2];
    JV[1] = KV[2]*IV[0] - IV[2]*KV[0];
    JV[2] = KV[0]*IV[1] - IV[0]*KV[1];
    TEMP = sqrt(JV[0]*JV[0]+JV[1]*JV[1]+JV[2]*JV[2]);
    for (int i=0;i<3;i++) JV[i] /= TEMP;

    R[0][3] = D1[0]*IV[0]+D1[1]*IV[1]+D1[2]*IV[2];
    R[0][2] = A1v[0]*IV[0]+A1v[1]*IV[1]+A1v[2]*IV[2];
    R[1][3] = D1[0]*JV[0]+D1[1]*JV[1]+D1[2]*JV[2];
    R[0][4] = D2[0]*IV[0]+D2[1]*IV[1]+D2[2]*IV[2] + R[0][2];
    R[1][4] = D2[0]*JV[0]+D2[1]*JV[1]+D2[2]*JV[2];

    VectorXd KSUM = VectorXd::Zero(36);

    for (int J = 1; J <= 4; J++) {
        if (J == JNOT) continue;
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0];
        int SUBSCB = M[KM+1];
        int SUBSCC = M[KM+2];

        double V[2], VV[2];
        for (int i=0;i<2;i++) {
            V[i]  = R[i][SUBSCB] - R[i][SUBSCA];
            VV[i] = R[i][SUBSCC] - R[i][SUBSCA];
        }
        double XSUBB = sqrt(V[0]*V[0]+V[1]*V[1]);
        double U1 = V[0]/XSUBB;
        double U2 = V[1]/XSUBB;

        double x1 = R[0][SUBSCA], y1 = R[1][SUBSCA];
        double x2 = R[0][SUBSCB], y2 = R[1][SUBSCB];
        double x3 = R[0][SUBSCC], y3 = R[1][SUBSCC];

        MatrixXd Aout;
        int nogo;
        ktrbsc(1, x1,y1,0.0, x2,y2,0.0, x3,y3,0.0,
               0.0, eye, E, nu, t, Aout, nogo);
        if (nogo) return 1;

        Matrix3d T;
        T << 1.0, 0.0, 0.0,
             0.0, U1,  U2,
             0.0, -U2, U1;

        int NBEGIN = 0;
        for (int i = 1; i <= 3; i++) {
            int NPOINT = KM + i;
            if (M[NPOINT-1] == NPIVOT) { NBEGIN = 27*i - 27; break; }
        }

        for (int i = 1; i <= 3; i++) {
            int NPOINT = NBEGIN + 9*i - 8;
            Matrix3d Ai;
            for (int r=0;r<3;r++) for(int c=0;c<3;c++) Ai(r,c) = Aout(NPOINT-1 + r*3 + c);
            Matrix3d TEMP9 = T.transpose() * Ai;
            Matrix3d PROD9 = TEMP9 * T;
            int NP = KM + i;
            NP = 9*M[NP-1] - 9;
            for (int k = 1; k <= 9; k++) {
                NP = NP + 1;
                KSUM(NP-1) += PROD9((k-1)/3, (k-1)%3)/2.0;
            }
        }
    }

    MatrixXd Emat = MatrixXd::Zero(6,3);
    Emat(0,0) = KV[0]; Emat(1,0) = KV[1]; Emat(2,0) = KV[2];
    Emat(3,0) = IV[0]; Emat(4,0) = IV[1]; Emat(5,0) = IV[2];
    Emat(3,1) = JV[0]; Emat(4,1) = JV[1]; Emat(5,1) = JV[2];

    MatrixXd TITE = Emat;
    MatrixXd KOUT = MatrixXd::Zero(12,12);

    for (int J = 1; J <= 4; J++) {
        MatrixXd TJTE = Emat;
        Matrix3d KSUMJ;
        for (int r=0;r<3;r++) for(int c=0;c<3;c++) KSUMJ(r,c) = KSUM(9*J-8-1 + r*3 + c);
        MatrixXd TEMP18 = KSUMJ * TJTE.transpose();
        MatrixXd KOUTb = TITE * TEMP18;
        int base = 3*(J-1);
        for (int r=0;r<6;r++) for(int c=0;c<6;c++) {
            int rr = (r<3) ? (3*(NPIVOT-1)+(r%3)) : (base+(r%3));
            int cc = (c<3) ? (3*(NPIVOT-1)+(c%3)) : (base+(c%3));
            KOUT(rr,cc) += KOUTb(r,c);
        }
    }

    cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<12;i++) {
        cout << "[";
        for (int j=0;j<12;j++) {
            cout << scientific << setprecision(6) << KOUT(i,j);
            if (j<11) cout << ",";
        }
        cout << "]";
        if (i<11) cout << ",";
    }
    cout << "]}" << endl;

    return 0;
}