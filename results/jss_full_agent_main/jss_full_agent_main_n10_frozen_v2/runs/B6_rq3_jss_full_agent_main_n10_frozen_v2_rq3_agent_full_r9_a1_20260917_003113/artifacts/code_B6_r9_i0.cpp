#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

static const double PI = 3.14159265358979323846;

struct MatOut {
    double G11, G12, G13, G22, G23, G33;
    double G2X211, G2X212, G2X222;
};
static MatOut matout;

static void MAT(int inflag, double E, double nu, double t) {
    if (inflag == 2) {
        double D0 = E * t * t * t / (12.0 * (1.0 - nu * nu));
        matout.G11 = D0;
        matout.G12 = D0 * nu;
        matout.G13 = 0.0;
        matout.G22 = D0;
        matout.G23 = 0.0;
        matout.G33 = D0 * (1.0 - nu) / 2.0;
    } else if (inflag == 3) {
        matout.G2X211 = 0.0;
        matout.G2X212 = 0.0;
        matout.G2X222 = 0.0;
    }
}

static int INVERD(MatrixXd& H) {
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) return 2;
    H = lu.inverse();
    return 0;
}

struct KTRBSCResult {
    MatrixXd KU;
    bool nogo;
};

static KTRBSCResult KTRBSC(double XSUBB, double XSUBC, double YSUBC,
                           double E, double nu, double t, double eye,
                           double T2) {
    KTRBSCResult res;
    res.nogo = false;
    res.KU = MatrixXd::Zero(9, 9);

    MAT(2, E, nu, t);
    Matrix3d G;
    G << matout.G11, matout.G12, matout.G13,
         matout.G12, matout.G22, matout.G23,
         matout.G13, matout.G23, matout.G33;

    Matrix3d D = G * eye;

    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;

    double XCSQ = XSUBC * XSUBC;
    double YCSQ = YSUBC * YSUBC;
    double XBSQ = XSUBB * XSUBB;
    double XCYC = XSUBC * YSUBC;
    double PX2 = (XBSQ + XSUBB * XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC * (XSUBB + 2.0 * XSUBC) / 12.0;
    double XBAR3 = 3.0 * XBAR;
    double YBAR3 = 3.0 * YBAR;
    double YBAR2 = 2.0 * YBAR;

    auto Dv = [&](int i) { return D(i - 1); };

    MatrixXd KX = MatrixXd::Zero(6, 6);
    KX(0,0) = Dv(1);
    KX(0,1) = Dv(3);
    KX(0,2) = Dv(2);
    KX(0,3) = Dv(1)*XBAR3;
    KX(0,4) = Dv(2)*XBAR + YBAR2*Dv(3);
    KX(0,5) = Dv(2)*YBAR3;
    KX(1,0) = Dv(2);
    KX(1,1) = Dv(9);
    KX(1,2) = Dv(6);
    KX(1,3) = Dv(3)*XBAR3;
    KX(1,4) = Dv(6)*XBAR + YBAR2*Dv(9);
    KX(1,5) = Dv(6)*YBAR3;
    KX(2,0) = Dv(3);
    KX(2,1) = Dv(6);
    KX(2,2) = Dv(5);
    KX(2,3) = Dv(2)*XBAR3;
    KX(2,4) = Dv(5)*XBAR + YBAR2*Dv(6);
    KX(2,5) = Dv(5)*YBAR3;
    KX(3,0) = KX(0,3);
    KX(3,1) = KX(1,3);
    KX(3,2) = KX(2,3);
    KX(3,3) = Dv(1)*9.0*PX2;
    KX(3,4) = Dv(2)*3.0*PX2 + 6.0*PXY2*Dv(3);
    KX(3,5) = Dv(2)*9.0*PXY2;
    KX(4,0) = KX(0,4);
    KX(4,1) = KX(1,4);
    KX(4,2) = KX(2,4);
    KX(4,3) = KX(3,4);
    KX(4,4) = Dv(5)*PX2 + 4.0*PXY2*Dv(6) + 4.0*PY2*Dv(9);
    KX(4,5) = Dv(5)*3.0*PXY2 + 6.0*PY2*Dv(6);
    KX(5,0) = KX(0,5);
    KX(5,1) = KX(1,5);
    KX(5,2) = KX(2,5);
    KX(5,3) = KX(3,5);
    KX(5,4) = KX(4,5);
    KX(5,5) = Dv(5)*9.0*PY2;
    KX *= 4.0 * AREA;

    MatrixXd HB = MatrixXd::Zero(6, 6);
    HB(0,0) = XBSQ;
    HB(0,3) = XBSQ*XSUBB;
    HB(1,1) = XSUBB;
    HB(2,0) = -2.0*XSUBB;
    HB(2,3) = -3.0*XBSQ;
    HB(3,0) = XCSQ;
    HB(3,1) = XCYC;
    HB(3,2) = YCSQ;
    HB(3,3) = XCSQ*XSUBC;
    HB(3,4) = YCSQ*XSUBC;
    HB(3,5) = YCSQ*YSUBC;
    HB(4,1) = XSUBC;
    HB(4,2) = YSUBC*2.0;
    HB(4,4) = XCYC*2.0;
    HB(4,5) = YCSQ*3.0;
    HB(5,0) = -2.0*XSUBC;
    HB(5,1) = -YSUBC;
    HB(5,3) = -3.0*XCSQ;
    HB(5,4) = -YCSQ;

    if (T2 != 0.0) {
        MAT(3, E, nu, t);
        if (!(matout.G2X211 == 0.0 && matout.G2X212 == 0.0 && matout.G2X222 == 0.0)) {
            Matrix2d G2X2;
            G2X2 << matout.G2X211*T2, matout.G2X212*T2,
                    matout.G2X212*T2, matout.G2X222*T2;
            double DETERM = G2X2(0,0)*G2X2(1,1) - G2X2(1,0)*G2X2(0,1);
            Matrix2d J2X2;
            J2X2 << G2X2(1,1)/DETERM, -G2X2(0,1)/DETERM,
                    -G2X2(1,0)/DETERM, G2X2(0,0)/DETERM;
            double TMP = 2.0*Dv(2) + 4.0*Dv(9);
            MatrixXd HYQ(2, 3);
            HYQ(0,0) = -6.0*(J2X2(0,0)*Dv(1) + J2X2(0,1)*Dv(3));
            HYQ(0,1) = -J2X2(0,0)*TMP - 6.0*J2X2(0,1)*Dv(6);
            HYQ(0,2) = -6.0*(J2X2(0,0)*Dv(6) + J2X2(0,1)*Dv(5));
            HYQ(1,0) = -6.0*(J2X2(1,0)*Dv(1) + J2X2(1,1)*Dv(3));
            HYQ(1,1) = -J2X2(1,0)*TMP - 6.0*J2X2(1,1)*Dv(6);
            HYQ(1,2) = -6.0*(J2X2(1,0)*Dv(6) + J2X2(1,1)*Dv(5));
            MatrixXd prod = G2X2 * HYQ;
            MatrixXd KYpart = HYQ.transpose() * prod;
            for (int i = 0; i < 3; i++) {
                KX(3+i, 3) += KYpart(i,0)*AREA;
                KX(3+i, 4) += KYpart(i,1)*AREA;
                KX(3+i, 5) += KYpart(i,2)*AREA;
            }
            for (int i = 0; i < 3; i++) {
                HB(1, 3+i) += XSUBB * HYQ(0, i);
                HB(3, 3+i) += XSUBC * HYQ(0, i) + YSUBC * HYQ(1, i);
            }
        }
    }

    MatrixXd HBinv = HB;
    int ising = INVERD(HBinv);
    if (ising == 2) { res.nogo = true; return res; }

    MatrixXd KQ = KX * HBinv;
    MatrixXd KII = HBinv.transpose() * KQ;

    MatrixXd S(6, 3);
    S << 1.0, 0.0, -XSUBB,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, YSUBC, -XSUBC,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

    MatrixXd KIA = KII * S;
    MatrixXd KAA = S.transpose() * KIA;
    KIA = -KIA;

    res.KU.block<3,3>(0,0) = KAA;
    res.KU.block<3,3>(3,0) = KIA.block<3,3>(0,0);
    res.KU.block<3,3>(6,0) = KIA.block<3,3>(3,0);
    res.KU.block<3,3>(0,3) = KIA.block<3,3>(0,0).transpose();
    res.KU.block<3,3>(0,6) = KIA.block<3,3>(3,0).transpose();
    res.KU.block<3,3>(3,3) = KII.block<3,3>(0,0);
    res.KU.block<3,3>(3,6) = KII.block<3,3>(0,3);
    res.KU.block<3,3>(6,3) = KII.block<3,3>(3,0);
    res.KU.block<3,3>(6,6) = KII.block<3,3>(3,3);

    return res;
}

int main() {
    double E = 200e9;
    double nu = 0.3;
    double t = 0.01;
    double eye = t*t*t/12.0;

    Vector3d A(0,0,0), B(1,0,0), C(1,1,0), D(0,1,0);

    int NPIVOT = 1;
    int JNOT = (NPIVOT - 2 <= 0) ? NPIVOT + 2 : NPIVOT - 2;

    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};

    Vector3d D1 = C - A;
    Vector3d D2 = D - B;
    Vector3d A1 = B - A;

    Vector3d KVECT = D1.cross(D2);
    KVECT /= KVECT.norm();

    double temp = A1.dot(KVECT) / 2.0;
    Vector3d IVECT = A1 - temp * KVECT;
    IVECT /= IVECT.norm();

    Vector3d JVECT = KVECT.cross(IVECT);
    JVECT /= JVECT.norm();

    MatrixXd R = MatrixXd::Zero(2, 4);
    R(0,2) = D1.dot(IVECT);
    R(0,1) = A1.dot(IVECT);
    R(1,2) = D1.dot(JVECT);
    R(0,3) = D2.dot(IVECT) + R(0,1);
    R(1,3) = D2.dot(JVECT);

    MatrixXd KSUM = MatrixXd::Zero(12, 3);

    for (int J = 1; J <= 4; J++) {
        if (J == JNOT) continue;
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0];
        int SUBSCB = M[KM+1];
        int SUBSCC = M[KM+2];

        Vector2d V, VV;
        for (int i = 0; i < 2; i++) {
            V(i) = R(i, SUBSCB-1) - R(i, SUBSCA-1);
            VV(i) = R(i, SUBSCC-1) - R(i, SUBSCA-1);
        }
        double XSUBB = V.norm();
        double U1 = V(0)/XSUBB;
        double U2 = V(1)/XSUBB;
        double XSUBC = U1*VV(0) + U2*VV(1);
        double YSUBC = U1*VV(1) - U2*VV(0);

        KTRBSCResult kr = KTRBSC(XSUBB, XSUBC, YSUBC, E, nu, t, eye, 0.0);
        if (kr.nogo) { cerr << "NOGO" << endl; return 1; }

        Matrix3d T;
        T << 1.0, 0.0, 0.0,
             0.0, U1, U2,
             0.0, -U2, U1;

        int NBEGIN = 0;
        for (int i = 1; i <= 3; i++) {
            int NPOINT = KM + i;
            if (M[NPOINT-1] == NPIVOT) {
                NBEGIN = 27*i - 27;
                break;
            }
        }
        int pivot_row = NBEGIN / 27;

        for (int i = 1; i <= 3; i++) {
            Matrix3d block = kr.KU.block<3,3>(pivot_row*3, (i-1)*3);
            Matrix3d TEMP9 = T.transpose() * block;
            Matrix3d PROD9 = TEMP9 * T;

            int NPOINT = KM + i;
            int node = M[NPOINT-1];
            KSUM.block<3,3>((node-1)*3, 0) += PROD9 / 2.0;
        }
    }

    MatrixXd E6 = MatrixXd::Zero(6, 3);
    E6(0,0) = KVECT(0); E6(1,0) = KVECT(1); E6(2,0) = KVECT(2);
    E6(3,1) = IVECT(0); E6(4,1) = IVECT(1); E6(5,1) = IVECT(2);
    E6(3,2) = JVECT(0); E6(4,2) = JVECT(1); E6(5,2) = JVECT(2);

    MatrixXd TITE = E6;
    MatrixXd K12 = MatrixXd::Zero(12, 12);

    for (int J = 1; J <= 4; J++) {
        MatrixXd TJTE = E6;
        Matrix3d KsumBlock = KSUM.block<3,3>((J-1)*3, 0);
        MatrixXd TEMP18 = KsumBlock * TJTE.transpose();
        MatrixXd KOUT = TITE * TEMP18;

        int pi = 0;
        int ji = (J-1)*3;

        K12.block<3,3>(pi, pi) += KOUT.block<3,3>(0,0);
        K12.block<3,3>(pi, ji) += KOUT.block<3,3>(0,3);
        K12.block<3,3>(ji, pi) += KOUT.block<3,3>(3,0);
        K12.block<3,3>(ji, ji) += KOUT.block<3,3>(3,3);
    }

    cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        cout << "[";
        for (int j = 0; j < 12; j++) {
            cout << scientific << setprecision(6) << K12(i,j);
            if (j < 11) cout << ",";
        }
        cout << "]";
        if (i < 11) cout << ",";
    }
    cout << "]}" << endl;

    return 0;
}