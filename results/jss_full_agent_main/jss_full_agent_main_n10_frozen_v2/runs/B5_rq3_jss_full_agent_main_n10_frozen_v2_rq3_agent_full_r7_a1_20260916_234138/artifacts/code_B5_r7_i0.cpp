#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

static const double PI = 3.14159265358979323846;
static const double DEGRA = PI/180.0;

static MatrixXd inverd(const MatrixXd& M, int& ising) {
    int n = M.rows();
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = 2; return MatrixXd::Zero(n,n); }
    ising = 1;
    return lu.inverse();
}

struct KTRBSCResult {
    MatrixXd KU[3][3];
    MatrixXd HINV;
    MatrixXd S;
};

static KTRBSCResult ktrbsc(double x1,double y1,double z1,
                           double x2,double y2,double z2,
                           double x3,double y3,double z3,
                           double E, double nu, double t,
                           double angle_deg, double eye,
                           double T2, double G2X211, double G2X212, double G2X222,
                           int iopt) {
    KTRBSCResult res;

    Vector3d d2(x2-x1, y2-y1, z2-z1);
    double XSUBB = d2.norm();
    Vector3d I = d2 / XSUBB;

    Vector3d d1(x3-x1, y3-y1, z3-z1);
    double XSUBC = I.dot(d1);

    Vector3d Kraw = I.cross(d1);
    double YSUBC = Kraw.norm();
    Vector3d K = Kraw / YSUBC;

    Vector3d J = K.cross(I);
    J.normalize();

    double THETA = angle_deg * DEGRA;
    double SINTH = sin(THETA);
    double COSTH = cos(THETA);
    if (fabs(SINTH) < 1.0e-6) SINTH = 0.0;

    double D0 = E * t*t*t / (12.0*(1.0 - nu*nu));
    double G11 = D0, G22 = D0, G33 = D0*(1.0-nu)/2.0;
    double G12 = nu*D0, G13 = 0.0, G23 = 0.0;

    Matrix3d G;
    G << G11, G12, G13,
         G12, G22, G23,
         G13, G23, G33;

    Matrix3d D = G * eye;

    double AREA = XSUBB*YSUBC/2.0;
    double XBAR = (XSUBB+XSUBC)/3.0;
    double YBAR = YSUBC/3.0;

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

    auto Dv = [&](int i,int j)->double{ return D(i-1,j-1); };

    MatrixXd KX = MatrixXd::Zero(6,6);
    KX(0,0)=Dv(1,1);
    KX(0,1)=Dv(3,1);
    KX(0,2)=Dv(2,1);
    KX(0,3)=Dv(1,1)*XBAR3;
    KX(0,4)=Dv(2,1)*XBAR + YBAR2*Dv(3,1);
    KX(0,5)=Dv(2,1)*YBAR3;
    KX(1,0)=KX(0,1);
    KX(1,1)=Dv(3,3);
    KX(1,2)=Dv(3,2);
    KX(1,3)=Dv(3,1)*XBAR3;
    KX(1,4)=Dv(3,2)*XBAR + YBAR2*Dv(3,3);
    KX(1,5)=Dv(3,2)*YBAR3;
    KX(2,0)=KX(0,2);
    KX(2,1)=KX(1,2);
    KX(2,2)=Dv(2,2);
    KX(2,3)=Dv(2,1)*XBAR3;
    KX(2,4)=Dv(2,2)*XBAR + YBAR2*Dv(3,2);
    KX(2,5)=Dv(2,2)*YBAR3;
    KX(3,0)=KX(0,3);
    KX(3,1)=KX(1,3);
    KX(3,2)=KX(2,3);
    KX(3,3)=Dv(1,1)*9.0*PX2;
    KX(3,4)=Dv(2,1)*3.0*PX2 + 6.0*PXY2*Dv(3,1);
    KX(3,5)=Dv(2,1)*9.0*PXY2;
    KX(4,0)=KX(0,4);
    KX(4,1)=KX(1,4);
    KX(4,2)=KX(2,4);
    KX(4,3)=KX(3,4);
    KX(4,4)=Dv(2,2)*PX2 + 4.0*PXY2*Dv(3,2) + 4.0*PY2*Dv(3,3);
    KX(4,5)=Dv(2,2)*3.0*PXY2 + 6.0*PY2*Dv(3,2);
    KX(5,0)=KX(0,5);
    KX(5,1)=KX(1,5);
    KX(5,2)=KX(2,5);
    KX(5,3)=KX(3,5);
    KX(5,4)=KX(4,5);
    KX(5,5)=Dv(2,2)*9.0*PY2;
    double TEMP = 4.0*AREA;
    KX *= TEMP;

    MatrixXd HB = MatrixXd::Zero(6,6);
    HB(0,0)=XBSQ;
    HB(0,3)=XBSQ*XSUBB;
    HB(0,4)=XSUBB;
    HB(0,5)=-2.0*XSUBB;
    HB(1,0)=-3.0*XBSQ;
    HB(2,2)=XCSQ;
    HB(2,3)=XCYC;
    HB(2,4)=YCSQ;
    HB(2,5)=XCSQ*XSUBC;
    HB(3,2)=YCSQ*XSUBC;
    HB(3,3)=YCSQ*YSUBC;
    HB(3,4)=XSUBC;
    HB(3,5)=YSUBC*2.0;
    HB(4,2)=XCYC*2.0;
    HB(4,3)=YCSQ*3.0;
    HB(4,4)=-2.0*XSUBC;
    HB(4,5)=-YSUBC;
    HB(5,2)=-3.0*XCSQ;
    HB(5,3)=-YCSQ;

    if (T2 != 0.0) {
        if (!(G2X211==0.0 && G2X212==0.0 && G2X222==0.0)) {
            Matrix2d G2X2;
            G2X2 << G2X211*T2, G2X212*T2,
                    G2X212*T2, G2X222*T2;
            double DETERM = G2X2(0,0)*G2X2(1,1) - G2X2(1,0)*G2X2(0,1);
            Matrix2d J2X2;
            J2X2 << G2X2(1,1)/DETERM, -G2X2(0,1)/DETERM,
                   -G2X2(1,0)/DETERM,  G2X2(0,0)/DETERM;
            double TEMP2 = 2.0*Dv(2,1) + 4.0*Dv(3,3);
            MatrixXd HYQ(2,3);
            HYQ(0,0) = -6.0*(J2X2(0,0)*Dv(1,1) + J2X2(0,1)*Dv(3,1));
            HYQ(0,1) = -J2X2(0,0)*TEMP2 - 6.0*J2X2(0,1)*Dv(3,2);
            HYQ(0,2) = -6.0*(J2X2(0,0)*Dv(3,2) + J2X2(0,1)*Dv(2,2));
            HYQ(1,0) = -6.0*(J2X2(1,0)*Dv(1,1) + J2X2(1,1)*Dv(3,1));
            HYQ(1,1) = -J2X2(1,0)*TEMP2 - 6.0*J2X2(1,1)*Dv(3,2);
            HYQ(1,2) = -6.0*(J2X2(1,0)*Dv(3,2) + J2X2(1,1)*Dv(2,2));
            MatrixXd prod = G2X2 * HYQ;
            MatrixXd lower = HYQ.transpose() * prod;
            for (int i=0;i<3;i++){
                KX(3+i,3) += lower(i,0)*AREA;
                KX(3+i,4) += lower(i,1)*AREA;
                KX(3+i,5) += lower(i,2)*AREA;
            }
            for (int i=0;i<3;i++){
                HB(1+i,0) += XSUBB*HYQ(0,i);
                HB(3+i,2) += XSUBC*HYQ(0,i) + YSUBC*HYQ(1,i);
            }
        }
    }

    int ising;
    MatrixXd HINV = inverd(HB, ising);

    MatrixXd KQ = KX * HINV;
    MatrixXd KII = HINV.transpose() * KQ;

    MatrixXd S(6,3);
    S << 1.0, 0.0, -XSUBB,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, YSUBC, -XSUBC,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

    MatrixXd KIA = -(KII * S);
    MatrixXd KAA = S.transpose() * (-KIA);

    MatrixXd KIA_top = KIA.topRows(3);
    MatrixXd KIA_bot = KIA.bottomRows(3);

    res.KU[0][0] = KAA;
    res.KU[0][1] = KIA_top.transpose();
    res.KU[0][2] = KIA_bot.transpose();
    res.KU[1][0] = KIA_top;
    res.KU[1][1] = KII.topLeftCorner(3,3);
    res.KU[1][2] = KII.topRightCorner(3,3);
    res.KU[2][0] = KIA_bot;
    res.KU[2][1] = KII.bottomLeftCorner(3,3);
    res.KU[2][2] = KII.bottomRightCorner(3,3);

    res.HINV = HINV;
    res.S = S;
    return res;
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double eye = t*t*t/12.0;

    double x1=0,y1=0,z1=0;
    double x2=1,y2=0,z2=0;
    double x3=0,y3=1,z3=0;

    int NPVT = 1;

    Vector3d V1(x1,y1,z1), V2(x2,y2,z2), V3(x3,y3,z3);
    Vector3d D2 = V2 - V1;
    Vector3d D1 = V3 - V1;

    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,1) = D2.norm();
    Vector3d IVECT = D2 / R(0,1);
    Vector3d KVECT = IVECT.cross(D1);
    R(1,2) = KVECT.norm();
    KVECT /= R(1,2);
    Vector3d JVECT = KVECT.cross(IVECT);
    JVECT.normalize();
    R(0,2) = D1.dot(IVECT);
    R(0,3) = (R(0,1)+R(0,2))/3.0;
    R(1,3) = R(1,2)/3.0;

    int M[9] = {1,2,4, 2,3,4, 3,1,4};

    MatrixXd KSUM = MatrixXd::Zero(7,9);
    MatrixXd G = MatrixXd::Zero(4,9);

    double THETA = 0.0;
    double SINANG = sin(THETA);
    double COSANG = cos(THETA);

    for (int J=1; J<=3; J++) {
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0];
        int SUBSCB = M[KM+1];
        int SUBSCC = M[KM+2];

        Vector2d V, VV;
        for (int i=0;i<2;i++){
            V(i)  = R(i,SUBSCB-1) - R(i,SUBSCA-1);
            VV(i) = R(i,SUBSCC-1) - R(i,SUBSCA-1);
        }
        double XSUBB = V.norm();
        double U1 = V(0)/XSUBB;
        double U2 = V(1)/XSUBB;
        double XSUBC = U1*VV(0) + U2*VV(1);
        double YSUBC = U1*VV(1) - U2*VV(0);

        double SINTH = SINANG*U1 - COSANG*U2;
        double COSTH = COSANG*U1 + SINANG*U2;
        if (fabs(SINTH) < 1.0e-6) SINTH = 0.0;

        KTRBSCResult r = ktrbsc(0,0,0, XSUBB,0,0, XSUBC,YSUBC,0,
                                E, nu, t, 0.0, eye, 0.0,
                                0.0,0.0,0.0, 2);

        Matrix3d T;
        T << 1.0, 0.0, 0.0,
             0.0, U1, U2,
             0.0, -U2, U1;

        for (int I=1; I<=3; I++) {
            MatrixXd prod = T.transpose() * r.KU[I-1][I-1] * T;
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] + 18;
            for (int k=0;k<9;k++){
                int row = (NPOINT-1+k)/9;
                int col = (NPOINT-1+k)%9;
                KSUM(row,col) += prod(k/3, k%3);
            }
        }
        for (int K=1; K<=2; K++) {
            int NPOINT = KM + K;
            if (M[NPOINT-1] != NPVT) continue;
            MatrixXd prod = T.transpose() * r.KU[K-1][2] * T;
            int NP = 9*NPVT - 9;
            for (int i=0;i<9;i++){
                int row=(NP+i)/9, col=(NP+i)%9;
                KSUM(row,col) += prod(i/3,i%3);
            }
            MatrixXd prod2 = T.transpose() * r.KU[2][K-1] * T;
            int NP2 = KM + 3 - K;
            NP2 = 9*M[NP2-1] - 9;
            for (int i=0;i<9;i++){
                int row=(NP2+i)/9, col=(NP2+i)%9;
                KSUM(row,col) += prod2(i/3,i%3);
            }
        }

        double TEMP1 = XSUBB - XSUBC;
        double TEMP2 = YSUBC*YSUBC;
        double L1 = sqrt(XSUBC*XSUBC + TEMP2);
        double L2 = sqrt(TEMP1*TEMP1 + TEMP2);
        double S1 = XSUBC/L1;
        double S2 = TEMP1/L2;
        double C1 = YSUBC/L1;
        double C2 = YSUBC/L2;
        double X1 = XSUBC/2.0;
        double Y1 = YSUBC/2.0;
        double X2 = (XSUBB+XSUBC)/2.0;
        double Y2 = Y1;

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
        HABC(0,0) += S1;
        HABC(0,1) += C1;
        HABC(1,0) += S2;
        HABC(1,1) -= C2;

        MatrixXd HpsiA = HABC;
        MatrixXd HpsiB(2,3), HpsiC(2,3);
        HpsiB << PROD12(0,0),PROD12(0,1),PROD12(0,2),
                 PROD12(1,0),PROD12(1,1),PROD12(1,2);
        HpsiC << PROD12(0,3),PROD12(0,4),PROD12(0,5),
                 PROD12(1,3),PROD12(1,4),PROD12(1,5);

        MatrixXd Hsub[3] = {HpsiA, HpsiB, HpsiC};

        for (int I=1; I<=3; I++) {
            MatrixXd TEMP9 = Hsub[I-1] * T;
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] - 9;
            if (J == 2) {
                NPOINT += 3;
                for (int k=0;k<6;k++){
                    int row=(NPOINT+k)/9, col=(NPOINT+k)%9;
                    G(row,col) += TEMP9(k/3,k%3);
                }
            } else if (J == 1) {
                for (int k=0;k<6;k++){
                    int row=(NPOINT+k)/9, col=(NPOINT+k)%9;
                    G(row,col) += TEMP9(k/3,k%3);
                }
            } else {
                G((NPOINT+7)/9,(NPOINT+7)%9) += TEMP9(0,0);
                G((NPOINT+8)/9,(NPOINT+8)%9) += TEMP9(0,1);
                G((NPOINT+9)/9,(NPOINT+9)%9) += TEMP9(0,2);
                G((NPOINT+1)/9,(NPOINT+1)%9) += TEMP9(1,0);
                G((NPOINT+2)/9,(NPOINT+2)%9) += TEMP9(1,1);
                G((NPOINT+3)/9,(NPOINT+3)%9) += TEMP9(1,2);
            }
        }
    }

    MatrixXd Emat = MatrixXd::Zero(6,3);
    Emat(0,0)=KVECT(0); Emat(3,0)=KVECT(1); Emat(6,0)=KVECT(2);
    Emat(4,1)=IVECT(0); Emat(7,1)=IVECT(1); Emat(10,1)=IVECT(2);
    Emat(5,2)=JVECT(0); Emat(8,2)=JVECT(1); Emat(11,2)=JVECT(2);

    MatrixXd TITE = Emat;

    MatrixXd G4m(3,3);
    for (int i=0;i<9;i++) G4m(i/3,i%3) = G((27+i)/9,(27+i)%9);

    int ising;
    MatrixXd G4inv = inverd(G4m, ising);

    MatrixXd Gpiv(3,3);
    int base = 9*NPVT - 9;
    for (int i=0;i<9;i++) Gpiv(i/3,i%3) = G((base+i)/9,(base+i)%9);
    MatrixXd PROD9 = G4inv * Gpiv;

    MatrixXd K44(3,3);
    for (int i=0;i<9;i++) K44(i/3,i%3) = KSUM((54+i)/9,(54+i)%9);
    MatrixXd TEMP9 = PROD9.transpose() * K44;

    MatrixXd KOUTs[3];

    for (int J=1; J<=3; J++) {
        MatrixXd KJ4(3,3);
        int b1 = 9*J+19-1;
        for (int i=0;i<9;i++) KJ4(i/3,i%3) = KSUM((b1+i)/9,(b1+i)%9);
        MatrixXd ARR9 = PROD9.transpose() * KJ4.transpose();

        int NBEGIN = 9*J-9;
        for (int i=0;i<9;i++){
            int row=(NBEGIN+i)/9, col=(NBEGIN+i)%9;
            KSUM(row,col) -= ARR9(i/3,i%3);
        }

        MatrixXd GJ(3,3);
        int b2 = 9*J-8-1;
        for (int i=0;i<9;i++) GJ(i/3,i%3) = G((b2+i)/9,(b2+i)%9);
        ARR9 = G4inv * GJ;

        MatrixXd KI4(3,3);
        int b3 = 9*NPVT+19-1;
        for (int i=0;i<9;i++) KI4(i/3,i%3) = KSUM((b3+i)/9,(b3+i)%9);
        MatrixXd ARRAY9 = KI4 * ARR9;

        for (int i=0;i<9;i++){
            int row=(NBEGIN+i)/9, col=(NBEGIN+i)%9;
            KSUM(row,col) -= ARRAY9(i/3,i%3);
        }

        ARRAY9 = TEMP9 * ARR9;
        for (int i=0;i<9;i++){
            int row=(NBEGIN+i)/9, col=(NBEGIN+i)%9;
            KSUM(row,col) += ARRAY9(i/3,i%3);
        }

        MatrixXd TJTE = Emat;
        MatrixXd KIJ(3,3);
        for (int i=0;i<9;i++) KIJ(i/3,i%3) = KSUM((NBEGIN+i)/9,(NBEGIN+i)%9);
        MatrixXd TEMP18 = KIJ * TJTE.transpose();
        MatrixXd KOUT = TITE * TEMP18;
        KOUTs[J-1] = KOUT;
    }

    // Assemble 9x9. Pivot is node NPVT=1 (A). Node J maps to grid NECPT(J+1).
    // For test: nodes A,B,C = 1,2,3. KOUT is 6x6 = (pivot 3 DOF, node J 3 DOF).
    // Place into 9x9: pivot block at rows 0-2, node J block at rows 3*(J-1)..3*(J-1)+2.
    MatrixXd K9 = MatrixXd::Zero(9,9);
    for (int J=1; J<=3; J++) {
        MatrixXd KOUT = KOUTs[J-1];
        int pbase = 3*(NPVT-1);
        int jbase = 3*(J-1);
        // KOUT(0:3,0:3) -> pivot-pivot
        K9.block(pbase, pbase, 3,3) += KOUT.block(0,0,3,3);
        // KOUT(0:3,3:6) -> pivot-nodeJ
        K9.block(pbase, jbase, 3,3) += KOUT.block(0,3,3,3);
        // KOUT(3:6,0:3) -> nodeJ-pivot
        K9.block(jbase, pbase, 3,3) += KOUT.block(3,0,3,3);
        // KOUT(3:6,3:6) -> nodeJ-nodeJ
        K9.block(jbase, jbase, 3,3) += KOUT.block(3,3,3,3);
    }

    cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<9;i++){
        cout << "[";
        for (int j=0;j<9;j++){
            cout << scientific << setprecision(6) << K9(i,j);
            if (j<8) cout << ",";
        }
        cout << "]";
        if (i<8) cout << ",";
    }
    cout << "]}" << endl;

    return 0;
}