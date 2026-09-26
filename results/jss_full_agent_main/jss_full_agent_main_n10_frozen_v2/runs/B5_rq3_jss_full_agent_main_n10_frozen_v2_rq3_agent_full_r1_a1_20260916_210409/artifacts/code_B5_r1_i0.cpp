#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

double G11, G12, G13, G22, G23, G33;

void MAT_setup(double E, double nu) {
    double factor = E / (1.0 - nu*nu);
    G11 = factor;
    G12 = factor * nu;
    G13 = 0.0;
    G22 = factor;
    G23 = 0.0;
    G33 = factor * (1.0 - nu) / 2.0;
}

struct KTRBSCResult {
    Matrix3d KAA;
    MatrixXd KIA;  // 6x3
    MatrixXd KII;  // 6x6
    MatrixXd Hinv; // 6x6
    MatrixXd S;    // 6x3
};

KTRBSCResult KTRBSC(double XSUBB, double XSUBC, double YSUBC, double EYE) {
    KTRBSCResult res;
    Matrix3d G;
    G << G11, G12, G13,
         G12, G22, G23,
         G13, G23, G33;
    Matrix3d D = G * EYE;

    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    double XCSQ = XSUBC*XSUBC, YCSQ = YSUBC*YSUBC, XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ + XSUBB*XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC * (XSUBB + 2.0*XSUBC) / 12.0;
    double XBAR3 = 3.0*XBAR, YBAR3 = 3.0*YBAR, YBAR2 = 2.0*YBAR;

    MatrixXd KX(6,6);
    KX(0,0)=D(0,0); KX(0,1)=D(0,2); KX(0,2)=D(0,1);
    KX(0,3)=D(0,0)*XBAR3; KX(0,4)=D(0,1)*XBAR+YBAR2*D(0,2); KX(0,5)=D(0,1)*YBAR3;
    KX(1,0)=KX(0,1); KX(1,1)=D(2,2); KX(1,2)=D(1,2);
    KX(1,3)=D(0,2)*XBAR3; KX(1,4)=D(1,2)*XBAR+YBAR2*D(2,2); KX(1,5)=D(1,2)*YBAR3;
    KX(2,0)=KX(0,2); KX(2,1)=KX(1,2); KX(2,2)=D(1,1);
    KX(2,3)=D(0,1)*XBAR3; KX(2,4)=D(1,1)*XBAR+YBAR2*D(1,2); KX(2,5)=D(1,1)*YBAR3;
    KX(3,0)=KX(0,3); KX(3,1)=KX(1,3); KX(3,2)=KX(2,3);
    KX(3,3)=D(0,0)*9.0*PX2; KX(3,4)=D(0,1)*3.0*PX2+6.0*PXY2*D(0,2); KX(3,5)=D(0,1)*9.0*PXY2;
    KX(4,0)=KX(0,4); KX(4,1)=KX(1,4); KX(4,2)=KX(2,4); KX(4,3)=KX(3,4);
    KX(4,4)=D(1,1)*PX2+4.0*PXY2*D(1,2)+4.0*PY2*D(2,2);
    KX(4,5)=D(1,1)*3.0*PXY2+6.0*PY2*D(1,2);
    KX(5,0)=KX(0,5); KX(5,1)=KX(1,5); KX(5,2)=KX(2,5); KX(5,3)=KX(3,5); KX(5,4)=KX(4,5);
    KX(5,5)=D(1,1)*9.0*PY2;
    KX *= 4.0*AREA;

    MatrixXd HBAR = MatrixXd::Zero(6,6);
    auto setA = [&](int idx, double val) {
        int i = idx - 37;
        HBAR(i % 6, i / 6) = val;
    };
    setA(37, XBSQ); setA(40, XBSQ*XSUBB); setA(44, XSUBB);
    setA(49, -2.0*XSUBB); setA(52, -3.0*XBSQ);
    setA(55, XCSQ); setA(56, XCYC); setA(57, YCSQ);
    setA(58, XCSQ*XSUBC); setA(59, YCSQ*XSUBC); setA(60, YCSQ*YSUBC);
    setA(62, XSUBC); setA(63, YSUBC*2.0);
    setA(65, XCYC*2.0); setA(66, YCSQ*3.0);
    setA(67, -2.0*XSUBC); setA(68, -YSUBC);
    setA(70, -3.0*XCSQ); setA(71, -YCSQ);

    MatrixXd Hinv = HBAR.inverse();
    MatrixXd KQHinv = KX * Hinv;
    MatrixXd KII = Hinv.transpose() * KQHinv;

    MatrixXd S(6,3);
    S << 1.0, 0.0, -XSUBB,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, YSUBC, -XSUBC,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

    MatrixXd KIA_raw = KII * S;
    Matrix3d KAA = S.transpose() * KIA_raw;
    MatrixXd KIA = -KIA_raw;

    res.KAA = KAA;
    res.KIA = KIA;
    res.KII = KII;
    res.Hinv = Hinv;
    res.S = S;
    return res;
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double I = t*t*t/12.0;
    MAT_setup(E, nu);

    Vector3d vA(0,0,0), vB(1,0,0), vC(0,1,0);
    Vector3d D2 = vB - vA, D1 = vC - vA;
    double R12 = D2.norm();
    Vector3d IVECT = D2 / R12;
    Vector3d KVECT;
    KVECT(0) = IVECT(1)*D1(2) - D1(1)*IVECT(2);
    KVECT(1) = IVECT(2)*D1(0) - D1(2)*IVECT(0);
    KVECT(2) = IVECT(0)*D1(1) - D1(0)*IVECT(1);
    double R23 = KVECT.norm();
    KVECT /= R23;
    Vector3d JVECT;
    JVECT(0) = KVECT(1)*IVECT(2) - IVECT(1)*KVECT(2);
    JVECT(1) = KVECT(2)*IVECT(0) - IVECT(2)*KVECT(0);
    JVECT(2) = KVECT(0)*IVECT(1) - IVECT(0)*KVECT(1);
    JVECT /= JVECT.norm();
    double R13 = D1.dot(IVECT);

    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,1)=R12; R(1,2)=R23; R(0,2)=R13;
    R(0,3)=(R12+R13)/3.0; R(1,3)=R23/3.0;

    int M[9] = {1,2,4, 2,3,4, 3,1,4};
    VectorXd KSUM = VectorXd::Zero(63);
    VectorXd G = VectorXd::Zero(36);

    double SINANG = 0.0, COSANG = 1.0;
    int NPIVOT = 1;

    for (int J = 1; J <= 3; J++) {
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0], SUBSCB = M[KM+1], SUBSCC = M[KM+2];
        double V[2], VV[2];
        for (int i = 0; i < 2; i++) {
            V[i] = R(i, SUBSCB-1) - R(i, SUBSCA-1);
            VV[i] = R(i, SUBSCC-1) - R(i, SUBSCA-1);
        }
        double XSUBB = sqrt(V[0]*V[0]+V[1]*V[1]);
        double U1 = V[0]/XSUBB, U2 = V[1]/XSUBB;
        double XSUBC = U1*VV[0]+U2*VV[1];
        double YSUBC = U1*VV[1]-U2*VV[0];

        KTRBSCResult r = KTRBSC(XSUBB, XSUBC, YSUBC, I);

        Matrix3d T;
        T << 1.0, 0.0, 0.0,
             0.0, U1, U2,
             0.0, -U2, U1;

        Matrix3d K11 = r.KAA;
        Matrix3d K21 = r.KIA.block<3,3>(0,0);
        Matrix3d K31 = r.KIA.block<3,3>(3,0);
        Matrix3d K12 = K21.transpose();
        Matrix3d K13 = K31.transpose();
        Matrix3d K22 = r.KII.block<3,3>(0,0);
        Matrix3d K23 = r.KII.block<3,3>(0,3);
        Matrix3d K32 = r.KII.block<3,3>(3,0);
        Matrix3d K33 = r.KII.block<3,3>(3,3);

        Matrix3d Kblk[3][3] = {{K11,K12,K13},{K21,K22,K23},{K31,K32,K33}};

        // Diagonal blocks: T^T * K_II * T -> KSUM at 9*M(KM+I)+18
        for (int I = 1; I <= 3; I++) {
            Matrix3d prod = T.transpose() * Kblk[I-1][I-1] * T;
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] + 18;
            for (int K = 1; K <= 9; K++) {
                KSUM(NPOINT + K - 1) += prod((K-1)%3, (K-1)/3);
            }
        }

        // Off-diagonal with pivot
        for (int K = 1; K <= 2; K++) {
            int NPOINT = KM + K;
            if (M[NPOINT-1] != NPIVOT) continue;
            Matrix3d prod = T.transpose() * Kblk[K-1][0] * T;
            int NP = 9*NPIVOT - 9;
            for (int i = 1; i <= 9; i++) {
                KSUM(NP + i - 1) += prod((i-1)%3, (i-1)/3);
            }
            Matrix3d prod2 = T.transpose() * Kblk[0][K-1] * T;
            int NP2 = KM + 3 - K;
            NP2 = 9*M[NP2-1] - 9;
            for (int i = 1; i <= 9; i++) {
                KSUM(NP2 + i - 1) += prod2((i-1)%3, (i-1)/3);
            }
        }

        // HQ matrix (2x6)
        double TEMP1 = XSUBB - XSUBC;
        double TEMP2 = YSUBC*YSUBC;
        double L1 = sqrt(XSUBC*XSUBC + TEMP2);
        double L2 = sqrt(TEMP1*TEMP1 + TEMP2);
        double S1 = XSUBC/L1, S2 = TEMP1/L2;
        double C1 = YSUBC/L1, C2 = YSUBC/L2;
        double X1 = XSUBC/2.0, Y1 = YSUBC/2.0;
        double X2 = (XSUBB+XSUBC)/2.0, Y2 = Y1;

        MatrixXd HQ(2,6);
        HQ(0,0)=-XSUBC*C1; HQ(0,1)=X1*S1-Y1*C1; HQ(0,2)=2.0*Y1*S1;
        HQ(0,3)=-3.0*X1*X1*C1; HQ(0,4)=Y1*(2.0*X1*S1-Y1*C1); HQ(0,5)=3.0*Y1*Y1*S1;
        HQ(1,0)=2.0*X2*C2; HQ(1,1)=X2*S2+Y2*C2; HQ(1,2)=2.0*Y2*S2;
        HQ(1,3)=3.0*X2*X2*C2; HQ(1,4)=Y2*(2.0*X2*S2+Y2*C2); HQ(1,5)=3.0*Y2*Y2*S2;

        MatrixXd PROD12 = HQ * r.Hinv;       // 2x6
        MatrixXd HABC = -PROD12 * r.S;       // 2x3
        HABC(0,1) += S1; HABC(0,2) += C1;
        HABC(1,1) += S2; HABC(1,2) -= C2;

        // Build HABCv2 (18 elements)
        VectorXd HABCv2(18);
        HABCv2(0) = -HABC(0,0);
        HABCv2(1) = -HABC(1,0) + S1;
        HABCv2(2) = -HABC(0,1) + C1;
        HABCv2(3) = -HABC(1,1);
        HABCv2(4) = -HABC(0,2) + S2;
        HABCv2(5) = -HABC(1,2) - C2;
        HABCv2(6) = PROD12(0,0); HABCv2(7) = PROD12(1,0); HABCv2(8) = PROD12(0,1);
        HABCv2(9) = PROD12(0,3); HABCv2(10) = PROD12(1,3); HABCv2(11) = PROD12(0,4);
        HABCv2(12) = PROD12(1,1); HABCv2(13) = PROD12(0,2); HABCv2(14) = PROD12(1,2);
        HABCv2(15) = PROD12(1,4); HABCv2(16) = PROD12(0,5); HABCv2(17) = PROD12(1,5);

        for (int I = 1; I <= 3; I++) {
            MatrixXd Hblk(2,3);
            Hblk(0,0)=HABCv2(6*I-6); Hblk(1,0)=HABCv2(6*I-5); Hblk(0,1)=HABCv2(6*I-4);
            Hblk(1,1)=HABCv2(6*I-3); Hblk(0,2)=HABCv2(6*I-2); Hblk(1,2)=HABCv2(6*I-1);
            MatrixXd TEMP9 = Hblk * T;  // 2x3

            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] - 9;

            if (J == 2) NPOINT += 3;
            if (J == 1 || J == 2) {
                for (int K = 1; K <= 6; K++) {
                    NPOINT += 1;
                    G(NPOINT-1) += TEMP9((K-1)%2, (K-1)/2);
                }
            } else {
                G(NPOINT+7-1) += TEMP9(0,0);
                G(NPOINT+8-1) += TEMP9(1,0);
                G(NPOINT+9-1) += TEMP9(0,1);
                G(NPOINT+1-1) += TEMP9(1,1);
                G(NPOINT+2-1) += TEMP9(0,2);
                G(NPOINT+3-1) += TEMP9(1,2);
            }
        }
    }

    // G4 = G(28..36)
    Matrix3d G4;
    for (int i = 0; i < 9; i++) G4(i%3, i/3) = G(27+i);
    Matrix3d G4inv = G4.inverse();

    // PROD9 = G4inv * G(9*NPIVOT-8) = G4inv * G_1
    Matrix3d G1;
    for (int i = 0; i < 9; i++) G1(i%3, i/3) = G(0+i);
    Matrix3d PROD9 = G4inv * G1;

    // TEMP9 = PROD9^T * KSUM(55..63)
    Matrix3d K44;
    for (int i = 0; i < 9; i++) K44(i%3, i/3) = KSUM(54+i);
    Matrix3d TEMP9 = PROD9.transpose() * K44;

    for (int J = 1; J <= 3; J++) {
        Matrix3d KJ4;
        for (int i = 0; i < 9; i++) KJ4(i%3, i/3) = KSUM(9*J+18+i);
        Matrix3d ARR9 = PROD9.transpose() * KJ4;

        int NBEGIN = 9*J - 9;
        for (int i = 0; i < 9; i++) KSUM(NBEGIN+i) -= ARR9(i%3, i/3);

        Matrix3d GJ;
        for (int i = 0; i < 9; i++) GJ(i%3, i/3) = G(9*J-9+i);
        ARR9 = G4inv * GJ;

        Matrix3d K4J;
        for (int i = 0; i < 9; i++) K4J(i%3, i/3) = KSUM(9*NPIVOT+18+i);
        Matrix3d ARRAY9 = K4J * ARR9;
        for (int i = 0; i < 9; i++) KSUM(NBEGIN+i) -= ARRAY9(i%3, i/3);

        ARRAY9 = TEMP9 * ARR9;
        for (int i = 0; i < 9; i++) KSUM(NBEGIN+i) += ARRAY9(i%3, i/3);
    }

    // Assemble 9x9
    auto getBlk = [&](int start) {
        Matrix3d B;
        for (int i = 0; i < 9; i++) B(i%3, i/3) = KSUM(start-1+i);
        return B;
    };
    Matrix3d K11 = getBlk(1), K12 = getBlk(10), K13 = getBlk(19);
    Matrix3d K22 = getBlk(37), K23 = getBlk(46);
    Matrix3d K33 = getBlk(28);

    MatrixXd K9 = MatrixXd::Zero(9,9);
    K9.block<3,3>(0,0) = K11;
    K9.block<3,3>(0,3) = K12;
    K9.block<3,3>(0,6) = K13;
    K9.block<3,3>(3,0) = K12.transpose();
    K9.block<3,3>(3,3) = K22;
    K9.block<3,3>(3,6) = K23;
    K9.block<3,3>(6,0) = K13.transpose();
    K9.block<3,3>(6,3) = K23.transpose();
    K9.block<3,3>(6,6) = K33;

    cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 9; i++) {
        cout << "[";
        for (int j = 0; j < 9; j++) {
            cout << scientific << setprecision(6) << K9(i,j);
            if (j < 8) cout << ",";
        }
        cout << "]";
        if (i < 8) cout << ",";
    }
    cout << "]}" << endl;

    return 0;
}