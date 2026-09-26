#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

double G11, G12, G13, G22, G23, G33;
double G2X211, G2X212, G2X222;
double E_young, nu_val, thickness;

void MAT_isotropic(double E, double nu) {
    double factor = E / (1.0 - nu*nu);
    G11 = factor; G12 = factor*nu; G13 = 0.0;
    G22 = factor; G23 = 0.0; G33 = factor*(1.0-nu)/2.0;
}

MatrixXd inverd(const MatrixXd& M, bool& ising) {
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = true; return MatrixXd::Zero(M.rows(), M.cols()); }
    ising = false;
    return lu.inverse();
}

void KTRBSC(double XSUBB, double XSUBC, double YSUBC, double EYE, double T2,
            MatrixXd KsuperU[3][3]) {
    Matrix3d G;
    G << G11, G12, G13, G12, G22, G23, G13, G23, G33;
    Matrix3d D = G * EYE;
    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    double XCSQ = XSUBC*XSUBC, YCSQ = YSUBC*YSUBC, XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ + XSUBB*XSUBC + XCSQ)/6.0;
    double PY2 = YCSQ/6.0;
    double PXY2 = YSUBC*(XSUBB + 2.0*XSUBC)/12.0;
    double XBAR3 = 3.0*XBAR, YBAR3 = 3.0*YBAR, YBAR2 = 2.0*YBAR;
    MatrixXd KX = MatrixXd::Zero(6,6);
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
    HBAR(0,0)=XBSQ; HBAR(0,3)=XBSQ*XSUBB;
    HBAR(1,1)=XSUBB;
    HBAR(2,0)=-2.0*XSUBB; HBAR(2,3)=-3.0*XBSQ;
    HBAR(3,1)=XCSQ; HBAR(3,2)=XCYC; HBAR(3,3)=YCSQ;
    HBAR(3,4)=XCSQ*XSUBC; HBAR(3,5)=YCSQ*XSUBC;
    HBAR(4,0)=YCSQ*YSUBC; HBAR(4,2)=XSUBC; HBAR(4,3)=YSUBC*2.0;
    HBAR(4,4)=XCYC*2.0; HBAR(4,5)=YCSQ*3.0;
    HBAR(5,0)=-2.0*XSUBC; HBAR(5,1)=-YSUBC; HBAR(5,3)=-3.0*XCSQ; HBAR(5,4)=-YCSQ;
    bool ising;
    MatrixXd HBARinv = inverd(HBAR, ising);
    if (ising) { for(int i=0;i<3;i++)for(int j=0;j<3;j++)KsuperU[i][j]=Matrix3d::Zero(); return; }
    MatrixXd KQ = KX * HBARinv;
    MatrixXd KII = HBARinv.transpose() * KQ;
    MatrixXd S = MatrixXd::Zero(6,3);
    S(0,0)=1.0; S(0,2)=-XSUBB;
    S(1,1)=1.0;
    S(2,2)=1.0;
    S(3,0)=1.0; S(3,1)=YSUBC; S(3,2)=-XSUBC;
    S(4,1)=1.0;
    S(5,2)=1.0;
    MatrixXd KIA = KII * S;
    MatrixXd KAA = S.transpose() * KIA;
    KIA = -KIA;
    KsuperU[0][0]=KAA;
    KsuperU[0][1]=KIA.topRows(3).transpose();
    KsuperU[0][2]=KIA.bottomRows(3).transpose();
    KsuperU[1][0]=KIA.topRows(3);
    KsuperU[1][1]=KII.topLeftCorner(3,3);
    KsuperU[1][2]=KII.topRightCorner(3,3);
    KsuperU[2][0]=KII.bottomLeftCorner(3,3);
    KsuperU[2][1]=KII.bottomLeftCorner(3,3);
    KsuperU[2][2]=KII.bottomRightCorner(3,3);
}

int main() {
    E_young = 200e9; nu_val = 0.3; thickness = 0.01;
    double EYE = thickness*thickness*thickness / 12.0;
    G2X211 = 0.0; G2X212 = 0.0; G2X222 = 0.0;
    MAT_isotropic(E_young, nu_val);
    Vector3d A(0,0,0), B(1,0,0), C(1,1,0), D(0,1,0);
    Vector3d D1 = C - A, D2 = D - B, A1 = B - A;
    Vector3d KVECT = D1.cross(D2); KVECT.normalize();
    double H = A1.dot(KVECT) / 2.0;
    Vector3d IVECT = A1 - H*KVECT; IVECT.normalize();
    Vector3d JVECT = KVECT.cross(IVECT); JVECT.normalize();
    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,2)=D1.dot(IVECT); R(0,1)=A1.dot(IVECT); R(1,2)=D1.dot(JVECT);
    R(0,3)=D2.dot(IVECT)+R(0,1); R(1,3)=D2.dot(JVECT);
    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    int NPIVOT = 1;
    int JNOT = (NPIVOT-2 <= 0) ? NPIVOT+2 : NPIVOT-2;
    MatrixXd KSUMblk = MatrixXd::Zero(4,9);
    for (int J = 1; J <= 4; J++) {
        if (J == JNOT) continue;
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0], SUBSCB = M[KM+1], SUBSCC = M[KM+2];
        Vector2d V, VV;
        for (int i = 0; i < 2; i++) {
            V(i) = R(i, SUBSCB-1) - R(i, SUBSCA-1);
            VV(i) = R(i, SUBSCC-1) - R(i, SUBSCA-1);
        }
        double XSUBB = V.norm();
        double U1 = V(0)/XSUBB, U2 = V(1)/XSUBB;
        double XSUBC = U1*VV(0) + U2*VV(1);
        double YSUBC = U1*VV(1) - U2*VV(0);
        MatrixXd KsuperU[3][3];
        KTRBSC(XSUBB, XSUBC, YSUBC, EYE, 0.0, KsuperU);
        Matrix3d T;
        T << 1.0,0.0,0.0, 0.0,U1,U2, 0.0,-U2,U1;
        int P = -1;
        for (int I = 1; I <= 3; I++) if (M[KM+I-1] == NPIVOT) { P = I; break; }
        for (int I = 1; I <= 3; I++) {
            Matrix3d KU = KsuperU[P-1][I-1];
            Matrix3d PROD9 = T.transpose() * KU * T;
            int node = M[KM+I-1];
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                    KSUMblk(node-1, 3*c+r) += PROD9(r,c) / 2.0;
        }
    }
    MatrixXd E = MatrixXd::Zero(6,3);
    E(0,0)=KVECT(0); E(1,0)=KVECT(1); E(2,0)=KVECT(2);
    E(3,0)=IVECT(0); E(4,0)=IVECT(1); E(5,0)=IVECT(2);
    E(3,1)=JVECT(0); E(4,1)=JVECT(1); E(5,1)=JVECT(2);
    MatrixXd TITE = E;
    MatrixXd K12 = MatrixXd::Zero(12,12);
    for (int J = 1; J <= 4; J++) {
        MatrixXd TJTE = E;
        Matrix3d KSUMJ;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++)
                KSUMJ(r,c) = KSUMblk(J-1, 3*c+r);
        MatrixXd TEMP18 = KSUMJ * TJTE.transpose(); // 3x6
        MatrixXd KOUT = TITE * TEMP18; // 6x6
        // Place top-left 3x3 at (pivot, J)
        int prow = (NPIVOT-1)*3, pcol = (J-1)*3;
        K12.block<3,3>(prow, pcol) = KOUT.topLeftCorner(3,3);
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