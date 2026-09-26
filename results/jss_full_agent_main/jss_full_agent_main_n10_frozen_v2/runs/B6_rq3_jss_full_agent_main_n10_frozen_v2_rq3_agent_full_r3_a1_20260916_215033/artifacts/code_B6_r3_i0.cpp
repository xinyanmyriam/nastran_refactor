#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct Material {
    double E, nu, t, I;
    double G11, G12, G13, G22, G23, G33;
};

void computeG(Material& mat) {
    double factor = mat.E / (1.0 - mat.nu*mat.nu);
    mat.G11 = factor;
    mat.G12 = factor * mat.nu;
    mat.G13 = 0.0;
    mat.G22 = factor;
    mat.G23 = 0.0;
    mat.G33 = factor * (1.0 - mat.nu) / 2.0;
}

MatrixXd KTRBSC(const Material& mat,
                double xsubb, double xsubc, double ysubc) {
    Matrix3d G;
    G << mat.G11, mat.G12, mat.G13,
         mat.G12, mat.G22, mat.G23,
         mat.G13, mat.G23, mat.G33;
    Matrix3d D = G * mat.I;

    double d1 = D(0,0), d2 = D(0,1), d3 = D(0,2);
    double d5 = D(1,1), d6 = D(1,2), d9 = D(2,2);

    double area = xsubb * ysubc / 2.0;
    double xbar = (xsubb + xsubc) / 3.0;
    double ybar = ysubc / 3.0;

    double xcsq = xsubc*xsubc;
    double ycsq = ysubc*ysubc;
    double xbsq = xsubb*xsubb;
    double xcyc = xsubc*ysubc;
    double px2 = (xbsq + xsubb*xsubc + xcsq) / 6.0;
    double py2 = ycsq / 6.0;
    double pxy2 = ysubc*(xsubb + 2.0*xsubc) / 12.0;
    double xbar3 = 3.0*xbar;
    double ybar3 = 3.0*ybar;
    double ybar2 = 2.0*ybar;

    MatrixXd KX = MatrixXd::Zero(6, 6);
    KX(0,0) = d1;
    KX(0,1) = d3;
    KX(0,2) = d2;
    KX(0,3) = d1*xbar3;
    KX(0,4) = d2*xbar + ybar2*d3;
    KX(0,5) = d2*ybar3;
    KX(1,0) = d3;
    KX(1,1) = d9;
    KX(1,2) = d6;
    KX(1,3) = d3*xbar3;
    KX(1,4) = d6*xbar + ybar2*d9;
    KX(1,5) = d6*ybar3;
    KX(2,0) = d2;
    KX(2,1) = d6;
    KX(2,2) = d5;
    KX(2,3) = d2*xbar3;
    KX(2,4) = d5*xbar + ybar2*d6;
    KX(2,5) = d5*ybar3;
    KX(3,0) = d1*xbar3;
    KX(3,1) = d3*xbar3;
    KX(3,2) = d2*xbar3;
    KX(3,3) = d1*9.0*px2;
    KX(3,4) = d2*3.0*px2 + 6.0*pxy2*d3;
    KX(3,5) = d2*9.0*pxy2;
    KX(4,0) = d2*xbar + ybar2*d3;
    KX(4,1) = d6*xbar + ybar2*d9;
    KX(4,2) = d5*xbar + ybar2*d6;
    KX(4,3) = d2*3.0*px2 + 6.0*pxy2*d3;
    KX(4,4) = d5*px2 + 4.0*pxy2*d6 + 4.0*py2*d9;
    KX(4,5) = d5*3.0*pxy2 + 6.0*py2*d6;
    KX(5,0) = d2*ybar3;
    KX(5,1) = d6*ybar3;
    KX(5,2) = d5*ybar3;
    KX(5,3) = d2*9.0*pxy2;
    KX(5,4) = d5*3.0*pxy2 + 6.0*py2*d6;
    KX(5,5) = d5*9.0*py2;
    KX *= 4.0*area;

    MatrixXd HB = MatrixXd::Zero(6, 6);
    HB(0,0) = xbsq;
    HB(0,3) = xbsq*xsubb;
    HB(1,1) = xsubb;
    HB(2,0) = -2.0*xsubb;
    HB(2,3) = -3.0*xbsq;
    HB(3,0) = xcsq;
    HB(3,1) = xcyc;
    HB(3,2) = ycsq;
    HB(3,3) = xcsq*xsubc;
    HB(3,4) = ycsq*xsubc;
    HB(3,5) = ycsq*ysubc;
    HB(4,1) = xsubc;
    HB(4,2) = ysubc*2.0;
    HB(4,4) = xcyc*2.0;
    HB(4,5) = ycsq*3.0;
    HB(5,0) = -2.0*xsubc;
    HB(5,1) = -ysubc;
    HB(5,3) = -3.0*xcsq;
    HB(5,4) = -ycsq;

    MatrixXd HBinv = HB.inverse();
    MatrixXd KQ = KX * HBinv;
    MatrixXd KII = HBinv.transpose() * KQ;

    MatrixXd S = MatrixXd::Zero(6, 3);
    S(0,0) = 1.0; S(0,2) = -xsubb;
    S(1,1) = 1.0;
    S(2,2) = 1.0;
    S(3,0) = 1.0; S(3,1) = ysubc; S(3,2) = -xsubc;
    S(4,1) = 1.0;
    S(5,2) = 1.0;

    MatrixXd KIA = KII * S;
    MatrixXd KAA = S.transpose() * KIA;
    KIA = -KIA;

    MatrixXd KU = MatrixXd::Zero(9, 9);
    KU.block<3,3>(0,0) = KAA;
    KU.block<3,3>(0,3) = KIA.transpose();
    KU.block<3,3>(0,6) = KIA.transpose();
    KU.block<3,3>(3,0) = KIA;
    KU.block<3,3>(6,0) = KIA;
    KU.block<3,3>(3,3) = KII.block<3,3>(0,0);
    KU.block<3,3>(3,6) = KII.block<3,3>(0,3);
    KU.block<3,3>(6,3) = KII.block<3,3>(3,0);
    KU.block<3,3>(6,6) = KII.block<3,3>(3,3);
    return KU;
}

int main() {
    Material mat;
    mat.E = 200e9;
    mat.nu = 0.3;
    mat.t = 0.01;
    mat.I = mat.t*mat.t*mat.t / 12.0;
    computeG(mat);

    double VQ1[3] = {0,0,0};
    double VQ2[3] = {1,0,0};
    double VQ3[3] = {1,1,0};
    double VQ4[3] = {0,1,0};

    double D1[3], D2[3], A1[3];
    for (int i=0;i<3;i++) {
        D1[i] = VQ3[i] - VQ1[i];
        D2[i] = VQ4[i] - VQ2[i];
        A1[i] = VQ2[i] - VQ1[i];
    }

    double KV[3];
    KV[0] = D1[1]*D2[2] - D2[1]*D1[2];
    KV[1] = D1[2]*D2[0] - D2[2]*D1[0];
    KV[2] = D1[0]*D2[1] - D2[0]*D1[1];
    double temp = sqrt(KV[0]*KV[0]+KV[1]*KV[1]+KV[2]*KV[2]);
    for (int i=0;i<3;i++) KV[i] /= temp;

    temp = (A1[0]*KV[0]+A1[1]*KV[1]+A1[2]*KV[2])/2.0;

    double IV[3];
    for (int i=0;i<3;i++) IV[i] = A1[i] - temp*KV[i];
    temp = sqrt(IV[0]*IV[0]+IV[1]*IV[1]+IV[2]*IV[2]);
    for (int i=0;i<3;i++) IV[i] /= temp;

    double JV[3];
    JV[0] = KV[1]*IV[2] - IV[1]*KV[2];
    JV[1] = KV[2]*IV[0] - IV[2]*KV[0];
    JV[2] = KV[0]*IV[1] - IV[1]*KV[0];
    temp = sqrt(JV[0]*JV[0]+JV[1]*JV[1]+JV[2]*JV[2]);
    for (int i=0;i<3;i++) JV[i] /= temp;

    double R[2][4];
    R[0][0] = 0; R[1][0] = 0;
    R[0][2] = D1[0]*IV[0]+D1[1]*IV[1]+D1[2]*IV[2];
    R[1][2] = D1[0]*JV[0]+D1[1]*JV[1]+D1[2]*JV[2];
    R[0][1] = A1[0]*IV[0]+A1[1]*IV[1]+A1[2]*IV[2];
    R[0][3] = D2[0]*IV[0]+D2[1]*IV[1]+D2[2]*IV[2] + R[0][1];
    R[1][3] = D2[0]*JV[0]+D2[1]*JV[1]+D2[2]*JV[2];

    int NPIVOT = 1;
    int JNOT = (NPIVOT - 2 <= 0) ? NPIVOT + 2 : NPIVOT - 2;

    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};

    MatrixXd KSUM = MatrixXd::Zero(12, 12);

    for (int J = 1; J <= 4; J++) {
        if (J == JNOT) continue;
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0];
        int SUBSCB = M[KM+1];
        int SUBSCC = M[KM+2];

        double V[2], VV[2];
        for (int i=0;i<2;i++) {
            V[i] = R[i][SUBSCB-1] - R[i][SUBSCA-1];
            VV[i] = R[i][SUBSCC-1] - R[i][SUBSCA-1];
        }
        double XSUBB = sqrt(V[0]*V[0]+V[1]*V[1]);
        double U1 = V[0]/XSUBB;
        double U2 = V[1]/XSUBB;
        double XSUBC = U1*VV[0] + U2*VV[1];
        double YSUBC = U1*VV[1] - U2*VV[0];

        MatrixXd KU = KTRBSC(mat, XSUBB, XSUBC, YSUBC);

        Matrix3d T;
        T << 1, 0, 0,
             0, U1, U2,
             0, -U2, U1;

        int NBEGIN = 0;
        for (int i = 1; i <= 3; i++) {
            int NPOINT = KM + i;
            if (M[NPOINT-1] == NPIVOT) {
                NBEGIN = 27*i - 27;
                break;
            }
        }

        for (int i = 1; i <= 3; i++) {
            int rowBlock = NBEGIN / 9;
            Matrix3d Ablock = KU.block<3,3>(rowBlock*3, (i-1)*3).eval();
            Matrix3d PROD9 = T.transpose() * Ablock * T;
            int NPOINT = KM + i;
            int targetRow = M[NPOINT-1] - 1;
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++)
                    KSUM(targetRow*3+r, (i-1)*3+c) += PROD9(r,c) / 2.0;
        }
    }

    MatrixXd E = MatrixXd::Zero(6, 3);
    E(0,0) = KV[0]; E(1,0) = KV[1]; E(2,0) = KV[2];
    E(3,0) = IV[0]; E(4,0) = IV[1]; E(5,0) = IV[2];
    E(3,1) = JV[0]; E(4,1) = JV[1]; E(5,1) = JV[2];

    MatrixXd TITE = E;

    MatrixXd K12 = MatrixXd::Zero(12, 12);
    for (int J = 1; J <= 4; J++) {
        MatrixXd TJTE = E;
        Matrix3d KSUMblock = KSUM.block<3,3>((J-1)*3, (J-1)*3).eval();
        MatrixXd TEMP18 = KSUMblock * TJTE.transpose();
        MatrixXd KOUT = TITE * TEMP18;

        int pr = (NPIVOT-1)*3;
        int jc = (J-1)*3;
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) {
                K12(pr+r, jc+c) += KOUT(r, c);
                K12(jc+r, pr+c) += KOUT(3+r, c);
                K12(jc+r, jc+c) += KOUT(3+r, 3+c);
            }
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