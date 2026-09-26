#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct Material {
    double E, nu, t;
    double G11, G12, G13, G22, G23, G33;
    double I;
};

void computeG(Material& mat) {
    double E = mat.E, nu = mat.nu;
    double denom = 1.0 - nu*nu;
    mat.G11 = E / denom;
    mat.G12 = nu * E / denom;
    mat.G13 = 0.0;
    mat.G22 = E / denom;
    mat.G23 = 0.0;
    mat.G33 = E / (2.0*(1.0+nu));
}

MatrixXd KTRBSC(double XSUBB, double XSUBC, double YSUBC, const Material& mat) {
    Matrix3d G;
    G << mat.G11, mat.G12, mat.G13,
         mat.G12, mat.G22, mat.G23,
         mat.G13, mat.G23, mat.G33;
    Matrix3d D = G * mat.I;
    
    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    
    double XCSQ = XSUBC*XSUBC;
    double YCSQ = YSUBC*YSUBC;
    double XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ + XSUBB*XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC*(XSUBB + 2.0*XSUBC) / 12.0;
    double XBAR3 = 3.0*XBAR;
    double YBAR3 = 3.0*YBAR;
    double YBAR2 = 2.0*YBAR;
    
    double d1 = D(0,0), d2 = D(0,1), d3 = D(0,2);
    double d5 = D(1,1), d6 = D(1,2), d9 = D(2,2);
    
    MatrixXd KX(6,6);
    KX(0,0) = d1;
    KX(0,1) = d3;
    KX(0,2) = d2;
    KX(0,3) = d1*XBAR3;
    KX(0,4) = d2*XBAR + YBAR2*d3;
    KX(0,5) = d2*YBAR3;
    KX(1,0) = d3;
    KX(1,1) = d9;
    KX(1,2) = d6;
    KX(1,3) = d3*XBAR3;
    KX(1,4) = d6*XBAR + YBAR2*d9;
    KX(1,5) = d6*YBAR3;
    KX(2,0) = d2;
    KX(2,1) = d6;
    KX(2,2) = d5;
    KX(2,3) = d2*XBAR3;
    KX(2,4) = d5*XBAR + YBAR2*d6;
    KX(2,5) = d5*YBAR3;
    KX(3,0) = d1*XBAR3;
    KX(3,1) = d3*XBAR3;
    KX(3,2) = d2*XBAR3;
    KX(3,3) = d1*9.0*PX2;
    KX(3,4) = d2*3.0*PX2 + 6.0*PXY2*d3;
    KX(3,5) = d2*9.0*PXY2;
    KX(4,0) = d2*XBAR + YBAR2*d3;
    KX(4,1) = d6*XBAR + YBAR2*d9;
    KX(4,2) = d5*XBAR + YBAR2*d6;
    KX(4,3) = d2*3.0*PX2 + 6.0*PXY2*d3;
    KX(4,4) = d5*PX2 + 4.0*PXY2*d6 + 4.0*PY2*d9;
    KX(4,5) = d5*3.0*PXY2 + 6.0*PY2*d6;
    KX(5,0) = d2*YBAR3;
    KX(5,1) = d6*YBAR3;
    KX(5,2) = d5*YBAR3;
    KX(5,3) = d2*9.0*PXY2;
    KX(5,4) = d5*3.0*PXY2 + 6.0*PY2*d6;
    KX(5,5) = d5*9.0*PY2;
    
    double TEMP = 4.0*AREA;
    KX = KX * TEMP;
    
    MatrixXd HBAR = MatrixXd::Zero(6,6);
    HBAR(0,0) = XBSQ;
    HBAR(3,3) = XBSQ*XSUBB;
    HBAR(1,1) = XSUBB;
    HBAR(2,0) = -2.0*XSUBB;
    HBAR(2,3) = -3.0*XBSQ;
    HBAR(3,0) = XCSQ;
    HBAR(3,1) = XCYC;
    HBAR(3,2) = YCSQ;
    HBAR(3,3) = XCSQ*XSUBC;
    HBAR(3,4) = YCSQ*XSUBC;
    HBAR(3,5) = YCSQ*YSUBC;
    HBAR(4,1) = XSUBC;
    HBAR(4,2) = YSUBC*2.0;
    HBAR(4,4) = XCYC*2.0;
    HBAR(4,5) = YCSQ*3.0;
    HBAR(5,0) = -2.0*XSUBC;
    HBAR(5,1) = -YSUBC;
    HBAR(5,3) = -3.0*XCSQ;
    HBAR(5,4) = -YCSQ;
    
    MatrixXd Hinv = HBAR.inverse();
    MatrixXd KQ = KX * Hinv;
    MatrixXd KII = Hinv.transpose() * KQ;
    
    MatrixXd S(6,3);
    S << 1.0, 0.0, -XSUBB,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         1.0, YSUBC, -XSUBC,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;
    
    MatrixXd KIA = KII * S;
    MatrixXd KAA = S.transpose() * KIA;
    KIA = -KIA;
    
    Matrix3d KIA_top = KIA.block<3,3>(0,0);
    Matrix3d KIA_bot = KIA.block<3,3>(3,0);
    
    MatrixXd KU = MatrixXd::Zero(9,9);
    KU.block<3,3>(0,0) = KAA;
    KU.block<3,3>(0,1) = KIA_top.transpose();
    KU.block<3,3>(0,2) = KIA_bot.transpose();
    KU.block<3,3>(1,0) = KIA_top;
    KU.block<3,3>(1,1) = KII.block<3,3>(0,0);
    KU.block<3,3>(1,2) = KII.block<3,3>(0,3);
    KU.block<3,3>(2,0) = KIA_bot;
    KU.block<3,3>(2,1) = KII.block<3,3>(3,0);
    KU.block<3,3>(2,2) = KII.block<3,3>(3,3);
    
    return KU;
}

int main() {
    Material mat;
    mat.E = 200e9;
    mat.nu = 0.3;
    mat.t = 0.01;
    mat.I = mat.t*mat.t*mat.t / 12.0;
    computeG(mat);
    
    Vector3d A(0,0,0), B(1,0,0), C(1,1,0), D(0,1,0);
    
    int NPIVOT = 1;
    int JNOT = NPIVOT + 2;
    
    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    
    Vector3d VQ[5];
    VQ[1] = A; VQ[2] = B; VQ[3] = C; VQ[4] = D;
    
    Vector3d D1 = VQ[3] - VQ[1];
    Vector3d D2 = VQ[4] - VQ[2];
    Vector3d A1 = VQ[2] - VQ[1];
    
    Vector3d KVECT = D1.cross(D2);
    KVECT.normalize();
    
    double H = A1.dot(KVECT) / 2.0;
    
    Vector3d IVECT = A1 - H*KVECT;
    IVECT.normalize();
    
    Vector3d JVECT = KVECT.cross(IVECT);
    JVECT.normalize();
    
    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,2) = D1.dot(IVECT);
    R(0,1) = A1.dot(IVECT);
    R(1,2) = D1.dot(JVECT);
    R(0,3) = D2.dot(IVECT) + R(0,1);
    R(1,3) = D2.dot(JVECT);
    
    Matrix3d KSUM[5];
    for (int i = 1; i <= 4; i++) KSUM[i] = Matrix3d::Zero();
    
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
        
        MatrixXd KU = KTRBSC(XSUBB, XSUBC, YSUBC, mat);
        
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
        
        for (int i = 1; i <= 3; i++) {
            int rowBlock = NBEGIN / 9;
            int colBlock = i - 1;
            Matrix3d Ablock = KU.block<3,3>(rowBlock, colBlock*3);
            Matrix3d PROD9 = T.transpose() * Ablock * T;
            int NPOINT = KM + i;
            int nodeIdx = M[NPOINT-1];
            KSUM[nodeIdx] += PROD9 / 2.0;
        }
    }
    
    MatrixXd E = MatrixXd::Zero(6,3);
    E(0,0) = KVECT(0); E(1,0) = KVECT(1); E(2,0) = KVECT(2);
    E(3,0) = IVECT(0); E(4,0) = IVECT(1); E(5,0) = IVECT(2);
    E(3,1) = JVECT(0); E(4,1) = JVECT(1); E(5,1) = JVECT(2);
    
    MatrixXd TITE = E;
    
    MatrixXd KOUT = MatrixXd::Zero(12,12);
    int pivotBase = (NPIVOT-1)*3;
    for (int J = 1; J <= 4; J++) {
        MatrixXd TJTE = E;
        Matrix3d Kblock = KSUM[J];
        MatrixXd TEMP18 = Kblock * TJTE.transpose();
        MatrixXd KOUTJ = TITE * TEMP18;
        int jBase = (J-1)*3;
        KOUT.block<3,3>(pivotBase, pivotBase) += KOUTJ.block<3,3>(0,0);
        KOUT.block<3,3>(pivotBase, jBase)    += KOUTJ.block<3,3>(0,3);
        KOUT.block<3,3>(jBase, pivotBase)    += KOUTJ.block<3,3>(3,0);
        KOUT.block<3,3>(jBase, jBase)        += KOUTJ.block<3,3>(3,3);
    }
    
    cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        cout << "[";
        for (int j = 0; j < 12; j++) {
            cout << scientific << setprecision(6) << KOUT(i,j);
            if (j < 11) cout << ",";
        }
        cout << "]";
        if (i < 11) cout << ",";
    }
    cout << "]}" << endl;
    
    return 0;
}