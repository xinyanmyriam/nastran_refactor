#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct MatData {
    double G11, G12, G13, G22, G23, G33;
    double G2X211, G2X212, G2X222;
};

MatData getMat(double Emod, double nu, double t) {
    MatData m;
    double I = t*t*t/12.0;
    double D0 = Emod*I/(1.0-nu*nu);
    m.G11 = D0; m.G12 = D0*nu; m.G13 = 0.0;
    m.G22 = D0; m.G23 = 0.0; m.G33 = D0*(1.0-nu)/2.0;
    m.G2X211 = 0.0; m.G2X212 = 0.0; m.G2X222 = 0.0;
    return m;
}

MatrixXd inverd(const MatrixXd& H, bool& ising) {
    int n = H.rows();
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) { ising = true; return MatrixXd::Zero(n,n); }
    ising = false;
    return lu.inverse();
}

MatrixXd ktrbsc(double XSUBB, double XSUBC, double YSUBC,
                const MatData& mat, double EYE, double T2, bool& nogo) {
    nogo = false;
    Matrix3d G;
    G << mat.G11, mat.G12, mat.G13,
         mat.G12, mat.G22, mat.G23,
         mat.G13, mat.G23, mat.G33;
    Matrix3d D = G * EYE;
    
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
    
    vector<double> a(200, 0.0);
    auto A_ = [&](int i) -> double& { return a[i]; };
    
    double d1=D(0,0), d2=D(0,1), d3=D(0,2);
    double d5=D(1,1), d6=D(1,2), d9=D(2,2);
    
    A_(1)=d1; A_(2)=d3; A_(3)=d2;
    A_(4)=d1*XBAR3; A_(5)=d2*XBAR+YBAR2*d3; A_(6)=d2*YBAR3;
    A_(7)=A_(2); A_(8)=d9; A_(9)=d6;
    A_(10)=d3*XBAR3; A_(11)=d6*XBAR+YBAR2*d9; A_(12)=d6*YBAR3;
    A_(13)=A_(3); A_(14)=A_(9); A_(15)=d5;
    A_(16)=d2*XBAR3; A_(17)=d5*XBAR+YBAR2*d6; A_(18)=d5*YBAR3;
    A_(19)=A_(4); A_(20)=A_(10); A_(21)=A_(16);
    A_(22)=d1*9.0*PX2; A_(23)=d2*3.0*PX2+6.0*PXY2*d3; A_(24)=d2*9.0*PXY2;
    A_(25)=A_(5); A_(26)=A_(11); A_(27)=A_(17);
    A_(28)=A_(23); A_(29)=d5*PX2+4.0*PXY2*d6+4.0*PY2*d9; A_(30)=d5*3.0*PXY2+6.0*PY2*d6;
    A_(31)=A_(6); A_(32)=A_(12); A_(33)=A_(18);
    A_(34)=A_(24); A_(35)=A_(30); A_(36)=d5*9.0*PY2;
    double TEMP = 4.0*AREA;
    for (int i=1;i<=36;i++) A_(i) *= TEMP;
    
    for (int i=37;i<=72;i++) A_(i)=0.0;
    A_(37)=XBSQ; A_(40)=XBSQ*XSUBB; A_(44)=XSUBB;
    A_(49)=-2.0*XSUBB; A_(52)=-3.0*XBSQ;
    A_(55)=XCSQ; A_(56)=XCYC; A_(57)=YCSQ;
    A_(58)=XCSQ*XSUBC; A_(59)=YCSQ*XSUBC; A_(60)=YCSQ*YSUBC;
    A_(62)=XSUBC; A_(63)=YSUBC*2.0;
    A_(65)=XCYC*2.0; A_(66)=YCSQ*3.0;
    A_(67)=-2.0*XSUBC; A_(68)=-YSUBC;
    A_(70)=-3.0*XCSQ; A_(71)=-YCSQ;
    
    MatrixXd H(6,6);
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) H(i,j)=A_(37+i*6+j);
    bool ising;
    MatrixXd Hinv = inverd(H, ising);
    if (ising) { nogo=true; return MatrixXd::Zero(9,9); }
    
    MatrixXd KX(6,6);
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) KX(i,j)=A_(1+i*6+j);
    MatrixXd KQ = KX * Hinv;
    MatrixXd KII = Hinv.transpose() * KQ;
    
    MatrixXd S(6,3);
    S << 1.0,0.0,-XSUBB, 0.0,1.0,0.0, 0.0,0.0,1.0,
         1.0,YSUBC,-XSUBC, 0.0,1.0,0.0, 0.0,0.0,1.0;
    MatrixXd KIA = KII * S;
    MatrixXd KAA = S.transpose() * KIA;
    KIA = -KIA;
    
    MatrixXd KU = MatrixXd::Zero(9,9);
    KU.block<3,3>(0,0) = KAA;
    KU.block<3,3>(0,3) = KIA.block<3,3>(0,0);
    KU.block<3,3>(0,6) = KIA.block<3,3>(3,0);
    KU.block<3,3>(3,0) = KIA.block<3,3>(0,0).transpose();
    KU.block<3,3>(6,0) = KIA.block<3,3>(3,0).transpose();
    KU.block<3,3>(3,3) = KII.block<3,3>(0,0);
    KU.block<3,3>(3,6) = KII.block<3,3>(0,3);
    KU.block<3,3>(6,3) = KII.block<3,3>(3,0);
    KU.block<3,3>(6,6) = KII.block<3,3>(3,3);
    return KU;
}

int main() {
    double Emod = 200e9, nu = 0.3, t = 0.01;
    double I = t*t*t/12.0;
    Vector3d A(0,0,0), B(1,0,0), C(1,1,0), D(0,1,0);
    MatData mat = getMat(Emod, nu, t);
    int NPIVOT = 1;
    int JNOT = (NPIVOT-2 <= 0) ? NPIVOT+2 : NPIVOT-2;
    Vector3d VQ1=A, VQ2=B, VQ3=C, VQ4=D;
    Vector3d D1 = VQ3-VQ1, D2 = VQ4-VQ2, A1 = VQ2-VQ1;
    Vector3d KV = D1.cross(D2);
    double temp = KV.norm(); KV /= temp;
    temp = A1.dot(KV)/2.0;
    Vector3d IV = A1 - temp*KV; IV /= IV.norm();
    Vector3d JV = KV.cross(IV); JV /= JV.norm();
    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,2)=D1.dot(IV); R(0,1)=A1.dot(IV); R(1,2)=D1.dot(JV);
    R(0,3)=D2.dot(IV)+R(0,1); R(1,3)=D2.dot(JV);
    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    
    vector<Matrix3d> KSUM(5, Matrix3d::Zero());
    
    for (int J=1; J<=4; J++) {
        if (J == JNOT) continue;
        int KM = 3*J-3;
        int SUBSCA=M[KM+0], SUBSCB=M[KM+1], SUBSCC=M[KM+2];
        Vector2d V, VV;
        for (int i=0;i<2;i++){ V(i)=R(i,SUBSCB-1)-R(i,SUBSCA-1); VV(i)=R(i,SUBSCC-1)-R(i,SUBSCA-1);}
        double XSUBB=V.norm();
        double U1=V(0)/XSUBB, U2=V(1)/XSUBB;
        double XSUBC=U1*VV(0)+U2*VV(1);
        double YSUBC=U1*VV(1)-U2*VV(0);
        bool nogo;
        MatrixXd KU = ktrbsc(XSUBB,XSUBC,YSUBC,mat,I,0.0,nogo);
        if (nogo){cerr<<"err"<<endl;return 1;}
        Matrix3d T;
        T << 1.0,0.0,0.0, 0.0,U1,U2, 0.0,-U2,U1;
        for (int i=1;i<=3;i++) {
            Matrix3d Ablock;
            for (int r=0;r<3;r++) for (int c=0;c<3;c++) Ablock(r,c)=KU((i-1)*3+r,c);
            Matrix3d PROD9 = T.transpose()*Ablock*T;
            int NPOINT=KM+i;
            int node = M[NPOINT-1];
            KSUM[node] += PROD9/2.0;
        }
    }
    
    MatrixXd Emat = MatrixXd::Zero(6,3);
    Emat(0,0)=KV(0); Emat(1,0)=KV(1); Emat(2,0)=KV(2);
    Emat(3,0)=IV(0); Emat(4,0)=IV(1); Emat(5,0)=IV(2);
    Emat(3,1)=JV(0); Emat(4,1)=JV(1); Emat(5,1)=JV(2);
    MatrixXd TITE = Emat;
    
    MatrixXd K12 = MatrixXd::Zero(12,12);
    int pivotDof = (NPIVOT-1)*3;
    for (int J=1; J<=4; J++) {
        Matrix3d KsumJ = KSUM[J];
        MatrixXd TJTE = Emat;
        MatrixXd TEMP18 = KsumJ * TJTE.transpose(); // 3x6
        MatrixXd KOUT = TITE * TEMP18; // 6x6
        int jDof = (J-1)*3;
        // KOUT rows/cols: [0:3]=pivot, [3:6]=node J
        K12.block<3,3>(pivotDof, pivotDof) += KOUT.block<3,3>(0,0);
        K12.block<3,3>(pivotDof, jDof)     += KOUT.block<3,3>(0,3);
        K12.block<3,3>(jDof, pivotDof)     += KOUT.block<3,3>(3,0);
        K12.block<3,3>(jDof, jDof)         += KOUT.block<3,3>(3,3);
    }
    
    cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<12;i++) {
        cout << "[";
        for (int j=0;j<12;j++) {
            cout << scientific << setprecision(6) << K12(i,j);
            if (j<11) cout << ",";
        }
        cout << "]";
        if (i<11) cout << ",";
    }
    cout << "]}" << endl;
    return 0;
}