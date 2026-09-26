#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct State {
    double G11, G12, G13, G22, G23, G33;
    double G2X211, G2X212, G2X222;
    double T2;
    double EYE;
    double XSUBB, XSUBC, YSUBC;
};

MatrixXd inverd(const MatrixXd& M, bool& ising) {
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = true; return MatrixXd::Zero(M.rows(), M.cols()); }
    ising = false;
    return lu.inverse();
}

MatrixXd ktrbsc(State& st) {
    double XSUBB = st.XSUBB, XSUBC = st.XSUBC, YSUBC = st.YSUBC;
    double EYE = st.EYE;
    
    Matrix3d G;
    G << st.G11, st.G12, st.G13,
         st.G12, st.G22, st.G23,
         st.G13, st.G23, st.G33;
    Matrix3d D = G * EYE;
    
    double AREA = XSUBB*YSUBC/2.0;
    double XBAR = (XSUBB+XSUBC)/3.0;
    double YBAR = YSUBC/3.0;
    double XCSQ = XSUBC*XSUBC, YCSQ = YSUBC*YSUBC, XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ+XSUBB*XSUBC+XCSQ)/6.0;
    double PY2 = YCSQ/6.0;
    double PXY2 = YSUBC*(XSUBB+2.0*XSUBC)/12.0;
    double XBAR3 = 3.0*XBAR, YBAR3 = 3.0*YBAR, YBAR2 = 2.0*YBAR;
    
    auto d = [&](int i){ return D((i-1)/3, (i-1)%3); };
    
    double a[37];
    a[1]=d(1); a[2]=d(3); a[3]=d(2);
    a[4]=d(1)*XBAR3; a[5]=d(2)*XBAR+YBAR2*d(3); a[6]=d(2)*YBAR3;
    a[7]=a[2]; a[8]=d(9); a[9]=d(6);
    a[10]=d(3)*XBAR3; a[11]=d(6)*XBAR+YBAR2*d(9); a[12]=d(6)*YBAR3;
    a[13]=a[3]; a[14]=a[9]; a[15]=d(5);
    a[16]=d(2)*XBAR3; a[17]=d(5)*XBAR+YBAR2*d(6); a[18]=d(5)*YBAR3;
    a[19]=a[4]; a[20]=a[10]; a[21]=a[16];
    a[22]=d(1)*9.0*PX2; a[23]=d(2)*3.0*PX2+6.0*PXY2*d(3); a[24]=d(2)*9.0*PXY2;
    a[25]=a[5]; a[26]=a[11]; a[27]=a[17];
    a[28]=a[23]; a[29]=d(5)*PX2+4.0*PXY2*d(6)+4.0*PY2*d(9); a[30]=d(5)*3.0*PXY2+6.0*PY2*d(6);
    a[31]=a[6]; a[32]=a[12]; a[33]=a[18];
    a[34]=a[24]; a[35]=a[30]; a[36]=d(5)*9.0*PY2;
    double TEMP = 4.0*AREA;
    for (int i=1;i<=36;i++) a[i]*=TEMP;
    
    double h[37];
    for (int i=1;i<=36;i++) h[i]=0.0;
    h[1]=XBSQ; h[4]=XBSQ*XSUBB; h[8]=XSUBB; h[13]=-2.0*XSUBB; h[16]=-3.0*XBSQ;
    h[19]=XCSQ; h[20]=XCYC; h[21]=YCSQ;
    h[22]=XCSQ*XSUBC; h[23]=YCSQ*XSUBC; h[24]=YCSQ*YSUBC;
    h[26]=XSUBC; h[27]=YSUBC*2.0;
    h[29]=XCYC*2.0; h[30]=YCSQ*3.0;
    h[31]=-2.0*XSUBC; h[32]=-YSUBC;
    h[34]=-3.0*XCSQ; h[35]=-YCSQ;
    
    if (st.T2 != 0.0 && !(st.G2X211==0.0 && st.G2X212==0.0 && st.G2X222==0.0)) {
        Matrix2d G2X2;
        G2X2 << st.G2X211*st.T2, st.G2X212*st.T2,
                st.G2X212*st.T2, st.G2X222*st.T2;
        double DETERM = G2X2(0,0)*G2X2(1,1)-G2X2(1,0)*G2X2(0,1);
        Matrix2d J2X2;
        J2X2 << G2X2(1,1)/DETERM, -G2X2(0,1)/DETERM,
                -G2X2(1,0)/DETERM, G2X2(0,0)/DETERM;
        double TEMP2 = 2.0*d(2)+4.0*d(9);
        MatrixXd HYQ(2,3);
        HYQ(0,0)=-6.0*(J2X2(0,0)*d(1)+J2X2(0,1)*d(3));
        HYQ(0,1)=-J2X2(0,0)*TEMP2-6.0*J2X2(0,1)*d(6);
        HYQ(0,2)=-6.0*(J2X2(0,0)*d(6)+J2X2(0,1)*d(5));
        HYQ(1,0)=-6.0*(J2X2(1,0)*d(1)+J2X2(1,1)*d(3));
        HYQ(1,1)=-J2X2(1,0)*TEMP2-6.0*J2X2(1,1)*d(6);
        HYQ(1,2)=-6.0*(J2X2(1,0)*d(6)+J2X2(1,1)*d(5));
        MatrixXd prod = G2X2*HYQ;
        Matrix3d KY = HYQ.transpose()*prod;
        a[22]+=KY(0,0)*AREA; a[23]+=KY(0,1)*AREA; a[24]+=KY(0,2)*AREA;
        a[28]+=KY(1,0)*AREA; a[29]+=KY(1,1)*AREA; a[30]+=KY(1,2)*AREA;
        a[34]+=KY(2,0)*AREA; a[35]+=KY(2,1)*AREA; a[36]+=KY(2,2)*AREA;
        h[4]+=XSUBB*HYQ(0,0); h[5]+=XSUBB*HYQ(0,1); h[6]+=XSUBB*HYQ(0,2);
        h[22]+=XSUBC*HYQ(0,0)+YSUBC*HYQ(1,0);
        h[23]+=XSUBC*HYQ(0,1)+YSUBC*HYQ(1,1);
        h[24]+=XSUBC*HYQ(0,2)+YSUBC*HYQ(1,2);
    }
    
    MatrixXd HBAR(6,6), KX(6,6);
    for (int i=0;i<6;i++) for (int j=0;j<6;j++) {
        HBAR(i,j)=h[i*6+j+1];
        KX(i,j)=a[i*6+j+1];
    }
    
    bool ising;
    MatrixXd Hinv = inverd(HBAR, ising);
    MatrixXd KH = KX*Hinv;
    MatrixXd KII = Hinv.transpose()*KH;
    
    MatrixXd S(6,3);
    S << 1.0,0.0,-XSUBB, 0.0,1.0,0.0, 0.0,0.0,1.0,
         1.0,YSUBC,-XSUBC, 0.0,1.0,0.0, 0.0,0.0,1.0;
    
    MatrixXd KIA = KII*S;
    MatrixXd KAA = S.transpose()*KIA;
    KIA = -KIA;
    
    MatrixXd KSU(9,9);
    KSU.block<3,3>(0,0)=KAA;
    KSU.block<3,3>(0,3)=KIA.topRows(3).transpose();
    KSU.block<3,3>(0,6)=KIA.bottomRows(3).transpose();
    KSU.block<3,3>(3,0)=KIA.topRows(3);
    KSU.block<3,3>(3,3)=KII.topLeftCorner(3,3);
    KSU.block<3,3>(3,6)=KII.topRightCorner(3,3);
    KSU.block<3,3>(6,0)=KIA.bottomRows(3);
    KSU.block<3,3>(6,3)=KII.bottomLeftCorner(3,3);
    KSU.block<3,3>(6,6)=KII.bottomRightCorner(3,3);
    return KSU;
}

int main() {
    double E=200e9, nu=0.3, t=0.01;
    double I = t*t*t/12.0;
    double G11 = E/(1.0-nu*nu);
    double G12 = nu*G11;
    double G22 = G11;
    double G33 = G11*(1.0-nu)/2.0;
    
    State st;
    st.G11=G11; st.G12=G12; st.G13=0.0;
    st.G22=G22; st.G23=0.0; st.G33=G33;
    st.G2X211=0.0; st.G2X212=0.0; st.G2X222=0.0;
    st.T2=0.0; st.EYE=I;
    
    Vector3d A(0,0,0), B(1,0,0), C(1,1,0), D(0,1,0);
    int NPIVOT=1;
    int JNOT = (NPIVOT-2 <= 0) ? NPIVOT+2 : NPIVOT-2;
    
    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    
    Vector3d D1 = C-A, D2 = D-B, A1 = B-A;
    Vector3d KVECT = D1.cross(D2); KVECT.normalize();
    double H = A1.dot(KVECT)/2.0;
    Vector3d IVECT = A1 - H*KVECT; IVECT.normalize();
    Vector3d JVECT = KVECT.cross(IVECT); JVECT.normalize();
    
    MatrixXd R = MatrixXd::Zero(2,4);
    R(0,2)=D1.dot(IVECT);
    R(0,1)=A1.dot(IVECT);
    R(1,2)=D1.dot(JVECT);
    R(0,3)=D2.dot(IVECT)+R(0,1);
    R(1,3)=D2.dot(JVECT);
    
    Matrix3d KSUM[4];
    for (int i=0;i<4;i++) KSUM[i]=Matrix3d::Zero();
    
    for (int J=1; J<=4; J++) {
        if (J==JNOT) continue;
        int KM = 3*J-3;
        int SUBSCA=M[KM+0], SUBSCB=M[KM+1], SUBSCC=M[KM+2];
        
        Vector2d V, VV;
        for (int i=0;i<2;i++) {
            V(i)=R(i,SUBSCB-1)-R(i,SUBSCA-1);
            VV(i)=R(i,SUBSCC-1)-R(i,SUBSCA-1);
        }
        double XSUBB=V.norm();
        double U1=V(0)/XSUBB, U2=V(1)/XSUBB;
        double XSUBC=U1*VV(0)+U2*VV(1);
        double YSUBC=U1*VV(1)-U2*VV(0);
        
        double SINANG=0.0, COSANG=1.0;
        double SINTH=SINANG*U1-COSANG*U2;
        if (fabs(SINTH)<1e-6) SINTH=0.0;
        
        st.XSUBB=XSUBB; st.XSUBC=XSUBC; st.YSUBC=YSUBC;
        
        MatrixXd KSU = ktrbsc(st);
        
        Matrix3d T;
        T << 1.0,0.0,0.0, 0.0,U1,U2, 0.0,-U2,U1;
        
        int pivot_pos = -1;
        for (int i=1;i<=3;i++) {
            if (M[KM+i-1]==NPIVOT) { pivot_pos=i; break; }
        }
        int block_row = pivot_pos-1;
        
        for (int i=1;i<=3;i++) {
            int block_col = i-1;
            Matrix3d Ablock = KSU.block<3,3>(block_row*3, block_col*3);
            Matrix3d TEMP9 = T.transpose()*Ablock;
            Matrix3d PROD9 = TEMP9*T;
            
            int NPOINT2 = KM+i;
            int node = M[NPOINT2-1];
            KSUM[node-1] += PROD9/2.0;
        }
    }
    
    MatrixXd Emat = MatrixXd::Zero(6,3);
    Emat(0,0)=KVECT(0); Emat(1,0)=KVECT(1); Emat(2,0)=KVECT(2);
    Emat(3,0)=IVECT(0); Emat(4,0)=IVECT(1); Emat(5,0)=IVECT(2);
    Emat(3,1)=JVECT(0); Emat(4,1)=JVECT(1); Emat(5,1)=JVECT(2);
    
    MatrixXd TITE = Emat;
    
    MatrixXd K12 = MatrixXd::Zero(12,12);
    int pivot_idx = NPIVOT-1;
    int r0 = 3*pivot_idx;
    
    for (int J=1; J<=4; J++) {
        MatrixXd TJTE = Emat;
        Matrix3d KsumJ = KSUM[J-1];
        MatrixXd TEMP18 = KsumJ * TJTE.transpose();
        MatrixXd KOUT = TITE * TEMP18;
        
        int node = J;
        int c0 = 3*(node-1);
        
        K12.block<3,3>(r0, r0) += KOUT.block<3,3>(0,0);
        K12.block<3,3>(r0, c0) += KOUT.block<3,3>(0,3);
        K12.block<3,3>(c0, r0) += KOUT.block<3,3>(3,0);
        K12.block<3,3>(c0, c0) += KOUT.block<3,3>(3,3);
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