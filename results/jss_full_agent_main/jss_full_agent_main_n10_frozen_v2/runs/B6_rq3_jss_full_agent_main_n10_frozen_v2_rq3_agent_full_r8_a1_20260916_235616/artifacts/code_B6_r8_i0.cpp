#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace Eigen;
using namespace std;

struct State {
    double G11, G12, G13, G22, G23, G33;
    double E, nu, t;
    double I;
    double DEGRA = M_PI/180.0;
    double SINTH, COSTH;
    bool NOGO = false;
};

static State S;

void MAT_iso() {
    double E = S.E, nu = S.nu;
    double coef = E/(1.0 - nu*nu);
    S.G11 = coef * 1.0;
    S.G12 = coef * nu;
    S.G13 = 0.0;
    S.G22 = coef * 1.0;
    S.G23 = 0.0;
    S.G33 = coef * (1.0 - nu)/2.0;
}

MatrixXd inverd(const MatrixXd& H, bool& ising) {
    FullPivLU<MatrixXd> lu(H);
    if (!lu.isInvertible()) { ising = true; return MatrixXd::Zero(H.rows(), H.cols()); }
    ising = false;
    return lu.inverse();
}

MatrixXd KTRBSC(double x1,double y1,double z1,
                double x2,double y2,double z2,
                double x3,double y3,double z3,
                double angle, double eye, double t2) {
    MatrixXd E = MatrixXd::Zero(3,6);
    auto eget = [&](int idx)->double& {
        int c = (idx-1)/3; int r = (idx-1)%3; return E(r,c);
    };
    
    eget(11) = x2 - x1;
    eget(14) = y2 - y1;
    eget(17) = z2 - z1;
    double XSUBB = sqrt(eget(11)*eget(11)+eget(14)*eget(14)+eget(17)*eget(17));
    if (XSUBB <= 1e-6) { S.NOGO = true; return MatrixXd::Zero(9,9); }
    eget(11)/=XSUBB; eget(14)/=XSUBB; eget(17)/=XSUBB;
    
    eget(2) = x3 - x1;
    eget(5) = y3 - y1;
    eget(8) = z3 - z1;
    double XSUBC = eget(11)*eget(2)+eget(14)*eget(5)+eget(17)*eget(8);
    
    eget(1) = eget(14)*eget(8) - eget(5)*eget(17);
    eget(4) = eget(2)*eget(17) - eget(11)*eget(8);
    eget(7) = eget(11)*eget(5) - eget(2)*eget(14);
    double YSUBC = sqrt(eget(1)*eget(1)+eget(4)*eget(4)+eget(7)*eget(7));
    if (YSUBC <= 1e-6) { S.NOGO = true; return MatrixXd::Zero(9,9); }
    eget(1)/=YSUBC; eget(4)/=YSUBC; eget(7)/=YSUBC;
    
    eget(12) = eget(4)*eget(17) - eget(14)*eget(7);
    eget(15) = eget(11)*eget(7) - eget(1)*eget(17);
    eget(18) = eget(1)*eget(14) - eget(11)*eget(4);
    double TEMP = sqrt(eget(12)*eget(12)+eget(15)*eget(15)+eget(18)*eget(18));
    eget(12)/=TEMP; eget(15)/=TEMP; eget(18)/=TEMP;
    eget(2)=0; eget(3)=0; eget(5)=0; eget(6)=0; eget(8)=0; eget(9)=0;
    eget(10)=0; eget(13)=0; eget(16)=0;
    
    double THETA = angle*S.DEGRA;
    S.SINTH = sin(THETA);
    S.COSTH = cos(THETA);
    if (fabs(S.SINTH) < 1e-6) S.SINTH = 0.0;
    
    MAT_iso();
    Matrix3d G;
    G << S.G11,S.G12,S.G13,
         S.G12,S.G22,S.G23,
         S.G13,S.G23,S.G33;
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
    
    vector<double> A(145, 0.0);
    auto d = [&](int i)->double { return D((i-1)%3, (i-1)/3); };
    
    A[1] = d(1);
    A[2] = d(3);
    A[3] = d(2);
    A[4] = d(1)*XBAR3;
    A[5] = d(2)*XBAR + YBAR2*d(3);
    A[6] = d(2)*YBAR3;
    A[7] = A[2];
    A[8] = d(9);
    A[9] = d(6);
    A[10] = d(3)*XBAR3;
    A[11] = d(6)*XBAR + YBAR2*d(9);
    A[12] = d(6)*YBAR3;
    A[13] = A[3];
    A[14] = A[9];
    A[15] = d(5);
    A[16] = d(2)*XBAR3;
    A[17] = d(5)*XBAR + YBAR2*d(6);
    A[18] = d(5)*YBAR3;
    A[19] = A[4];
    A[20] = A[10];
    A[21] = A[16];
    A[22] = d(1)*9.0*PX2;
    A[23] = d(2)*3.0*PX2 + 6.0*PXY2*d(3);
    A[24] = d(2)*9.0*PXY2;
    A[25] = A[5];
    A[26] = A[11];
    A[27] = A[17];
    A[28] = A[23];
    A[29] = d(5)*PX2 + 4.0*PXY2*d(6) + 4.0*PY2*d(9);
    A[30] = d(5)*3.0*PXY2 + 6.0*PY2*d(6);
    A[31] = A[6];
    A[32] = A[12];
    A[33] = A[18];
    A[34] = A[24];
    A[35] = A[30];
    A[36] = d(5)*9.0*PY2;
    TEMP = 4.0*AREA;
    for (int i=1;i<=36;i++) A[i]*=TEMP;
    
    for (int i=37;i<=72;i++) A[i]=0.0;
    A[37]=XBSQ;
    A[40]=XBSQ*XSUBB;
    A[44]=XSUBB;
    A[49]=-2.0*XSUBB;
    A[52]=-3.0*XBSQ;
    A[55]=XCSQ;
    A[56]=XCYC;
    A[57]=YCSQ;
    A[58]=XCSQ*XSUBC;
    A[59]=YCSQ*XSUBC;
    A[60]=YCSQ*YSUBC;
    A[62]=XSUBC;
    A[63]=YSUBC*2.0;
    A[65]=XCYC*2.0;
    A[66]=YCSQ*3.0;
    A[67]=-2.0*XSUBC;
    A[68]=-YSUBC;
    A[70]=-3.0*XCSQ;
    A[71]=-YCSQ;
    
    MatrixXd H(6,6);
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) H(i,j)=A[37+i+6*j];
    bool ising;
    MatrixXd Hinv = inverd(H, ising);
    if (ising) { S.NOGO=true; return MatrixXd::Zero(9,9); }
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) A[37+i+6*j]=Hinv(i,j);
    
    MatrixXd KX(6,6);
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) KX(i,j)=A[1+i+6*j];
    MatrixXd KQ = KX*Hinv;
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) A[73+i+6*j]=KQ(i,j);
    
    MatrixXd KII = Hinv.transpose()*KQ;
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) A[109+i+6*j]=KII(i,j);
    
    MatrixXd Sm(6,3);
    Sm << 1,0,-XSUBB,
          0,1,0,
          0,0,1,
          1,YSUBC,-XSUBC,
          0,1,0,
          0,0,1;
    
    MatrixXd KIA = KII*Sm;
    for(int i=0;i<6;i++) for(int j=0;j<3;j++) A[46+i+6*j]=KIA(i,j);
    
    MatrixXd KAA = Sm.transpose()*KIA;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) A[1+i+3*j]=KAA(i,j);
    
    for(int i=46;i<=63;i++) A[i]=-A[i];
    
    for(int i=28;i<=36;i++) A[i]=A[i+18];
    A[10]=A[46]; A[11]=A[49]; A[12]=A[52];
    A[13]=A[47]; A[14]=A[50]; A[15]=A[53];
    A[16]=A[48]; A[17]=A[51]; A[18]=A[54];
    A[19]=A[55]; A[20]=A[58]; A[21]=A[61];
    A[22]=A[56]; A[23]=A[59]; A[24]=A[62];
    A[25]=A[57]; A[26]=A[60]; A[27]=A[63];
    A[37]=A[109]; A[38]=A[110]; A[39]=A[111];
    A[40]=A[115]; A[41]=A[116]; A[42]=A[117];
    A[43]=A[121]; A[44]=A[122]; A[45]=A[123];
    A[46]=A[112]; A[47]=A[113]; A[48]=A[114];
    A[49]=A[118]; A[50]=A[119]; A[51]=A[120];
    A[52]=A[124]; A[53]=A[125]; A[54]=A[126];
    A[64]=A[127]; A[65]=A[128]; A[66]=A[129];
    A[67]=A[133]; A[68]=A[134]; A[69]=A[135];
    A[70]=A[139]; A[71