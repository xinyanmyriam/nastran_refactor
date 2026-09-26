#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct MatOut { double G11,G12,G13,G22,G23,G33; };

MatOut computeMat(double E, double nu, double t) {
    MatOut m;
    double D = E*t*t*t/(12.0*(1.0-nu*nu));
    m.G11 = D; m.G22 = D; m.G33 = D*(1.0-nu)/2.0;
    m.G12 = nu*D; m.G13 = 0.0; m.G23 = 0.0;
    return m;
}

void ktrbsc(double x1,double y1,double z1,
            double x2,double y2,double z2,
            double x3,double y3,double z3,
            double E, double nu, double t, double EYE, double T2,
            double angle, vector<double>& Aout)
{
    vector<double> A(226, 0.0);
    vector<double> Evec(19, 0.0);
    
    Evec[11] = x2 - x1; Evec[14] = y2 - y1; Evec[17] = z2 - z1;
    double XSUBB = sqrt(Evec[11]*Evec[11]+Evec[14]*Evec[14]+Evec[17]*Evec[17]);
    Evec[11]/=XSUBB; Evec[14]/=XSUBB; Evec[17]/=XSUBB;
    
    Evec[2] = x3 - x1; Evec[5] = y3 - y1; Evec[8] = z3 - z1;
    double XSUBC = Evec[11]*Evec[2]+Evec[14]*Evec[5]+Evec[17]*Evec[8];
    
    Evec[1] = Evec[14]*Evec[8] - Evec[5]*Evec[17];
    Evec[4] = Evec[2]*Evec[17] - Evec[11]*Evec[8];
    Evec[7] = Evec[11]*Evec[5] - Evec[2]*Evec[14];
    double YSUBC = sqrt(Evec[1]*Evec[1]+Evec[4]*Evec[4]+Evec[7]*Evec[7]);
    Evec[1]/=YSUBC; Evec[4]/=YSUBC; Evec[7]/=YSUBC;
    
    Evec[12] = Evec[4]*Evec[17] - Evec[14]*Evec[7];
    Evec[15] = Evec[11]*Evec[7] - Evec[1]*Evec[17];
    Evec[18] = Evec[1]*Evec[14] - Evec[11]*Evec[4];
    double TEMP = sqrt(Evec[12]*Evec[12]+Evec[15]*Evec[15]+Evec[18]*Evec[18]);
    Evec[12]/=TEMP; Evec[15]/=TEMP; Evec[18]/=TEMP;
    Evec[2]=0;Evec[3]=0;Evec[5]=0;Evec[6]=0;Evec[8]=0;Evec[9]=0;Evec[10]=0;Evec[13]=0;Evec[16]=0;
    
    MatOut mo = computeMat(E, nu, t);
    vector<double> G(10,0.0);
    G[1]=mo.G11;G[2]=mo.G12;G[3]=mo.G13;
    G[4]=mo.G12;G[5]=mo.G22;G[6]=mo.G23;
    G[7]=mo.G13;G[8]=mo.G23;G[9]=mo.G33;
    
    vector<double> D(10,0.0);
    for(int i=1;i<=9;i++) D[i]=G[i]*EYE;
    
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
    
    A[1]=D[1]; A[2]=D[3]; A[3]=D[2];
    A[4]=D[1]*XBAR3; A[5]=D[2]*XBAR+YBAR2*D[3]; A[6]=D[2]*YBAR3;
    A[7]=A[2]; A[8]=D[9]; A[9]=D[6];
    A[10]=D[3]*XBAR3; A[11]=D[6]*XBAR+YBAR2*D[9]; A[12]=D[6]*YBAR3;
    A[13]=A[3]; A[14]=A[9]; A[15]=D[5];
    A[16]=D[2]*XBAR3; A[17]=D[5]*XBAR+YBAR2*D[6]; A[18]=D[5]*YBAR3;
    A[19]=A[4]; A[20]=A[10]; A[21]=A[16];
    A[22]=D[1]*9.0*PX2; A[23]=D[2]*3.0*PX2+6.0*PXY2*D[3]; A[24]=D[2]*9.0*PXY2;
    A[25]=A[5]; A[26]=A[11]; A[27]=A[17];
    A[28]=A[23]; A[29]=D[5]*PX2+4.0*PXY2*D[6]+4.0*PY2*D[9]; A[30]=D[5]*3.0*PXY2+6.0*PY2*D[6];
    A[31]=A[6]; A[32]=A[12]; A[33]=A[18];
    A[34]=A[24]; A[35]=A[30]; A[36]=D[5]*9.0*PY2;
    TEMP = 4.0*AREA;
    for(int i=1;i<=36;i++) A[i]*=TEMP;
    
    for(int i=37;i<=72;i++) A[i]=0.0;
    A[37]=XBSQ; A[40]=XBSQ*XSUBB; A[44]=XSUBB;
    A[49]=-2.0*XSUBB; A[52]=-3.0*XBSQ;
    A[55]=XCSQ; A[56]=XCYC; A[57]=YCSQ;
    A[58]=XCSQ*XSUBC; A[59]=YCSQ*XSUBC; A[60]=YCSQ*YSUBC;
    A[62]=XSUBC; A[63]=YSUBC*2.0;
    A[65]=XCYC*2.0; A[66]=YCSQ*3.0;
    A[67]=-2.0*XSUBC; A[68]=-YSUBC;
    A[70]=-3.0*XCSQ; A[71]=-YCSQ;
    
    if(T2 != 0.0) {
        MatOut mo2 = computeMat(E, nu, t);
        double g211=mo2.G11, g212=mo2.G12, g222=mo2.G22;
        if(!(g211==0 && g212==0 && g222==0)) {
            double G2X2[4];
            G2X2[0]=g211*T2; G2X2[1]=g212*T2; G2X2[2]=g212*T2; G2X2[3]=g222*T2;
            double DETERM = G2X2[0]*G2X2[3]-G2X2[2]*G2X2[1];
            double J2X2[4];
            J2X2[0]=G2X2[3]/DETERM; J2X2[1]=-G2X2[1]/DETERM;
            J2X2[2]=-G2X2[2]/DETERM; J2X2[3]=G2X2[0]/DETERM;
            
            TEMP = 2.0*D[2]+4.0*D[9];
            A[73]=-6.0*(J2X2[0]*D[1]+J2X2[1]*D[3]);
            A[74]=-J2X2[0]*TEMP-6.0*J2X2[1]*D[6];
            A[75]=-6.0*(J2X2[0]*D[6]+J2X2[1]*D[5]);
            A[76]=-6.0*(J2X2[1]*D[1]+J2X2[3]*D[3]);
            A[77]=-J2X2[1]*TEMP-6.0*J2X2[3]*D[6];
            A[78]=-6.0*(J2X2[1]*D[6]+J2X2[3]*D[5]);
            
            MatrixXd G2(2,2), H(2,3);
            G2<<G2X2[0],G2X2[1],G2X2[2],G2X2[3];
            H<<A[73],A[74],A[75],A[76],A[77],A[78];
            MatrixXd prod = G2*H;
            for(int i=0;i<6;i++) A[79+i]=prod(i/3,i%3);
            
            MatrixXd prodT = H.transpose()*prod;
            for(int i=0;i<9;i++) A[85+i]=prodT(i/3,i%3);
            
            for(int i=1;i<=3;i++) {
                A[i+21]+=A[i+84]*AREA;
                A[i+27]+=A[i+87]*AREA;
                A[i+33]+=A[i+90]*AREA;
            }
            for(int i=1;i<=3;i++) {
                A[i+39]+=XSUBB*A[i+72];
                A[i+57]+=XSUBC*A[i+72]+YSUBC*A[i+75];
            }
        }
    }
    
    MatrixXd Hm(6,6);
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) Hm(i,j)=A[37+i*6+j];
    MatrixXd Hinv = Hm.inverse();
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) A[37+i*6+j]=Hinv(i,j);
    
    MatrixXd KX(6,6);
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) KX(i,j)=A[1+i*6+j];
    MatrixXd KQHinv = KX*Hinv;
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) A[73+i*6+j]=KQHinv(i,j);
    
    MatrixXd KII = Hinv.transpose()*KQHinv;
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) A[109+i*6+j]=KII(i,j);
    
    A[82]=1.0; A[83]=0.0; A[84]=-XSUBB;
    A[85]=0.0; A[86]=1.0; A[87]=0.0;
    A[88]=0.0; A[89]=0.0; A[90]=1.0;
    A[91]=1.0; A[92]=YSUBC; A[93]=-XSUBC;
    A[94]=0.0; A[95]=1.0; A[96]=0.0;
    A[97]=0.0; A[98]=0.0; A[99]=1.0;
    
    MatrixXd Sm(6,3);
    for(int i=0;i<6;i++) for(int j=0;j<3;j++) Sm(i,j)=A[82+i*3+j];
    MatrixXd KIA = KII*Sm;
    for(int i=0;i<6;i++) for(int j=0;j<3;j++) A[46+i*3+j]=KIA(i,j);
    
    MatrixXd KAA = Sm.transpose()*KIA;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++) A[1+i*3+j]=KAA(i,j);
    
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
    A[70]=A[139]; A[71]=A[140]; A[72]=A[141];
    A[73]=A[130]; A[74]=A[131]; A[75]=A[132];
    A[76]=A[136]; A[77]=A[137]; A[78]=A[138];
    A[79]=A[142]; A[80]=A[143]; A[81]=A[144];
    
    Aout = A;
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double I = t*t*t/12.0;
    
    double X[4] = {0,1,1,0};
    double Y[4] = {0,0,1,1};
    double Z[4] = {0,0,0,0};
    
    int M[12] = {2,4,1, 3,1,2, 4,2,3, 1,3,4};
    
    int NPIVOT = 1;
    int JNOT = NPIVOT + 2;
    
    double D1[3], D2[3], A1[3];
    for(int i=0;i<3;i++) {
        D1[i] = (i==0?X[2]:i==1?Y[2]:Z[2]) - (i==0?X[0]:i==1?Y[0]:Z[0]);
        D2[i] = (i==0?X[3]:i==1?Y[3]:Z[3]) - (i==0?X[1]:i==1?Y[1]:Z[1]);
        A1[i] = (i==0?X[1]:i==1?Y[1]:Z[1]) - (i==0?X[0]:i==1?Y[0]:Z[0]);
    }
    
    double KV[3], IV[3], JV[3];
    KV[0] = D1[1]*D2[2]-D2[1]*D1[2];
    KV[1] = D1[2]*D2[0]-D2[2]*D1[0];
    KV[2] = D1[0]*D2[1]-D2[0]*D1[1];
    double TEMP = sqrt(KV[0]*KV[0]+KV[1]*KV[1]+KV[2]*KV[2]);
    for(int i=0;i<3;i++) KV[i]/=TEMP;
    
    TEMP = (A1[0]*KV[0]+A1[1]*KV[1]+A1[2]*KV[2])/2.0;
    for(int i=0;i<3;i++) IV[i] = A1[i]-TEMP*KV[i];
    TEMP = sqrt(IV[0]*IV[0]+IV[1]*IV[1]+IV[2]*IV[2]);
    for(int i=0;i<3;i++) IV[i]/=TEMP;
    
    JV[0] = KV[1]*IV[2]-IV[1]*KV[2];
    JV[1] = KV[2]*IV[0]-IV[2]*KV[0];
    JV[2] = KV[0]*IV[1]-IV[0]*KV[1];
    TEMP = sqrt(JV[0]*JV[0]+JV[1]*JV[1]+JV[2]*JV[2]);
    for(int i=0;i<3;i++) JV[i]/=TEMP;
    
    double R[2][5];
    for(int i=0;i<2;i++) for(int j=0;j<5;j++) R[i][j]=0;
    R[0][3] = D1[0]*IV[0]+D1[1]*IV[1]+D1[2]*IV[2];
    R[0][2] = A1[0]*IV[0]+A1[1]*IV[1]+A1[2]*IV[2];
    R[1][3] = D1[0]*JV[0]+D1[1]*JV[1]+D1[2]*JV[2];
    R[0][4] = D2[0]*IV[0]+D2[1]*IV[1]+D2[2]*IV[2]+R[0][2];
    R[1][4] = D2[0]*JV[0]+D2[1]*JV[1]+D2[2]*JV[2];
    
    vector<double> KSUM(37, 0.0);
    
    for(int J=1; J<=4; J++) {
        if(J == JNOT) continue;
        int KM = 3*J-3;
        int SUBSCA = M[KM+0];
        int SUBSCB = M[KM+1];
        int SUBSCC = M[KM+2];
        
        double V[2], VV[2];
        for(int i=0;i<2;i++) {
            V[i] = R[i][SUBSCB]-R[i][SUBSCA];
            VV[i] = R[i][SUBSCC]-R[i][SUBSCA];
        }
        double XSUBB = sqrt(V[0]*V[0]+V[1]*V[1]);
        double U1 = V[0]/XSUBB;
        double U2 = V[1]/XSUBB;
        double XSUBC = U1*VV[0]+U2*VV[1];
        double YSUBC = U1*VV[1]-U2*VV[0];
        
        vector<double> Aout;
        ktrbsc(0,0,0, XSUBB,0,0, XSUBC,YSUBC,0, E, nu, t, I, 0.0, 0.0, Aout);
        
        Matrix3d T;
        T << 1,0,0, 0,U1,U2, 0,-U2,U1;
        
        int NBEGIN = 0;
        for(int i=1;i<=3;i++) {
            int NPOINT = KM+i;
            if(M[NPOINT-1] == NPIVOT) { NBEGIN = 27*i-27; break; }
        }
        
        for(int i=1;i<=3;i++) {
            int NPOINT = NBEGIN + 9*i - 8;
            Matrix3d Kij;
            for(int r=0;r<3;r++) for(int c=0;c<3;c++) Kij(r,c)=Aout[NPOINT+r*3+c];
            Matrix3d PROD9 = T.transpose()*Kij*T;
            
            int NPOINT2 = KM+i;
            int gnode = M[NPOINT2-1];
            int base = 9*gnode-9;
            for(int k=0;k<9;k++) KSUM[base+k+1] += PROD9(k/3,k%3)/2.0;
        }
    }
    
    // Build full 12x12 from KSUM (pivot row) - but we need full matrix
    // KSUM is pivot row (3x12). For full matrix, use symmetry.
    MatrixXd Kfull = MatrixXd::Zero(12,12);
    for(int r=0;r<3;r++) for(int c=0;c<12;c++) {
        Kfull(r,c) = KSUM[r*12+c+1];
    }
    // Symmetrize
    for(int i=0;i<12;i++) for(int j=0;j<i;j++) Kfull(i,j)=Kfull(j,i);
    
    cout << "{\"stiffness_matrix\":[";
    for(int i=0;i<12;i++) {
        cout << "[";
        for(int j=0;j<12;j++) {
            cout << std::scientific << std::setprecision(6) << Kfull(i,j);
            if(j<11) cout << ",";
        }
        cout << "]";
        if(i<11) cout << ",";
    }
    cout << "]}" << endl;
    
    return 0;
}