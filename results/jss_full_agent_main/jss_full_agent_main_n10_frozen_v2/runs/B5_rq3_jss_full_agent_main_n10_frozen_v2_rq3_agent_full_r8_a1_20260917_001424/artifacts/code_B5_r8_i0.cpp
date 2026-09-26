#define _USE_MATH_DEFINES
#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

struct State {
    double DEGRA = M_PI / 180.0;
    double G11, G12, G13, G22, G23, G33;
    double E, nu, t;
    double ECPT[26];
    int NECPT[26];
    double A[226];
    double T[10];
    double TEMP9[10];
    double PROD9[10];
    double ARR9[10];
    double ARRAY9[10];
    double KSUM[64];
    double XSUBB, XSUBC, YSUBC;
    double Evec[19];
    double TEMP, L1, L2, S1, S2, C1, C2, X1, X2, Y1, Y2;
    double TEMP1, TEMP2;
    double DETERM;
    int NPIVOT;
    double THETA;
    int ISING;
    int SUBSCA, SUBSCB, SUBSCC;
    double SINTH, COSTH;
    double SINANG, COSANG;
    double U1, U2;
    double IVECT[4], JVECT[4], KVECT[4];
    double V[3], VV[3];
    double R[3][5];
    double D1[4], D2[4];
    double HQ[13];
    double PROD12[13];
    double HABC[19];
    double G[37];
    double TITE[19];
    double TJTE[19];
    double KOUT[37];
    double TEMP18[19];
    double V1[4], V2[4], V3[4];
    double REQUIV[9];
    int M[10];
    double G2X211, G2X212, G2X222;
    double GSUBE;
    int IOPT4;
    int IFKGG, IF4GG;
    int K4GGSW;
    int NOGO;
    int NPVT;
};

State st;

void GMMATD(const double* A, int nra, int nca, int tra,
            const double* B, int nrb, int ncb, int trb,
            double* C) {
    int ra = tra ? nca : nra;
    int ca = tra ? nra : nca;
    int cb = trb ? nrb : ncb;
    for (int i = 0; i < ra; i++) {
        for (int j = 0; j < cb; j++) {
            double sum = 0.0;
            for (int k = 0; k < ca; k++) {
                double av, bv;
                if (tra) av = A[k + i*nca];
                else    av = A[i + k*nra];
                if (trb) bv = B[j + k*ncb];
                else     bv = B[k + j*nrb];
                sum += av * bv;
            }
            C[i + j*ra] = sum;
        }
    }
}

void INVERD(int n, double* A, int lda, double* , int ,
            double& determ, int& ising, double* ) {
    MatrixXd M(n, 2*n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            M(i,j) = A[i + j*lda];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            M(i, n+j) = (i==j) ? 1.0 : 0.0;
    determ = 1.0;
    ising = -1;
    for (int col = 0; col < n; col++) {
        int piv = col;
        double maxv = fabs(M(col,col));
        for (int r = col+1; r < n; r++)
            if (fabs(M(r,col)) > maxv) { maxv = fabs(M(r,col)); piv = r; }
        if (maxv < 1e-300) { ising = 2; return; }
        if (piv != col) { M.row(piv).swap(M.row(col)); determ = -determ; }
        determ *= M(col,col);
        double d = M(col,col);
        M.row(col) /= d;
        for (int r = 0; r < n; r++)
            if (r != col) {
                double f = M(r,col);
                if (f != 0.0) M.row(r) -= f * M.row(col);
            }
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i + j*lda] = M(i, n+j);
}

void MAT() {
    double E = st.E, nu = st.nu;
    double f = E / (1.0 - nu*nu);
    st.G11 = f; st.G12 = f*nu; st.G13 = 0.0;
    st.G22 = f; st.G23 = 0.0; st.G33 = f*(1.0-nu)/2.0;
}

void TRANSD(int , double* T) {
    for (int i = 0; i < 9; i++) T[i] = 0.0;
    T[0]=1.0; T[4]=1.0; T[8]=1.0;
}

void KTRBSC(int IOPT) {
    int NTYPE = (IOPT > 0) ? 1 : 0;
    if (NTYPE != 1) {
        st.Evec[11] = st.ECPT[18] - st.ECPT[14];
        st.Evec[14] = st.ECPT[19] - st.ECPT[15];
        st.Evec[17] = st.ECPT[20] - st.ECPT[16];
        st.XSUBB = sqrt(st.Evec[11]*st.Evec[11]+st.Evec[14]*st.Evec[14]+st.Evec[17]*st.Evec[17]);
        if (st.XSUBB <= 1e-6) { st.NOGO=1; return; }
        st.Evec[11]/=st.XSUBB; st.Evec[14]/=st.XSUBB; st.Evec[17]/=st.XSUBB;
        st.Evec[2]=st.ECPT[22]-st.ECPT[14];
        st.Evec[5]=st.ECPT[23]-st.ECPT[15];
        st.Evec[8]=st.ECPT[24]-st.ECPT[16];
        st.XSUBC=st.Evec[11]*st.Evec[2]+st.Evec[14]*st.Evec[5]+st.Evec[17]*st.Evec[8];
        st.Evec[1]=st.Evec[14]*st.Evec[8]-st.Evec[5]*st.Evec[17];
        st.Evec[4]=st.Evec[2]*st.Evec[17]-st.Evec[11]*st.Evec[8];
        st.Evec[7]=st.Evec[11]*st.Evec[5]-st.Evec[2]*st.Evec[14];
        st.YSUBC=sqrt(st.Evec[1]*st.Evec[1]+st.Evec[4]*st.Evec[4]+st.Evec[7]*st.Evec[7]);
        if (st.YSUBC <= 1e-6) { st.NOGO=1; return; }
        st.Evec[1]/=st.YSUBC; st.Evec[4]/=st.YSUBC; st.Evec[7]/=st.YSUBC;
        st.Evec[12]=st.Evec[4]*st.Evec[17]-st.Evec[14]*st.Evec[7];
        st.Evec[15]=st.Evec[11]*st.Evec[7]-st.Evec[1]*st.Evec[17];
        st.Evec[18]=st.Evec[1]*st.Evec[14]-st.Evec[11]*st.Evec[4];
        st.TEMP=sqrt(st.Evec[12]*st.Evec[12]+st.Evec[15]*st.Evec[15]+st.Evec[18]*st.Evec[18]);
        st.Evec[12]/=st.TEMP; st.Evec[15]/=st.TEMP; st.Evec[18]/=st.TEMP;
        st.Evec[2]=0; st.Evec[3]=0; st.Evec[5]=0; st.Evec[6]=0;
        st.Evec[8]=0; st.Evec[9]=0; st.Evec[10]=0; st.Evec[13]=0; st.Evec[16]=0;
        st.THETA=st.ECPT[5]*st.DEGRA;
        st.SINTH=sin(st.THETA); st.COSTH=cos(st.THETA);
        if (fabs(st.SINTH)<1e-6) st.SINTH=0.0;
    }
    MAT();
    double Gloc[10];
    Gloc[1]=st.G11; Gloc[2]=st.G12; Gloc[3]=st.G13;
    Gloc[4]=st.G12; Gloc[5]=st.G22; Gloc[6]=st.G23;
    Gloc[7]=st.G13; Gloc[8]=st.G23; Gloc[9]=st.G33;
    double D[10];
    double EYE = st.ECPT[7];
    for (int i=1;i<=9;i++) D[i]=Gloc[i]*EYE;
    double XSUBB=st.XSUBB, XSUBC=st.XSUBC, YSUBC=st.YSUBC;
    double AREA=XSUBB*YSUBC/2.0;
    double XBAR=(XSUBB+XSUBC)/3.0, YBAR=YSUBC/3.0;
    double XCSQ=XSUBC*XSUBC, YCSQ=YSUBC*YSUBC, XBSQ=XSUBB*XSUBB;
    double XCYC=XSUBC*YSUBC;
    double PX2=(XBSQ+XSUBB*XSUBC+XCSQ)/6.0;
    double PY2=YCSQ/6.0;
    double PXY2=YSUBC*(XSUBB+2.0*XSUBC)/12.0;
    double XBAR3=3.0*XBAR, YBAR3=3.0*YBAR, YBAR2=2.0*YBAR;
    double* A=st.A;
    A[1]=D[1]; A[2]=D[3]; A[3]=D[2];
    A[4]=D[1]*XBAR3; A[5]=D[2]*XBAR+YBAR2*D[3]; A[6]=D[2]*YBAR3;
    A[7]=A[2]; A[8]=D[9]; A[9]=D[6];
    A[10]=D[3]*XBAR3; A[11]=D[6]*XBAR+YBAR2*D[9]; A[12]=D[6]*YBAR3;
    A[13]=A[3]; A[14]=A[9]; A[15]=D[5];
    A[16]=D[2]*XBAR3; A[17]=D[5]*XBAR+YBAR2*D[6]; A[18]=D[5]*YBAR3;
    A[19]=A[4]; A[20]=A[10]; A[21]=A[16];
    A[22]=D[1]*9.0*PX2;
    A[23]=D[2]*3.0*PX2+6.0*PXY2*D[3];
    A[24]=D[2]*9.0*PXY2;
    A[25]=A[5]; A[26]=A[11]; A[27]=A[17]; A[28]=A[23];
    A[29]=D[5]*PX2+4.0*PXY2*D[6]+4.0*PY2*D[9];
    A[30]=D[5]*3.0*PXY2+6.0*PY2*D[6];
    A[31]=A[6]; A[32]=A[12]; A[33]=A[18]; A[34]=A[24]; A[35]=A[30];
    A[36]=D[5]*9.0*PY2;
    double TEMP=4.0*AREA;
    for (int i=1;i<=36;i++) A[i]*=TEMP;
    for (int i=37;i<=72;i++) A[i]=0.0;
    A[37]=XBSQ; A[40]=XBSQ*XSUBB; A[44]=XSUBB;
    A[49]=-2.0*XSUBB; A[52]=-3.0*XBSQ;
    A[55]=XCSQ; A[56]=XCYC; A[57]=YCSQ;
    A[58]=XCSQ*XSUBC; A[59]=YCSQ*XSUBC; A[60]=YCSQ*YSUBC;
    A[62]=XSUBC; A[63]=YSUBC*2.0;
    A[65]=XCYC*2.0; A[66]=YCSQ*3.0;
    A[67]=-2.0*XSUBC; A[68]=-YSUBC;
    A[70]=-3.0*XCSQ; A[71]=-YCSQ;
    double T2=st.ECPT[9];
    if (T2 != 0.0) {
        MAT();
        if (!(st.G2X211==0.0 && st.G2X212==0.0 && st.G2X222==0.0)) {
            double G2X2[5];
            G2X2[1]=st.G2X211*T2; G2X2[2]=st.G2X212*T2;
            G2X2[3]=st.G2X212*T2; G2X2[4]=st.G2X222*T2;
            double DETERM=G2X2[1]*G2X2[4]-G2X2[3]*G2X2[2];
            double J2X2[5];
            J2X2[1]=G2X2[4]/DETERM; J2X2[2]=-G2X2[2]/DETERM;
            J2X2[3]=-G2X2[3]/DETERM; J2X2[4]=G2X2[1]/DETERM;
            double TEMP2=2.0*D[2]+4.0*D[9];
            A[73]=-6.0*(J2X2[1]*D[1]+J2X2[2]*D[3]);
            A[74]=-J2X2[1]*TEMP2-6.0*J2X2[2]*D[6];
            A[75]=-6.0*(J2X2[1]*D[6]+J2X2[2]*D[5]);
            A[76]=-6.0*(J2X2[2]*D[1]+J2X2[4]*D[3]);
            A[77]=-J2X2[2]*TEMP2-6.0*J2X2[4]*D[6];
            A[78]=-6.0*(J2X2[2]*D[6]+J2X2[4]*D[5]);
            GMMATD(&G2X2[1],2,2,0, &A[73],2,3,0, &A[79]);
            GMMATD(&A[73],2,3,1, &A[79],2,3,0, &A[85]);
            for (int i=1;i<=3;i++) {
                A[i+21]+=A[i+84]*AREA;
                A[i+27]+=A[i+87]*AREA;
                A[i+33]+=A[i+90]*AREA;
            }
            for (int i=1;i<=3;i++) {
                A[i+39]+=XSUBB*A[i+72];
                A[i+57]+=XSUBC*A[i+72]+YSUBC*A[i+75];
            }
        }
    }
    st.ISING=-1;
    INVERD(6, &A[37], 6, &A[73], 0, st.DETERM, st.ISING, &A[79]);
    if (st.ISING==2) { st.NOGO=1; return; }
    GMMATD(&A[1],6,6,0, &A[37],6,6,0, &A[73]);
    GMMATD(&A[37],6,6,1, &A[73],6,6,0, &A[109]);
    double Sloc[19];
    if (IOPT==2) {
        for (int i=37;i<=72;i++) A[i+108]=A[i];
    }
    Sloc[1]=1.0; Sloc[2]=0.0; Sloc[3]=-XSUBB;
    Sloc[4]=0.0; Sloc[5]=1.0; Sloc[6]=0.0;
    Sloc[7]=0.0; Sloc[8]=0.0; Sloc[9]=1.0;
    Sloc[10]=1.0; Sloc[11]=YSUBC; Sloc[12]=-XSUBC;
    Sloc[13]=0.0; Sloc[14]=1.0; Sloc[15]=0.0;
    Sloc[16]=0.0; Sloc[17]=0.0; Sloc[18]=1.0;
    GMMATD(&A[109],6,6,0, &Sloc[1],6,3,0, &A[46]);
    GMMATD(&Sloc[1],6,3,1, &A[46],6,3,0, &A[1]);
    for (int i=46;i<=63;i++) A[i]=-A[i];
    for (int i=28;i<=36;i++) A[i]=A[i+18];
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
    if (NTYPE==1) return;
    for (int i=1;i<=3;i++) {
        if (st.NECPT[i]==st.NPVT) { st.NPIVOT=i; goto found; }
    }
    st.NOGO=1; return;
found:
    {
        int NPIVOT=st.NPIVOT;
        if (st.NECPT[4*NPIVOT+9]!=0) {
            double TI[10];
            TRANSD(st.NECPT[4*NPIVOT+9], TI);
            GMMATD(&TI[1],3,3,1, &st.Evec[1],3,3,0, &st.TITE[1]);
            GMMATD(&TI[1],3,3,1, &st.Evec[10],3,3,0, &st.TITE[10]);
        } else {
            for (int i=1;i<=18;i++) st.TITE[i]=st.Evec[i];
        }
        int NPT1=(NPIVOT==1)?28:1;
        for (int J=1;J<=3;J++) {
            if (st.NECPT[4*J+9]!=0) {
                double TI[10];
                TRANSD(st.NECPT[4*J+9], TI);
                GMMATD(&TI[1],3,3,1, &st.Evec[1],3,3,0, &st.TJTE[1]);
                GMMATD(&TI[1],3,3,1, &st.Evec[10],3,3,0, &st.TJTE[10]);
            } else {
                for (int i=1;i<=18;i++) st.TJTE[i]=st.Evec[i];
            }
            int NPT2=27*NPIVOT+9*J-35;
            GMMATD(&A[NPT2],3,3,0, &st.TJTE[1],6,3,1, &A[NPT1]);
            GMMATD(&st.TITE[1],6,3,0, &A[NPT1],3,6,0, &A[100]);
        }
    }
}

int main() {
    for (int i=0;i<226;i++) st.A[i]=0.0;
    for (int i=0;i<64;i++) st.KSUM[i]=0.0;
    for (int i=0;i<37;i++) st.G[i]=0.0;
    for (int i=0;i<10;i++) st.M[i]=0;
    st.M[1]=1; st.M[2]=2; st.M[3]=4;
    st.M[4]=2; st.M[5]=3; st.M[6]=4;
    st.M[7]=3; st.M[8]=1; st.M[9]=4;

    st.E=200e9; st.nu=0.3; st.t=0.01;
    double I=st.t*st.t*st.t/12.0;

    st.ECPT[1]=1; st.ECPT[2]=1; st.ECPT[3]=2; st.ECPT[4]=3;
    st.ECPT[5]=0.0; st.ECPT[6]=1; st.ECPT[7]=I; st.ECPT[8]=0; st.ECPT[9]=0.0;
    st.ECPT[10]=0.0; st.ECPT[11]=0.0; st.ECPT[12]=0.0;
    st.ECPT[13]=0;
    st.ECPT[14]=0.0; st.ECPT[15]=0.0; st.ECPT[16]=0.0;
    st.ECPT[17]=0;
    st.ECPT[18]=1.0; st.ECPT[19]=0.0; st.ECPT[20]=0.0;
    st.ECPT[21]=0;
    st.ECPT[22]=0.0; st.ECPT[23]=1.0; st.ECPT[24]=0.0;
    st.ECPT[25]=0.0;
    for (int i=1;i<=25;i++) st.NECPT[i]=(int)st.ECPT[i];

    st.NPVT=1; st.IOPT4=0; st.IFKGG=0; st.IF4GG=0;
    st.K4GGSW=0; st.NOGO=0; st.GSUBE=0.0;
    st.G2X211=0.0; st.G2X212=0.0; st.G2X222=0.0;

    st.THETA=st.ECPT[5]*st.DEGRA;
    st.SINANG=sin(st.THETA); st.COSANG=cos(st.THETA);

    for (int i=0;i<9;i++) st.REQUIV[i]=0.0;
    for (int i=1;i<=2;i++) for (int j=1;j<=4;j++) st.R[i][j]=0.0;

    st.V1[1]=st.ECPT[14]; st.V1[2]=st.ECPT[15]; st.V1[3]=st.ECPT[16];
    st.V2[1]=st.ECPT[18]; st.V2[2]=st.ECPT[19]; st.V2[3]=st.ECPT[20];
    st.V3[1]=st.ECPT[22]; st.V3[2]=st.ECPT[23]; st.V3[3]=st.ECPT[24];

    for (int i=1;i<=3;i++) {
        st.D2[i]=st.V2[i]-st.V1[i];
        st.D1[i]=st.V3[i]-st.V1[i];
    }
    st.R[1][2]=sqrt(st.D2[1]*st.D2[1]+st.D2[2]*st.D2[2]+st.D2[3]*st.D2[3]);
    for (int i=1;i<=3;i++) st.IVECT[i]=st.D2[i]/st.R[1][2];
    st.KVECT[1]=st.IVECT[2]*st.D1[3]-st.D1[2]*st.IVECT[3];
    st.KVECT[2]=st.IVECT[3]*st.D1[1]-st.D1[3]*st.IVECT[1];
    st.KVECT[3]=st.IVECT[1]*st.D1[2]-st.D1[1]*st.IVECT[2];
    st.R[2][3]=sqrt(st.KVECT[1]*st.KVECT[1]+st.KVECT[2]*st.KVECT[2]+st.KVECT[3]*st.KVECT[3]);
    for (int i=1;i<=3;i++) st.KVECT[i]/=st.R[2][3];
    st.JVECT[1]=st.KVECT[2]*st.IVECT[3]-st.IVECT[2]*st.KVECT[3];
    st.JVECT[2]=st.KVECT[3]*st.IVECT[1]-st.IVECT[3]*st.KVECT[1];
    st.JVECT[