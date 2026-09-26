#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

static const double PI = 3.14159265358979323846;

double G11, G12, G13, G22, G23, G33;

MatrixXd inverd(const MatrixXd& M, int& ising) {
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = 2; return MatrixXd::Zero(M.rows(), M.cols()); }
    ising = 1;
    return lu.inverse();
}

struct KTRBSCResult {
    MatrixXd KU[3][3];
    MatrixXd Hinv;
    MatrixXd S;
    double XSUBB, XSUBC, YSUBC;
    bool ok;
};

KTRBSCResult ktrbsc(double X1,double Y1,double Z1,
                    double X2,double Y2,double Z2,
                    double X3,double Y3,double Z3,
                    double EYE, double ANGLE, double T2,
                    double g11,double g12,double g13,
                    double g22,double g23,double g33) {
    KTRBSCResult res;
    res.ok = true;
    Vector3d Ivec, Jvec, Kvec;
    Vector3d d2(X2-X1, Y2-Y1, Z2-Z1);
    double XSUBB = d2.norm();
    if (XSUBB <= 1.0e-6) { res.ok=false; return res; }
    Ivec = d2 / XSUBB;
    Vector3d d1(X3-X1, Y3-Y1, Z3-Z1);
    double XSUBC = Ivec.dot(d1);
    Kvec = Ivec.cross(d1);
    double YSUBC = Kvec.norm();
    if (YSUBC <= 1.0e-6) { res.ok=false; return res; }
    Kvec = Kvec / YSUBC;
    Jvec = Kvec.cross(Ivec);
    Jvec.normalize();
    res.XSUBB = XSUBB; res.XSUBC = XSUBC; res.YSUBC = YSUBC;
    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    double XCSQ = XSUBC*XSUBC, YCSQ = YSUBC*YSUBC, XBSQ = XSUBB*XSUBB;
    double XCYC = XSUBC*YSUBC;
    double PX2 = (XBSQ + XSUBB*XSUBC + XCSQ)/6.0;
    double PY2 = YCSQ/6.0;
    double PXY2 = YSUBC*(XSUBB + 2.0*XSUBC)/12.0;
    double XBAR3 = 3.0*XBAR, YBAR3 = 3.0*YBAR, YBAR2 = 2.0*YBAR;
    double Dv[10];
    Dv[1]=g11*EYE; Dv[2]=g12*EYE; Dv[3]=g13*EYE;
    Dv[4]=g12*EYE; Dv[5]=g22*EYE; Dv[6]=g23*EYE;
    Dv[7]=g13*EYE; Dv[8]=g23*EYE; Dv[9]=g33*EYE;
    double A[226];
    for (int i=0;i<226;i++) A[i]=0.0;
    A[ 1]=Dv[1]; A[ 2]=Dv[3]; A[ 3]=Dv[2];
    A[ 4]=Dv[1]*XBAR3; A[ 5]=Dv[2]*XBAR+YBAR2*Dv[3]; A[ 6]=Dv[2]*YBAR3;
    A[ 7]=A[2]; A[ 8]=Dv[9]; A[ 9]=Dv[6];
    A[10]=Dv[3]*XBAR3; A[11]=Dv[6]*XBAR+YBAR2*Dv[9]; A[12]=Dv[6]*YBAR3;
    A[13]=A[3]; A[14]=A[9]; A[15]=Dv[5];
    A[16]=Dv[2]*XBAR3; A[17]=Dv[5]*XBAR+YBAR2*Dv[6]; A[18]=Dv[5]*YBAR3;
    A[19]=A[4]; A[20]=A[10]; A[21]=A[16];
    A[22]=Dv[1]*9.0*PX2; A[23]=Dv[2]*3.0*PX2+6.0*PXY2*Dv[3]; A[24]=Dv[2]*9.0*PXY2;
    A[25]=A[5]; A[26]=A[11]; A[27]=A[17]; A[28]=A[23];
    A[29]=Dv[5]*PX2+4.0*PXY2*Dv[6]+4.0*PY2*Dv[9];
    A[30]=Dv[5]*3.0*PXY2+6.0*PY2*Dv[6];
    A[31]=A[6]; A[32]=A[12]; A[33]=A[18]; A[34]=A[24]; A[35]=A[30];
    A[36]=Dv[5]*9.0*PY2;
    double TEMP = 4.0*AREA;
    for (int i=1;i<=36;i++) A[i]*=TEMP;
    for (int i=37;i<=72;i++) A[i]=0.0;
    A[37]=XBSQ; A[40]=XBSQ*XSUBB; A[44]=XSUBB; A[49]=-2.0*XSUBB; A[52]=-3.0*XBSQ;
    A[55]=XCSQ; A[56]=XCYC; A[57]=YCSQ;
    A[58]=XCSQ*XSUBC; A[59]=YCSQ*XSUBC; A[60]=YCSQ*YSUBC;
    A[62]=XSUBC; A[63]=YSUBC*2.0; A[65]=XCYC*2.0; A[66]=YCSQ*3.0;
    A[67]=-2.0*XSUBC; A[68]=-YSUBC; A[70]=-3.0*XCSQ; A[71]=-YCSQ;
    MatrixXd H(6,6);
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) H(r,c)=A[37+c*6+r];
    int ising;
    MatrixXd Hinv = inverd(H, ising);
    if (ising==2) { res.ok=false; return res; }
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) A[37+c*6+r]=Hinv(r,c);
    MatrixXd KX(6,6);
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) KX(r,c)=A[1+c*6+r];
    MatrixXd KQHinv = KX * Hinv;
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) A[73+c*6+r]=KQHinv(r,c);
    MatrixXd KII = Hinv.transpose() * KQHinv;
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) A[109+c*6+r]=KII(r,c);
    for (int i=37;i<=72;i++) A[i+108]=A[i];
    double S[19];
    S[1]=1; S[2]=0; S[3]=-XSUBB; S[4]=0; S[5]=1; S[6]=0; S[7]=0; S[8]=0; S[9]=1;
    S[10]=1; S[11]=YSUBC; S[12]=-XSUBC; S[13]=0; S[14]=1; S[15]=0; S[16]=0; S[17]=0; S[18]=1;
    for (int i=1;i<=18;i++) A[81+i]=S[i];
    MatrixXd Smat(6,3);
    for (int c=0;c<3;c++) for (int r=0;r<6;r++) Smat(r,c)=S[1+c*6+r];
    MatrixXd KIA = -KII * Smat;
    for (int c=0;c<3;c++) for (int r=0;r<6;r++) A[46+c*6+r]=KIA(r,c);
    MatrixXd KAA = Smat.transpose() * KIA;
    for (int c=0;c<3;c++) for (int r=0;r<3;r++) A[1+c*3+r]=KAA(r,c);
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
    for (int i=0;i<3;i++) for (int j=0;j<3;j++) {
        int base = 1 + (i*3+j)*9;
        Matrix3d M;
        for (int c=0;c<3;c++) for (int r=0;r<3;r++) M(r,c)=A[base+c*3+r];
        res.KU[i][j]=M;
    }
    MatrixXd HinvOut(6,6);
    for (int c=0;c<6;c++) for (int r=0;r<6;r++) HinvOut(r,c)=A[145+c*6+r];
    res.Hinv=HinvOut;
    MatrixXd Sout(6,3);
    for (int c=0;c<3;c++) for (int r=0;r<6;r++) Sout(r,c)=A[82+c*6+r];
    res.S=Sout;
    return res;
}

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double I = t*t*t/12.0;
    double denom = 1.0 - nu*nu;
    G11 = E/denom; G12 = nu*E/denom; G13 = 0.0;
    G22 = E/denom; G23 = 0.0; G33 = E/(2.0*(1.0+nu));
    double X1=0,Y1=0,Z1=0, X2=1,Y2=0,Z2=0, X3=0,Y3=1,Z3=0;
    double ANGLE=0.0, T2=0.0;
    KTRBSCResult r = ktrbsc(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,I,ANGLE,T2,G11,G12,G13,G22,G23,G33);
    if (!r.ok) { cerr << "KTRBSC failed\n"; return 1; }
    Vector3d V1(X1,Y1,Z1), V2(X2,Y2,Z2), V3(X3,Y3,Z3);
    Vector3d D2 = V2 - V1, D1 = V3 - V1;
    double R[2][5];
    for (int i=0;i<2;i++) for (int j=0;j<5;j++) R[i][j]=0.0;
    R[0][1] = D2.norm();
    Vector3d IVECT = D2 / R[0][1];
    Vector3d KVECT = IVECT.cross(D1);
    R[1][2] = KVECT.norm();
    KVECT = KVECT / R[1][2];
    Vector3d JVECT = KVECT.cross(IVECT);
    JVECT.normalize();
    R[0][2] = D1.dot(IVECT);
    R[0][3] = (R[0][1] + R[0][2]) / 3.0;
    R[1][3] = R[1][2] / 3.0;
    int M[10];
    M[1]=1; M[2]=2; M[3]=4; M[4]=2; M[5]=3; M[6]=4; M[7]=3; M[8]=1; M[9]=4;
    double KSUM[64]; for (int i=0;i<64;i++) KSUM[i]=0.0;
    double Gmat[37]; for (int i=0;i<37;i++) Gmat[i]=0.0;
    double THETA = ANGLE * (PI/180.0);
    double SINANG = sin(THETA), COSANG = cos(THETA);
    for (int J=1; J<=3; J++) {
        int KM = 3*J - 3;
        int SUBSCA = M[KM+1], SUBSCB = M[KM+2], SUBSCC = M[KM+3];
        double Vv[3], VVv[3];
        for (int i=0;i<2;i++) {
            Vv[i] = R[i][SUBSCB-1] - R[i][SUBSCA-1];
            VVv[i]= R[i][SUBSCC-1] - R[i][SUBSCA-1];
        }
        double XSUBB = sqrt(Vv[0]*Vv[0] + Vv[1]*Vv[1]);
        double U1 = Vv[0]/XSUBB, U2 = Vv[1]/XSUBB;
        double XSUBC = U1*VVv[0] + U2*VVv[1];
        double YSUBC = U1*VVv[1] - U2*VVv[0];
        double SINTH = SINANG*U1 - COSANG*U2;
        double COSTH = COSANG*U1 + SINANG*U2;
        if (fabs(SINTH) < 1.0e-6) SINTH = 0.0;
        KTRBSCResult sub = ktrbsc(0,0,0, XSUBB,0,0, XSUBC,YSUBC,0,I,ANGLE,T2,G11,G12,G13,G22,G23,G33);
        if (!sub.ok) { cerr << "sub failed\n"; return 1; }
        Matrix3d T;
        T << 1.0, 0.0, 0.0, 0.0, U1, -U2, 0.0, U2, U1;
        for (int i=0;i<3;i++) {
            Matrix3d prod = T.transpose() * sub.KU[i][i] * T;
            int NPOINT = KM + (i+1);
            NPOINT = 9*M[NPOINT] + 18;
            for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++)
                KSUM[NPOINT + c*3 + rr + 1] += prod(rr,c);
        }
        int NPIVOT = 1;
        for (int K=1; K<=2; K++) {
            int NPOINT = KM + K;
            if (M[NPOINT] != NPIVOT) continue;
            Matrix3d prod1 = T.transpose() * sub.KU[K-1][2] * T;
            int NP = 9*NPIVOT - 9;
            for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++)
                KSUM[NP + c*3 + rr + 1] += prod1(rr,c);
            Matrix3d prod2 = T.transpose() * sub.KU[2][K-1] * T;
            int NP2 = KM + 3 - K;
            NP2 = 9*M[NP2] - 9;
            for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++)
                KSUM[NP2 + c*3 + rr + 1] += prod2(rr,c);
        }
        double TEMP1 = XSUBB - XSUBC;
        double TEMP2 = YSUBC*YSUBC;
        double L1 = sqrt(XSUBC*XSUBC + TEMP2);
        double L2 = sqrt(TEMP1*TEMP1 + TEMP2);
        double S1 = XSUBC/L1, S2 = TEMP1/L2;
        double C1 = YSUBC/L1, C2 = YSUBC/L2;
        double X1h = XSUBC/2.0, Y1h = YSUBC/2.0;
        double X2h = (XSUBB + XSUBC)/2.0, Y2h = Y1h;
        double HQ[13];
        HQ[ 1] = -XSUBC*C1;
        HQ[ 2] = X1h*S1 - Y1h*C1;
        HQ[ 3] = 2.0*Y1h*S1;
        HQ[ 4] = -3.0*X1h*X1h*C1;
        HQ[ 5] = Y1h*(2.0*X1h*S1 - Y1h*C1);
        HQ[ 6] = 3.0*Y1h*Y1h*S1;
        HQ[ 7] = 2.0*X2h*C2;
        HQ[ 8] = X2h*S2 + Y2h*C2;
        HQ[ 9] = 2.0*Y2h*S2;
        HQ[10] = 3.0*X2h*X2h*C2;
        HQ[11] = Y2h*(2.0*X2h*S2 + Y2h*C2);
        HQ[12] = 3.0*Y2h*Y2h*S2;
        MatrixXd HQm(2,6);
        for (int c=0;c<6;c++) for (int rr=0;rr<2;rr++) HQm(rr,c)=HQ[1+c*2+rr];
        MatrixXd PROD12 = HQm * sub.Hinv;
        MatrixXd PS = PROD12 * sub.S;
        double HABC[19];
        HABC[1] = -PS(0,0);
        HABC[2] = -PS(0,1) + S1;
        HABC[3] = -PS(0,2) + C1;
        HABC[4] = -PS(1,0);
        HABC[5] = -PS(1,1) + S2;
        HABC[6] = -PS(1,2) - C2;
        double P12[13];
        for (int c=0;c<6;c++) for (int rr=0;rr<2;rr++) P12[1+c*2+rr]=PROD12(rr,c);
        HABC[ 7]=P12[ 1]; HABC[ 8]=P12[ 2]; HABC[ 9]=P12[ 3];
        HABC[10]=P12[ 7]; HABC[11]=P12[ 8]; HABC[12]=P12[ 9];
        HABC[13]=P12[ 4]; HABC[14]=P12[ 5]; HABC[15]=P12[ 6];
        HABC[16]=P12[10]; HABC[17]=P12[11]; HABC[18]=P12[12];
        for (int i=1;i<=3;i++) {
            MatrixXd Hi(2,3);
            for (int c=0;c<3;c++) for (int rr=0;rr<2;rr++) Hi(rr,c)=HABC[6*i-5 + c*2 + rr];
            MatrixXd TEMP9m = Hi * T;
            int NPOINT = KM + i;
            NPOINT = 9*M[NPOINT] - 9;
            if (J == 2) {
                NPOINT += 3;
                for (int k=0;k<6;k++) {
                    NPOINT += 1;
                    int cc = k/2, rr = k%2;
                    Gmat[NPOINT] += TEMP9m(rr,cc);
                }
            } else if (J == 1) {
                for (int k=0;k<6;k++) {
                    NPOINT += 1;
                    int cc = k/2, rr = k%2;
                    Gmat[NPOINT] += TEMP9m(rr,cc);
                }
            } else {
                Gmat[NPOINT + 7] += TEMP9m(0,0);
                Gmat[NPOINT + 8] += TEMP9m(0,1);
                Gmat[NPOINT + 9] += TEMP9m(0,2);
                Gmat[NPOINT + 1] += TEMP9m(1,0);
                Gmat[NPOINT + 2] += TEMP9m(1,1);
                Gmat[NPOINT + 3] += TEMP9m(1,2);
            }
        }
    }
    double Emat[19];
    for (int i=1;i<=18;i++) Emat[i]=0.0;
    Emat[ 1]=KVECT(0); Emat[ 4]=KVECT(1); Emat[ 7]=KVECT(2);
    Emat[11]=IVECT(0); Emat[14]=IVECT(1); Emat[17]=IVECT(2);
    Emat[12]=JVECT(0); Emat[15]=JVECT(1); Emat[18]=JVECT(2);
    double TITE[19];
    for (int i=1;i<=18;i++) TITE[i]=Emat[i];
    Matrix3d G4;
    for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) G4(rr,c)=Gmat[28+c*3+rr];
    int ising;
    Matrix3d G4inv = inverd(G4, ising);
    if (ising==2) { cerr << "G4 singular\n"; return 1; }
    int NPIVOT = 1;
    Matrix3d Gpiv;
    for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) Gpiv(rr,c)=Gmat[9*NPIVOT-8+c*3+rr];
    Matrix3d PROD9 = G4inv * Gpiv;
    Matrix3d K44;
    for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) K44(rr,c)=KSUM[55+c*3+rr];
    Matrix3d TEMP9 = PROD9.transpose() * K44;
    for (int J=1;J<=3;J++) {
        Matrix3d KJ4;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KJ4(rr,c)=KSUM[9*J+19+c*3+rr];
        Matrix3d ARR9 = PROD9.transpose() * KJ4.transpose();
        int NBEGIN = 9*J-9;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++)
            KSUM[NBEGIN + c*3 + rr + 1] -= ARR9(rr,c);
        Matrix3d GJ;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) GJ(rr,c)=Gmat[9*J-8+c*3+rr];
        ARR9 = G4inv * GJ;
        Matrix3d KI4;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) KI4(rr,c)=KSUM[9*NPIVOT+19+c*3+rr];
        Matrix3d ARRAY9 = KI4 * ARR9;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++)
            KSUM[NBEGIN + c*3 + rr + 1] -= ARRAY9(rr,c);
        ARRAY9 = TEMP9 * ARR9;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++)
            KSUM[NBEGIN + c*3 + rr + 1] += ARRAY9(rr,c);
    }
    auto getBlock = [&](int base) {
        Matrix3d M;
        for (int c=0;c<3;c++) for (int rr=0;rr<3;rr++) M(rr,c)=KSUM[base+c*3+rr];
        return M;
    };
    Matrix3d K11 = getBlock(1);
    Matrix3d K12 = getBlock(10);
    Matrix3d K13 = getBlock(19);
    Matrix3d K22 = getBlock(37);
    Matrix3d K23 = getBlock(46);
    Matrix3d K33 = getBlock(55);
    MatrixXd K9 = MatrixXd::Zero(9,9);
    K9.block<3,3>(0,0)=K11; K9.block<3,3>(0,3)=K12; K9.block<3,3>(0,6)=K13;
    K9.block<3,3>(3,0)=K12.transpose(); K9.block<3,3>(3,3)=K22; K9.block<3,3>(3,6)=K23;
    K9.block<3,3>(6,0)=K13.transpose(); K9.block<3,3>(6,3)=K23.transpose(); K9.block<3,3>(6,6)=K33;
    cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<9;i++) {
        cout << "[";
        for (int j=0;j<9;j++) {
            cout << scientific << setprecision(6) << K9(i,j);
            if (j<8) cout << ",";
        }
        cout << "]";
        if (i<8) cout << ",";
    }
    cout << "]}" << endl;
    return 0;
}