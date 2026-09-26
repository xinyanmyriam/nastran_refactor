#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct MatOut {
    double G11, G12, G13, G22, G23, G33;
    double RHO, ALPHA1, ALPHA2, ALP12, T0, GSUBE;
    double SIGTEN, SIGCOM, SIGSHE, G2X211, G2X212, G2X222;
};

static MatOut matout;
static double DEGRA = M_PI / 180.0;

static void GMMATD(const double* A, int ra, int ca, int transa,
                   const double* B, int rb, int cb, int transb,
                   double* C) {
    int ar = transa ? ca : ra;
    int ac = transa ? ra : ca;
    int br = transb ? cb : rb;
    int bc = transb ? rb : cb;
    for (int i = 0; i < ar; ++i) {
        for (int j = 0; j < bc; ++j) {
            double s = 0.0;
            for (int k = 0; k < ac; ++k) {
                double a = transa ? A[k * ra + i] : A[i * ca + k];
                double b = transb ? B[j * rb + k] : B[k * cb + j];
                s += a * b;
            }
            C[i * bc + j] = s;
        }
    }
}

static int INVERD(int n, const double* A, int lda, double* B) {
    MatrixXd M(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M(i, j) = A[i * lda + j];
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) return 2;
    MatrixXd Minv = lu.inverse();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            B[i * lda + j] = Minv(i, j);
    return 1;
}

static double g_E = 200e9, g_nu = 0.3;

static void MAT(int inflag) {
    double E = g_E, nu = g_nu;
    double f = E / (1.0 - nu * nu);
    matout.G11 = f;
    matout.G22 = f;
    matout.G33 = f * (1.0 - nu) / 2.0;
    matout.G12 = nu * f;
    matout.G13 = 0.0;
    matout.G23 = 0.0;
    matout.G2X211 = 0.0;
    matout.G2X212 = 0.0;
    matout.G2X222 = 0.0;
    matout.GSUBE = 0.0;
}

static void KTRBSC(int IOPT, double XSUBB, double XSUBC, double YSUBC,
                   double EYE, double T2, double SINTH, double COSTH,
                   double* A) {
    double G[9], D[9];
    MAT(2);
    G[0] = matout.G11; G[1] = matout.G12; G[2] = matout.G13;
    G[3] = matout.G12; G[4] = matout.G22; G[5] = matout.G23;
    G[6] = matout.G13; G[7] = matout.G23; G[8] = matout.G33;
    for (int i = 0; i < 9; ++i) D[i] = G[i] * EYE;

    double AREA = XSUBB * YSUBC / 2.0;
    double XBAR = (XSUBB + XSUBC) / 3.0;
    double YBAR = YSUBC / 3.0;
    double XCSQ = XSUBC * XSUBC;
    double YCSQ = YSUBC * YSUBC;
    double XBSQ = XSUBB * XSUBB;
    double XCYC = XSUBC * YSUBC;
    double PX2 = (XBSQ + XSUBB * XSUBC + XCSQ) / 6.0;
    double PY2 = YCSQ / 6.0;
    double PXY2 = YSUBC * (XSUBB + 2.0 * XSUBC) / 12.0;
    double XBAR3 = 3.0 * XBAR;
    double YBAR3 = 3.0 * YBAR;
    double YBAR2 = 2.0 * YBAR;

    A[0] = D[0];
    A[1] = D[2];
    A[2] = D[1];
    A[3] = D[0] * XBAR3;
    A[4] = D[1] * XBAR + YBAR2 * D[2];
    A[5] = D[1] * YBAR3;
    A[6] = A[1];
    A[7] = D[8];
    A[8] = D[5];
    A[9] = D[2] * XBAR3;
    A[10] = D[5] * XBAR + YBAR2 * D[8];
    A[11] = D[5] * YBAR3;
    A[12] = A[2];
    A[13] = A[8];
    A[14] = D[4];
    A[15] = D[1] * XBAR3;
    A[16] = D[4] * XBAR + YBAR2 * D[5];
    A[17] = D[4] * YBAR3;
    A[18] = A[3];
    A[19] = A[9];
    A[20] = A[15];
    A[21] = D[0] * 9.0 * PX2;
    A[22] = D[1] * 3.0 * PX2 + 6.0 * PXY2 * D[2];
    A[23] = D[1] * 9.0 * PXY2;
    A[24] = A[4];
    A[25] = A[10];
    A[26] = A[16];
    A[27] = A[22];
    A[28] = D[4] * PX2 + 4.0 * PXY2 * D[5] + 4.0 * PY2 * D[8];
    A[29] = D[4] * 3.0 * PXY2 + 6.0 * PY2 * D[5];
    A[30] = A[5];
    A[31] = A[11];
    A[32] = A[17];
    A[33] = A[23];
    A[34] = A[29];
    A[35] = D[4] * 9.0 * PY2;
    double TEMP = 4.0 * AREA;
    for (int i = 0; i < 36; ++i) A[i] *= TEMP;

    for (int i = 36; i < 72; ++i) A[i] = 0.0;
    A[36] = XBSQ;
    A[39] = XBSQ * XSUBB;
    A[43] = XSUBB;
    A[48] = -2.0 * XSUBB;
    A[51] = -3.0 * XBSQ;
    A[54] = XCSQ;
    A[55] = XCYC;
    A[56] = YCSQ;
    A[57] = XCSQ * XSUBC;
    A[58] = YCSQ * XSUBC;
    A[59] = YCSQ * YSUBC;
    A[61] = XSUBC;
    A[62] = YSUBC * 2.0;
    A[64] = XCYC * 2.0;
    A[65] = YCSQ * 3.0;
    A[66] = -2.0 * XSUBC;
    A[67] = -YSUBC;
    A[69] = -3.0 * XCSQ;
    A[70] = -YCSQ;

    if (T2 != 0.0) {
        MAT(3);
        if (!(matout.G2X211 == 0.0 && matout.G2X212 == 0.0 && matout.G2X222 == 0.0)) {
            double G2X2[4], J2X2[4];
            G2X2[0] = matout.G2X211 * T2;
            G2X2[1] = matout.G2X212 * T2;
            G2X2[2] = matout.G2X212 * T2;
            G2X2[3] = matout.G2X222 * T2;
            double DETERM = G2X2[0] * G2X2[3] - G2X2[2] * G2X2[1];
            J2X2[0] = G2X2[3] / DETERM;
            J2X2[1] = -G2X2[1] / DETERM;
            J2X2[2] = -G2X2[2] / DETERM;
            J2X2[3] = G2X2[0] / DETERM;
            double T2v = 2.0 * D[1] + 4.0 * D[8];
            A[72] = -6.0 * (J2X2[0] * D[0] + J2X2[1] * D[2]);
            A[73] = -J2X2[0] * T2v - 6.0 * J2X2[1] * D[5];
            A[74] = -6.0 * (J2X2[0] * D[5] + J2X2[1] * D[4]);
            A[75] = -6.0 * (J2X2[1] * D[0] + J2X2[3] * D[2]);
            A[76] = -J2X2[1] * T2v - 6.0 * J2X2[3] * D[5];
            A[77] = -6.0 * (J2X2[1] * D[5] + J2X2[3] * D[4]);
            GMMATD(G2X2, 2, 2, 0, &A[72], 2, 3, 0, &A[78]);
            GMMATD(&A[72], 2, 3, 1, &A[78], 2, 3, 0, &A[84]);
            for (int i = 0; i < 3; ++i) {
                A[i + 21] += A[i + 84] * AREA;
                A[i + 27] += A[i + 87] * AREA;
                A[i + 33] += A[i + 90] * AREA;
            }
            for (int i = 0; i < 3; ++i) {
                A[i + 39] += XSUBB * A[i + 72];
                A[i + 57] += XSUBC * A[i + 72] + YSUBC * A[i + 75];
            }
        }
    }

    {
        double Hinv[36];
        int ising = INVERD(6, &A[36], 6, Hinv);
        for (int i = 0; i < 36; ++i) A[36 + i] = Hinv[i];
        if (ising == 2) return;
    }

    GMMATD(&A[0], 6, 6, 0, &A[36], 6, 6, 0, &A[72]);
    GMMATD(&A[36], 6, 6, 1, &A[72], 6, 6, 0, &A[108]);

    if (IOPT == 2) {
        for (int i = 36; i < 72; ++i) A[i + 108] = A[i];
    }

    double S[18];
    S[0] = 1.0; S[1] = 0.0; S[2] = -XSUBB;
    S[3] = 0.0; S[4] = 1.0; S[5] = 0.0;
    S[6] = 0.0; S[7] = 0.0; S[8] = 1.0;
    S[9] = 1.0; S[10] = YSUBC; S[11] = -XSUBC;
    S[12] = 0.0; S[13] = 1.0; S[14] = 0.0;
    S[15] = 0.0; S[16] = 0.0; S[17] = 1.0;
    for (int i = 0; i < 18; ++i) A[81 + i] = S[i];

    GMMATD(&A[108], 6, 6, 0, S, 6, 3, 0, &A[45]);
    GMMATD(S, 6, 3, 1, &A[45], 6, 3, 0, &A[0]);
    for (int i = 45; i < 63; ++i) A[i] = -A[i];

    for (int i = 27; i < 36; ++i) A[i] = A[i + 18];
    A[9] = A[45]; A[10] = A[48]; A[11] = A[51];
    A[12] = A[46]; A[13] = A[49]; A[14] = A[52];
    A[15] = A[47]; A[16] = A[50]; A[17] = A[53];
    A[18] = A[54]; A[19] = A[57]; A[20] = A[60];
    A[21] = A[55]; A[22] = A[58]; A[23] = A[61];
    A[24] = A[56]; A[25] = A[59]; A[26] = A[62];
    A[36] = A[108]; A[37] = A[109]; A[38] = A[110];
    A[39] = A[114]; A[40] = A[115]; A[41] = A[116];
    A[42] = A[120]; A[43] = A[121]; A[44] = A[122];
    A[45] = A[111]; A[46] = A[112]; A[47] = A[113];
    A[48] = A[117]; A[49] = A[118]; A[50] = A[119];
    A[51] = A[123]; A[52] = A[124]; A[53] = A[125];
    A[63] = A[126]; A[64] = A[127]; A[65] = A[128];
    A[66] = A[132]; A[67] = A[133]; A[68] = A[134];
    A[69] = A[138]; A[70] = A[139]; A[71] = A[140];
    A[72] = A[129]; A[73] = A[130]; A[74] = A[131];
    A[75] = A[135]; A[76] = A[136]; A[77] = A[137];
    A[78] = A[141]; A[79] = A[142]; A[80] = A[143];
}

int main() {
    double t = 0.01;
    double EYE = t * t * t / 12.0;
    g_E = 200e9;
    g_nu = 0.3;

    double V1[3] = {0.0, 0.0, 0.0};
    double V2[3] = {1.0, 0.0, 0.0};
    double V3[3] = {0.0, 1.0, 0.0};

    double R[2][4];
    for (int i = 0; i < 8; ++i) (&R[0][0])[i] = 0.0;

    double D2[3], D1[3];
    for (int i = 0; i < 3; ++i) {
        D2[i] = V2[i] - V1[i];
        D1[i] = V3[i] - V1[i];
    }

    R[0][1] = sqrt(D2[0]*D2[0] + D2[1]*D2[1] + D2[2]*D2[2]);
    double IVECT[3];
    for (int i = 0; i < 3; ++i) IVECT[i] = D2[i] / R[0][1];

    double KVECT[3];
    KVECT[0] = IVECT[1]*D1[2] - D1[1]*IVECT[2];
    KVECT[1] = IVECT[2]*D1[0] - D1[2]*IVECT[0];
    KVECT[2] = IVECT[0]*D1[1] - D1[0]*IVECT[1];

    R[1][2] = sqrt(KVECT[0]*KVECT[0] + KVECT[1]*KVECT[1] + KVECT[2]*KVECT[2]);
    for (int i = 0; i < 3; ++i) KVECT[i] /= R[1][2];

    double JVECT[3];
    JVECT[0] = KVECT[1]*IVECT[2] - IVECT[1]*KVECT[2];
    JVECT[1] = KVECT[2]*IVECT[0] - IVECT[2]*KVECT[0];
    JVECT[2] = KVECT[0]*IVECT[1] - IVECT[0]*KVECT[1];
    double TEMP = sqrt(JVECT[0]*JVECT[0] + JVECT[1]*JVECT[1] + JVECT[2]*JVECT[2]);
    for (int i = 0; i < 3; ++i) JVECT[i] /= TEMP;

    R[0][2] = D1[0]*IVECT[0] + D1[1]*IVECT[1] + D1[2]*IVECT[2];
    R[0][3] = (R[0][1] + R[0][2]) / 3.0;
    R[1][3] = R[1][2] / 3.0;

    int M[9] = {1,2,4, 2,3,4, 3,1,4};

    double KSUM[63];
    double G[36];
    for (int i = 0; i < 63; ++i) KSUM[i] = 0.0;
    for (int i = 0; i < 36; ++i) G[i] = 0.0;

    double THETA = 0.0;
    double SINANG = sin(THETA);
    double COSANG = cos(THETA);

    double A[226];

    for (int J = 1; J <= 3; ++J) {
        int KM = 3*J - 3;
        int SUBSCA = M[KM+0];
        int SUBSCB = M[KM+1];
        int SUBSCC = M[KM+2];

        double V[2], VV[2];
        for (int i = 0; i < 2; ++i) {
            V[i] = R[i][SUBSCB-1] - R[i][SUBSCA-1];
            VV[i] = R[i][SUBSCC-1] - R[i][SUBSCA-1];
        }
        double XSUBB = sqrt(V[0]*V[0] + V[1]*V[1]);
        double U1 = V[0] / XSUBB;
        double U2 = V[1] / XSUBB;
        double XSUBC = U1*VV[0] + U2*VV[1];
        double YSUBC = U1*VV[1] - U2*VV[0];

        double SINTH = SINANG*U1 - COSANG*U2;
        double COSTH = COSANG*U1 + SINANG*U2;
        if (fabs(SINTH) < 1.0e-6) SINTH = 0.0;

        for (int i = 0; i < 226; ++i) A[i] = 0.0;
        KTRBSC(2, XSUBB, XSUBC, YSUBC, EYE, 0.0, SINTH, COSTH, A);

        double T[9];
        T[0] = 1.0; T[1] = 0.0; T[2] = 0.0;
        T[3] = 0.0; T[4] = U1;  T[5] = U2;
        T[6] = 0.0; T[7] = -U2; T[8] = U1;

        for (int I = 1; I <= 3; ++I) {
            double TEMP9[9], PROD9[9];
            GMMATD(T, 3, 3, 1, &A[27*I - 9], 3, 3, 0, TEMP9);
            GMMATD(TEMP9, 3, 3, 0, T, 3, 3, 0, PROD9);
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] + 18;
            for (int K = 1; K <= 9; ++K) {
                int NSUBC = NPOINT + K;
                KSUM[NSUBC-1] += PROD9[K-1];
            }
        }

        for (int K = 1; K <= 2; ++K) {
            int NPOINT = KM + K;
            if (M[NPOINT-1] != 1) continue;
            double TEMP9[9], PROD9[9];
            GMMATD(T, 3, 3, 1, &A[36*K - 36], 3, 3, 0, TEMP9);
            GMMATD(TEMP9, 3, 3, 0, T, 3, 3, 0, PROD9);
            NPOINT = 9*1 - 9;
            for (int I = 1; I <= 9; ++I) {
                int NSUBC = NPOINT + I;
                KSUM[NSUBC-1] += PROD9[I-1];
            }
            GMMATD(T, 3, 3, 1, &A[18*K - 9], 3, 3, 0, TEMP9);
            GMMATD(TEMP9, 3, 3, 0, T, 3, 3, 0, PROD9);
            NPOINT = KM + 3 - K;
            NPOINT = 9*M[NPOINT-1] - 9;
            for (int I = 1; I <= 9; ++I) {
                int NSUBC = NPOINT + I;
                KSUM[NSUBC-1] += PROD9[I-1];
            }
        }

        double TEMP1 = XSUBB - XSUBC;
        double TEMP2 = YSUBC * YSUBC;
        double L1 = sqrt(XSUBC*XSUBC + TEMP2);
        double L2 = sqrt(TEMP1*TEMP1 + TEMP2);
        double S1 = XSUBC / L1;
        double S2 = TEMP1 / L2;
        double C1 = YSUBC / L1;
        double C2 = YSUBC / L2;
        double X1 = XSUBC / 2.0;
        double Y1 = YSUBC / 2.0;
        double X2 = (XSUBB + XSUBC) / 2.0;
        double Y2 = Y1;
        double HQ[12];
        HQ[0] = -XSUBC * C1;
        HQ[1] = X1*S1 - Y1*C1;
        HQ[2] = 2.0*Y1*S1;
        HQ[3] = -3.0*X1*X1*C1;
        HQ[4] = Y1*(2.0*X1*S1 - Y1*C1);
        HQ[5] = 3.0*Y1*Y1*S1;
        HQ[6] = 2.0*X2*C2;
        HQ[7] = X2*S2 + Y2*C2;
        HQ[8] = 2.0*Y2*S2;
        HQ[9] = 3.0*X2*X2*C2;
        HQ[10] = Y2*(2.0*X2*S2 + Y2*C2);
        HQ[11] = 3.0*Y2*Y2*S2;

        double PROD12[12], HABC[18];
        GMMATD(HQ, 2, 6, 0, &A[144], 6, 6, 0, PROD12);
        GMMATD(PROD12, 2, 6, 0, &A[81], 6, 3, 0, HABC);

        HABC[0] = -HABC[0];
        HABC[1] = -HABC[1] + S1;
        HABC[2] = -HABC[2] + C1;
        HABC[3] = -HABC[3];
        HABC[4] = -HABC[4] + S2;
        HABC[5] = -HABC[5] - C2;

        HABC[6] = PROD12[0];
        HABC[7] = PROD12[1];
        HABC[8] = PROD12[2];
        HABC[9] = PROD12[6];
        HABC[10] = PROD12[7];
        HABC[11] = PROD12[8];
        HABC[12] = PROD12[3];
        HABC[13] = PROD12[4];
        HABC[14] = PROD12[5];
        HABC[15] = PROD12[9];
        HABC[16] = PROD12[10];
        HABC[17] = PROD12[11];

        for (int I = 1; I <= 3; ++I) {
            double TEMP9[9];
            GMMATD(&HABC[6*I - 6], 2, 3, 0, T, 3, 3, 0, TEMP9);
            int NPOINT = KM + I;
            NPOINT = 9*M[NPOINT-1] - 9;
            if (J == 2) {
                NPOINT += 3;
                for (int K = 1; K <= 6; ++K) {
                    NPOINT += 1;
                    G[NPOINT-1] += TEMP9[K-1];
                }
            } else if (J == 1) {
                for (int K = 1; K <= 6; ++K) {
                    NPOINT += 1;
                    G[NPOINT-1] += TEMP9[K-1];
                }
            } else {
                G[NPOINT + 7 - 1] += TEMP9[0];
                G[NPOINT + 8 - 1] += TEMP9[1];
                G[NPOINT + 9 - 1] += TEMP9[2];
                G[NPOINT + 1 - 1] += TEMP9[3];
                G[NPOINT + 2 - 1] += TEMP9[4];
                G[NPOINT + 3 - 1] += TEMP9[5];
            }
        }
    }

    double E[18];
    for (int i = 0; i < 18; ++i) E[i] = 0.0;
    E[0] = KVECT[0]; E[3] = KVECT[1]; E[6] = KVECT[2];
    E[10] = IVECT[0]; E[13] = IVECT[1]; E[16] = IVECT[2];
    E[11] = JVECT[0]; E[14] = JVECT[1]; E[17] = JVECT[2];

    double TITE[18];
    for (int i = 0; i < 18; ++i) TITE[i] = E[i];

    double PROD9[9], TEMP9[9];
    int ising = INVERD(3, &G[27], 3, PROD9);
    if (ising == 2) { cerr << "singular G4\n"; return 1; }

    GMMATD(&G[27], 3, 3, 0, &G[0], 3, 3, 0, PROD9);
    GMMATD(PROD9, 3, 3, 1, &KSUM[54], 3, 3, 0, TEMP9);

    for (int J = 1; J <= 3; ++J) {
        double ARR9[9], ARRAY9[9];
        GMMATD(PROD9, 3, 3, 1, &KSUM[9*J + 18], 3, 3, 1, ARR9);
        int NBEGIN = 9*J - 9;
        for (int I = 1; I <= 9; ++I) {
            int NPOINT = NBEGIN + I;
            KSUM[NPOINT-1] -= ARR9[I-1];
        }
        GMMATD(&G[27], 3, 3, 0, &G[9*J - 9], 3, 3, 0, ARR9);
        GMMATD(&KSUM[9*1 + 18], 3, 3, 0, ARR9, 3, 3, 0, ARRAY9);
        for (int I = 1; I <= 9; ++I) {
            int NPOINT = NBEGIN + I;
            KSUM[NPOINT-1] -= ARRAY9[I-1];
        }
        GMMATD(TEMP9, 3, 3, 0, ARR9, 3, 3, 0, ARRAY9);
        for (int I = 1; I <= 9; ++I) {
            int NPOINT = NBEGIN + I;
            KSUM[NPOINT-1] += ARRAY9[I-1];
        }
    }

    MatrixXd K = MatrixXd::Zero(9, 9);
    for (int J = 1; J <= 3; ++J) {
        for (int I = 1; I <= 3; ++I) {
            int b = M[3*(J-1) + (I-1)];
            int base = 9*b + 18;
            for (int r = 0; r < 3; ++r) {
                for (int c = 0; c < 3; ++c) {
                    K((I-1)*3 + r, (J-1)*3 + c) = KSUM[base + r*3 + c - 1];
                }
            }
        }
    }

    cout << "{\"stiffness_matrix\":[";
    cout << setprecision(6) << scientific;
    for (int i = 0; i < 9; ++i) {
        cout << "[";
        for (int j = 0; j < 9; ++j) {
            cout << K(i, j);
            if (j < 8) cout << ",";
        }
        cout << "]";
        if (i < 8) cout << ",";
    }
    cout << "]}" << endl;

    return 0;
}