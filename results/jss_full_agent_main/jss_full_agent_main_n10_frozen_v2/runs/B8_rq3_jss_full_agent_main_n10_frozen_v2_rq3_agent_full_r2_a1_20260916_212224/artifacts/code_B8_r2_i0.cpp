#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Eigen;

struct State {
    double E, nu, GG;
    int NPVT;
    int IOPT4;
    double GSUBE;
    double VOLUME, SURFAC;
    bool HEAT, HYDRO;
    int NGPT, DIREC, KOUNT;
    double TVOL;
    int IDFLAG;
    int IDELEM;
    int IGFLAG;
    bool NOGOO;
};

static MatrixXd Kglobal;
static std::vector<int> silToLocal;

static void SMA1B_insert(const Matrix3d& kij, int silRow, int silCol) {
    int li = -1, lj = -1;
    for (int i = 0; i < 6; i++) {
        if (silToLocal[i] == silRow) li = i;
        if (silToLocal[i] == silCol) lj = i;
    }
    if (li < 0 || lj < 0) return;
    Kglobal.block<3,3>(3*li, 3*lj) += kij;
}

static void KTETRA(State& st, std::vector<double>& ECPT, int IOPT, int JTYPE) {
    auto NECPT = [&](int i) -> int { return (int)std::round(ECPT[i-1]); };

    Matrix4d H;
    H(0,0)=1; H(0,1)=ECPT[7]; H(0,2)=ECPT[8]; H(0,3)=ECPT[9];
    H(1,0)=1; H(1,1)=ECPT[11]; H(1,2)=ECPT[12]; H(1,3)=ECPT[13];
    H(2,0)=1; H(2,1)=ECPT[15]; H(2,2)=ECPT[16]; H(2,3)=ECPT[17];
    H(3,0)=1; H(3,1)=ECPT[19]; H(3,2)=ECPT[20]; H(3,3)=ECPT[21];

    double HDETER = H.determinant();
    Matrix4d Hinv = H.inverse();

    if (std::abs(HDETER) < 1e-30) { st.NOGOO = true; return; }
    if (IOPT >= 100) IOPT -= 100;

    bool pivotFound = false;
    for (int i = 3; i <= 6; i++) if (NECPT(i) == st.NPVT) { pivotFound = true; break; }
    if (!pivotFound) return;

    HDETER = std::abs(HDETER);

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0+st.nu)*(1.0-2.0*st.nu);
    double G1 = st.E*(1.0-st.nu)/TEMP1;
    double G2 = st.E*st.nu/TEMP1;
    G(0,0)=G1; G(1,1)=G1; G(2,2)=G1;
    G(0,1)=G2; G(0,2)=G2;
    G(1,0)=G2; G(1,2)=G2;
    G(2,0)=G2; G(2,1)=G2;
    G(3,3)=st.GG; G(4,4)=st.GG; G(5,5)=st.GG;

    std::vector<Matrix<double,6,3>> C(4);
    for (int i = 0; i < 4; i++) {
        C[i].setZero();
        C[i](0,0) = Hinv(i, 1);
        C[i](4,0) = Hinv(i, 2);
        C[i](2,1) = Hinv(i, 3);
        C[i](4,1) = Hinv(i, 3);
        C[i](5,1) = Hinv(i, 2);
        C[i](0,2) = Hinv(i, 3);
        C[i](2,2) = Hinv(i, 1);
        C[i](3,2) = Hinv(i, 2);
        C[i](4,2) = Hinv(i, 1);
    }

    if (IOPT >= 11 && IOPT <= 22) {
        HDETER = HDETER/36.0;
        if (IOPT <= 16) HDETER = HDETER*2.0;
    } else {
        HDETER = HDETER/12.0;
    }

    int pivotIdx = -1;
    for (int i = 2; i <= 5; i++) if (NECPT(i+1) == st.NPVT) { pivotIdx = i - 2; break; }
    if (pivotIdx < 0) return;

    Matrix<double,6,3> Cp = C[pivotIdx];
    int KA = 4*(pivotIdx+2) - 1;
    int csidPivot = NECPT(KA);

    Matrix<double,3,6> GCT;
    if (csidPivot != 0) {
        Matrix3d T = Matrix3d::Identity();
        Matrix<double,3,6> CT = T.transpose() * Cp.transpose();
        GCT = CT * G;
    } else {
        GCT = Cp.transpose() * G;
    }
    GCT *= HDETER;

    for (int i = 1; i <= 4; i++) {
        int csidNode = NECPT(4*i+3);
        Matrix3d KIJ;
        if (csidNode != 0) {
            Matrix3d T = Matrix3d::Identity();
            Matrix<double,6,3> CT = C[i-1] * T;
            KIJ = GCT * CT;
        } else {
            KIJ = GCT * C[i-1];
        }
        int sil = NECPT(i+2);
        SMA1B_insert(KIJ, st.NPVT, sil);
    }
}

int main() {
    std::vector<double> ECPT(100, 0.0);
    ECPT[0] = 1.0; ECPT[1] = 1.0;
    ECPT[2] = 1.0; ECPT[3] = 2.0; ECPT[4] = 3.0;
    ECPT[5] = 4.0; ECPT[6] = 5.0; ECPT[7] = 6.0;
    ECPT[8] = 0.0;  ECPT[9] = 0.0;  ECPT[10] = 0.0; ECPT[11] = 0.0;
    ECPT[12] = 0.0; ECPT[13] = 1.0; ECPT[14] = 0.0; ECPT[15] = 0.0;
    ECPT[16] = 0.0; ECPT[17] = 0.0; ECPT[18] = 1.0; ECPT[19] = 0.0;
    ECPT[20] = 0.0; ECPT[21] = 0.0; ECPT[22] = 0.0; ECPT[23] = 1.0;
    ECPT[24] = 0.0; ECPT[25] = 1.0; ECPT[26] = 0.0; ECPT[27] = 1.0;
    ECPT[28] = 0.0; ECPT[29] = 0.0; ECPT[30] = 1.0; ECPT[31] = 1.0;
    ECPT[32] = 0.0;

    State st;
    st.E = 200e9; st.nu = 0.3;
    st.GG = st.E / (2.0*(1.0+st.nu));
    st.NOGOO = false; st.IOPT4 = 0; st.GSUBE = 0.0;
    st.VOLUME = 0.0; st.SURFAC = 0.0;
    st.HEAT = false; st.HYDRO = false;
    st.NGPT = 99; st.DIREC = 0; st.KOUNT = 0;
    st.TVOL = 0.0; st.IDFLAG = 0; st.IDELEM = 0; st.IGFLAG = 0;

    silToLocal = {1,2,3,4,5,6};
    Kglobal = MatrixXd::Zero(18,18);

    for (int j = 0; j < 50; j++) ECPT[j+50] = ECPT[j];

    int M[12][4] = {
        {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
        {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
        {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
    };

    int NGRIDS = 6;
    int ITET = 1, NTET = 12;

    for (int pivot = 1; pivot <= 6; pivot++) {
        st.NPVT = pivot;
        for (int I = ITET; I <= NTET; I++) {
            int JTYPE = 1;
            if (I == NTET) JTYPE = -1;
            int IOPT = I + 10;
            for (int J = 1; J <= 4; J++) {
                int KPOINT = M[I-1][J-1];
                ECPT[J+1] = ECPT[KPOINT+51];
                KPOINT = 4*KPOINT + NGRIDS - 3;
                int JPOINT = 4*J + 2;
                ECPT[JPOINT]   = ECPT[KPOINT+51];
                ECPT[JPOINT+1] = ECPT[KPOINT+52];
                ECPT[JPOINT+2] = ECPT[KPOINT+53];
                ECPT[JPOINT+3] = ECPT[KPOINT+54];
            }
            KTETRA(st, ECPT, IOPT, JTYPE);
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; i++) {
        std::cout << "[";
        for (int j = 0; j < 18; j++) {
            std::cout << std::scientific << std::setprecision(6) << Kglobal(i,j);
            if (j < 17) std::cout << ",";
        }
        std::cout << "]";
        if (i < 17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}