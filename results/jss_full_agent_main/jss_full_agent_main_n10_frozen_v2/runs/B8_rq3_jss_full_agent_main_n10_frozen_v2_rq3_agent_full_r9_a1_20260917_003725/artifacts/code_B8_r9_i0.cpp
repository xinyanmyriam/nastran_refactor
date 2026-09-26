#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Eigen;

static double E = 200e9, NU = 0.3, GG;

Matrix<double,6,6> buildG() {
    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double temp1 = (1.0+NU)*(1.0-2.0*NU);
    double g1 = E*(1.0-NU)/temp1;
    double g2 = E*NU/temp1;
    G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
    G(0,1)=g2; G(0,2)=g2;
    G(1,0)=g2; G(1,2)=g2;
    G(2,0)=g2; G(2,1)=g2;
    G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;
    return G;
}

double invert4(const Matrix4d& H, Matrix4d& Hinv) {
    Hinv = H.inverse();
    return H.determinant();
}

struct TetResult {
    Matrix3d K[4];
    bool valid;
};

TetResult ktetra(const Matrix<double,4,3>& coords, const Matrix<double,6,6>& Gmat,
                 double hdeterDiv, double hdeterMul, int pivotLocal) {
    TetResult res;
    res.valid = true;
    for (int i=0;i<4;i++) res.K[i] = Matrix3d::Zero();

    Matrix4d H;
    for (int i=0;i<4;i++) {
        H(i,0) = 1.0;
        H(i,1) = coords(i,0);
        H(i,2) = coords(i,1);
        H(i,3) = coords(i,2);
    }
    Matrix4d Hinv;
    double hdeter = invert4(H, Hinv);
    if (std::abs(hdeter) < 1e-30) { res.valid=false; return res; }
    hdeter = std::abs(hdeter);

    Matrix<double,6,3> C[4];
    for (int i=0;i<4;i++) {
        C[i] = Matrix<double,6,3>::Zero();
        C[i](0,0) = Hinv(i,1);
        C[i](4,0) = Hinv(i,2);
        C[i](2,1) = Hinv(i,3);
        C[i](4,1) = Hinv(i,3);
        C[i](5,1) = Hinv(i,2);
        C[i](0,2) = Hinv(i,3);
        C[i](2,2) = Hinv(i,1);
        C[i](3,2) = Hinv(i,2);
        C[i](4,2) = Hinv(i,1);
    }

    double hd = hdeter / hdeterDiv * hdeterMul;

    Matrix<double,3,6> GCT = C[pivotLocal].transpose() * Gmat;
    GCT = GCT * hd;

    for (int i=0;i<4;i++) {
        res.K[i] = GCT * C[i];
    }
    return res;
}

int main() {
    Matrix<double,6,3> nodes;
    nodes << 0,0,0,
             1,0,0,
             0,1,0,
             0,0,1,
             1,0,1,
             0,1,1;

    GG = E / (2.0*(1.0+NU));
    Matrix<double,6,6> Gmat = buildG();

    int M[12][4] = {
        {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
        {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
        {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
    };

    Matrix<double,18,18> K = Matrix<double,18,18>::Zero();

    for (int pivotNode = 0; pivotNode < 6; pivotNode++) {
        for (int t = 0; t < 12; t++) {
            Matrix<double,4,3> coords;
            int localPivot = -1;
            for (int j=0;j<4;j++) {
                int nodeIdx = M[t][j] - 1;
                coords.row(j) = nodes.row(nodeIdx);
                if (nodeIdx == pivotNode) localPivot = j;
            }
            if (localPivot < 0) continue;

            int IOPT = t + 11;
            double div = 36.0;
            double mul = (IOPT <= 16) ? 2.0 : 1.0;

            TetResult r = ktetra(coords, Gmat, div, mul, localPivot);
            if (!r.valid) continue;

            for (int i=0;i<4;i++) {
                int gi = M[t][i] - 1;
                K.block<3,3>(pivotNode*3, gi*3) += r.K[i];
            }
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<18;i++) {
        std::cout << "[";
        for (int j=0;j<18;j++) {
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if (j<17) std::cout << ",";
        }
        std::cout << "]";
        if (i<17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}