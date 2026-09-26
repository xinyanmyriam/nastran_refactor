#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Eigen;

static double invert4(const Matrix4d& H, Matrix4d& Hinv) {
    double det = H.determinant();
    if (std::abs(det) < 1e-300) return 0.0;
    Hinv = H.inverse();
    return det;
}

int main() {
    double coords[6][3] = {
        {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}
    };
    double E = 200e9, nu = 0.3;
    double GG = E/(2.0*(1.0+nu));

    int M[12][4] = {
        {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
        {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
        {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
    };

    MatrixXd K = MatrixXd::Zero(18,18);

    for (int itet = 0; itet < 12; ++itet) {
        int IOPT = (itet+1) + 10;
        Matrix4d H;
        for (int i = 0; i < 4; ++i) {
            int node = M[itet][i] - 1;
            H(i,0) = 1.0;
            H(i,1) = coords[node][0];
            H(i,2) = coords[node][1];
            H(i,3) = coords[node][2];
        }
        Matrix4d Hinv;
        double HDETER = invert4(H, Hinv);
        HDETER = std::abs(HDETER);

        Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
        double TEMP1 = (1.0+nu)*(1.0-2.0*nu);
        double g1 = E*(1.0-nu)/TEMP1;
        double g2 = E*nu/TEMP1;
        G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
        G(0,1)=g2; G(0,2)=g2;
        G(1,0)=g2; G(1,2)=g2;
        G(2,0)=g2; G(2,1)=g2;
        G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

        std::vector<Matrix<double,6,3>> Cblocks(4);
        for (int I = 1; I <= 4; ++I) {
            Matrix<double,6,3> Cb = Matrix<double,6,3>::Zero();
            double hI4  = Hinv(I-1, 1);
            double hI8  = Hinv(I-1, 2);
            double hI12 = Hinv(I-1, 3);
            Cb(0,0) = hI4;
            Cb(4,0) = hI8;
            Cb(2,1) = hI12;
            Cb(4,1) = hI12;
            Cb(5,1) = hI8;
            Cb(0,2) = hI12;
            Cb(2,2) = hI4;
            Cb(3,2) = hI8;
            Cb(4,2) = hI4;
            Cblocks[I-1] = Cb;
        }

        double hd = HDETER/36.0;
        if (IOPT <= 16) hd *= 2.0;

        Matrix<double,12,12> Kloc = Matrix<double,12,12>::Zero();
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                Matrix3d kij = Cblocks[i].transpose() * G * Cblocks[j] * hd;
                Kloc.block<3,3>(3*i, 3*j) = kij;
            }
        }

        for (int i = 0; i < 4; ++i) {
            int gi = M[itet][i] - 1;
            for (int j = 0; j < 4; ++j) {
                int gj = M[itet][j] - 1;
                K.block<3,3>(3*gi, 3*gj) += Kloc.block<3,3>(3*i, 3*j);
            }
        }
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; ++i) {
        std::cout << "[";
        for (int j = 0; j < 18; ++j) {
            std::cout << K(i,j);
            if (j < 17) std::cout << ",";
        }
        std::cout << "]";
        if (i < 17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}