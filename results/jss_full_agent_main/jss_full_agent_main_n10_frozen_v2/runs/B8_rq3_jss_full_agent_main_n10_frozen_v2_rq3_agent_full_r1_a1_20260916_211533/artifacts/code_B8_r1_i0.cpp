#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Eigen;

int main() {
    double coords[6][3] = {
        {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}
    };
    double E = 200e9, nu = 0.3;
    double GG = E/(2*(1+nu));

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double temp1 = (1.0+nu)*(1.0-2.0*nu);
    double g1 = E*(1.0-nu)/temp1;
    double g2 = E*nu/temp1;
    G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
    G(0,1)=g2; G(0,2)=g2;
    G(1,0)=g2; G(1,2)=g2;
    G(2,0)=g2; G(2,1)=g2;
    G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

    int M[12][4] = {
        {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
        {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
        {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
    };

    Matrix<double,18,18> K = Matrix<double,18,18>::Zero();

    for (int it = 0; it < 12; ++it) {
        Matrix4d H;
        for (int i = 0; i < 4; ++i) {
            int n = M[it][i] - 1;
            H(i,0) = 1.0;
            H(i,1) = coords[n][0];
            H(i,2) = coords[n][1];
            H(i,3) = coords[n][2];
        }
        double HDETER = H.determinant();
        Matrix4d Hinv = H.inverse();

        int IOPT = it + 1 + 10;
        double hd = std::abs(HDETER) / 36.0;
        if (IOPT <= 16) hd *= 2.0;

        Matrix<double,6,3> C[4];
        for (int i = 0; i < 4; ++i) {
            C[i].setZero();
            double b = Hinv(i,1);
            double c = Hinv(i,2);
            double d = Hinv(i,3);
            C[i](0,0) = b;
            C[i](4,0) = c;
            C[i](2,1) = d;
            C[i](4,1) = d;
            C[i](5,1) = c;
            C[i](0,2) = d;
            C[i](2,2) = b;
            C[i](3,2) = c;
            C[i](4,2) = b;
        }

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                Matrix3d Kij = hd * C[i].transpose() * G * C[j];
                int ni = M[it][i] - 1;
                int nj = M[it][j] - 1;
                K.block<3,3>(3*ni, 3*nj) += Kij;
            }
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; ++i) {
        std::cout << "[";
        for (int j = 0; j < 18; ++j) {
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if (j < 17) std::cout << ",";
        }
        std::cout << "]";
        if (i < 17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}