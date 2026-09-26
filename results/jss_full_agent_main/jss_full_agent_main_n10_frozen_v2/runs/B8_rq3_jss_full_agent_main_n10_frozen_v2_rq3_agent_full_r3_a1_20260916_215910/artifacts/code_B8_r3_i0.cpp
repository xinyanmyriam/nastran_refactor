#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Eigen;

static const int M[12][4] = {
    {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
    {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
    {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
};

int main() {
    double coords[6][3] = {
        {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}
    };
    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0*(1.0+nu));

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0+nu)*(1.0-2.0*nu);
    double g1 = E*(1.0-nu)/TEMP1;
    double g2 = E*nu/TEMP1;
    G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
    G(0,1)=g2; G(0,2)=g2;
    G(1,0)=g2; G(1,2)=g2;
    G(2,0)=g2; G(2,1)=g2;
    G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

    Matrix<double,18,18> K = Matrix<double,18,18>::Zero();

    for (int itet = 0; itet < 12; itet++) {
        int nodes[4];
        for (int j = 0; j < 4; j++) nodes[j] = M[itet][j] - 1;

        Matrix4d H;
        for (int i = 0; i < 4; i++) {
            H(i,0) = 1.0;
            H(i,1) = coords[nodes[i]][0];
            H(i,2) = coords[nodes[i]][1];
            H(i,3) = coords[nodes[i]][2];
        }

        double det = H.determinant();
        if (std::abs(det) < 1e-300) continue;
        Matrix4d Hinv = H.inverse();
        double HDETER = std::abs(det);

        Matrix<double,6,3> Cmat[4];
        for (int I = 0; I < 4; I++) {
            Cmat[I] = Matrix<double,6,3>::Zero();
            double h_i4  = Hinv(I,1);
            double h_i8  = Hinv(I,2);
            double h_i12 = Hinv(I,3);
            Cmat[I](0,0) = h_i4;
            Cmat[I](4,0) = h_i8;
            Cmat[I](2,1) = h_i12;
            Cmat[I](4,1) = h_i12;
            Cmat[I](5,1) = h_i8;
            Cmat[I](0,2) = h_i12;
            Cmat[I](2,2) = h_i4;
            Cmat[I](3,2) = h_i8;
            Cmat[I](4,2) = h_i4;
        }

        int IOPT = itet + 11;
        HDETER = HDETER / 36.0;
        if (IOPT <= 16) HDETER = HDETER * 2.0;

        // Full 12x12 tetra stiffness: K_ab = C_a^T * G * C_b * HDETER
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                Matrix3d Tab = Cmat[a].transpose() * G * Cmat[b] * HDETER;
                int na = nodes[a], nb = nodes[b];
                for (int p = 0; p < 3; p++)
                    for (int q = 0; q < 3; q++)
                        K(3*na+p, 3*nb+q) += Tab(p,q);
            }
        }
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 18; i++) {
        std::cout << "[";
        for (int j = 0; j < 18; j++) {
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if (j < 17) std::cout << ",";
        }
        std::cout << "]";
        if (i < 17) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}