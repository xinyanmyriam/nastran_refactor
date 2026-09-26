#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Eigen;

double INVERD4(Matrix4d& A, int& ISING) {
    ISING = 0;
    Matrix4d M = A;
    double det = M.determinant();
    if (std::abs(det) < 1e-300) { ISING = 2; return det; }
    A = M.inverse();
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

    for (int itet = 0; itet < 12; itet++) {
        int IOPT = itet + 1 + 10;
        Matrix4d H;
        for (int j = 0; j < 4; j++) {
            int node = M[itet][j] - 1;
            H(j,0) = 1.0;
            H(j,1) = coords[node][0];
            H(j,2) = coords[node][1];
            H(j,3) = coords[node][2];
        }
        int ISING = 0;
        double HDETER = INVERD4(H, ISING);
        if (ISING == 2) continue;
        HDETER = std::abs(HDETER);

        MatrixXd G = MatrixXd::Zero(6,6);
        double TEMP1 = (1.0+nu)*(1.0-2.0*nu);
        double g1 = E*(1.0-nu)/TEMP1;
        double g2 = E*nu/TEMP1;
        G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
        G(0,1)=g2; G(0,2)=g2;
        G(1,0)=g2; G(1,2)=g2;
        G(2,0)=g2; G(2,1)=g2;
        G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

        HDETER = HDETER/36.0;
        if (IOPT <= 16) HDETER = HDETER*2.0;

        MatrixXd Ci[4];
        for (int i = 0; i < 4; i++) {
            Ci[i] = MatrixXd::Zero(6,3);
            double dNdx = H(i,1);
            double dNdy = H(i,2);
            double dNdz = H(i,3);
            Ci[i](0,0) = dNdx;
            Ci[i](4,0) = dNdy;
            Ci[i](2,1) = dNdz;
            Ci[i](4,1) = dNdz;
            Ci[i](5,1) = dNdy;
            Ci[i](0,2) = dNdz;
            Ci[i](2,2) = dNdx;
            Ci[i](3,2) = dNdy;
            Ci[i](4,2) = dNdx;
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Matrix3d Kij = Ci[i].transpose() * G * Ci[j] * HDETER;
                int ni = M[itet][i] - 1;
                int nj = M[itet][j] - 1;
                for (int a = 0; a < 3; a++)
                    for (int b = 0; b < 3; b++)
                        K(ni*3+a, nj*3+b) += Kij(a,b);
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