#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double coords[6][3] = {
        {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}
    };

    int M[12][4] = {
        {1,2,3,4},{1,2,3,5},{1,2,3,6},{1,4,5,6},
        {2,4,5,6},{3,4,5,6},{2,1,4,6},{2,3,4,6},
        {1,3,4,5},{2,3,4,5},{3,1,5,6},{2,1,5,6}
    };

    MatrixXd K = MatrixXd::Zero(18,18);
    int NPVT = 1;

    for (int itet = 0; itet < 12; itet++) {
        std::vector<double> tecpt(23, 0.0);
        tecpt[0] = 1.0; tecpt[1] = 1.0;
        for (int j = 0; j < 4; j++) tecpt[2+j] = M[itet][j];
        for (int j = 0; j < 4; j++) {
            int kp = M[itet][j];
            tecpt[6 + 4*j] = 0.0;
            tecpt[7 + 4*j] = coords[kp-1][0];
            tecpt[8 + 4*j] = coords[kp-1][1];
            tecpt[9 + 4*j] = coords[kp-1][2];
        }
        tecpt[22] = 0.0;

        int nodeMap[4];
        for (int j = 0; j < 4; j++) nodeMap[j] = M[itet][j]-1;

        Matrix4d H;
        H(0,0)=1.0; H(0,1)=tecpt[7]; H(0,2)=tecpt[8]; H(0,3)=tecpt[9];
        H(1,0)=1.0; H(1,1)=tecpt[11]; H(1,2)=tecpt[12]; H(1,3)=tecpt[13];
        H(2,0)=1.0; H(2,1)=tecpt[15]; H(2,2)=tecpt[16]; H(2,3)=tecpt[17];
        H(3,0)=1.0; H(3,1)=tecpt[19]; H(3,2)=tecpt[20]; H(3,3)=tecpt[21];

        double hdeter = H.determinant();
        if (std::abs(hdeter) < 1e-300) continue;
        Matrix4d Hinv = H.inverse();
        hdeter = std::abs(hdeter);

        double E = 200e9, nu = 0.3;
        double GG = E/(2.0*(1.0+nu));
        double temp1 = (1.0+nu)*(1.0-2.0*nu);
        double g1 = E*(1.0-nu)/temp1;
        double g2 = E*nu/temp1;
        Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
        G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
        G(0,1)=g2; G(0,2)=g2;
        G(1,0)=g2; G(1,2)=g2;
        G(2,0)=g2; G(2,1)=g2;
        G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

        std::vector<Matrix<double,6,3>> C(4);
        for (int i = 0; i < 4; i++) {
            C[i] = Matrix<double,6,3>::Zero();
            double h4 = Hinv(i,0);
            double h8 = Hinv(i,1);
            double h12= Hinv(i,2);
            C[i](0,0) = h4;
            C[i](4,0) = h8;
            C[i](2,1) = h12;
            C[i](4,1) = h12;
            C[i](5,1) = h8;
            C[i](0,2) = h12;
            C[i](2,2) = h4;
            C[i](3,2) = h8;
            C[i](4,2) = h4;
        }

        double hd;
        int iopt = (itet+1) + 10;
        if (iopt >= 11 && iopt <= 22) {
            hd = hdeter/36.0;
            if (iopt <= 16) hd = hd*2.0;
        } else {
            hd = hdeter/6.0;
        }

        int pivotLocal = -1;
        for (int i = 0; i < 4; i++) {
            if (nodeMap[i]+1 == NPVT) { pivotLocal = i; break; }
        }
        if (pivotLocal < 0) continue;

        Matrix<double,3,6> GCT = C[pivotLocal].transpose() * G;
        GCT *= hd;

        for (int i = 0; i < 4; i++) {
            Matrix3d KIJ = GCT * C[i];
            int ni = nodeMap[i];
            int np = nodeMap[pivotLocal];
            for (int a = 0; a < 3; a++)
                for (int b = 0; b < 3; b++)
                    K(3*ni+a, 3*np+b) += KIJ(a,b);
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