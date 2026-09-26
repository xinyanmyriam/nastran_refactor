#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double x1=0,y1=0,z1=0;
    double x2=1,y2=0,z2=0;
    double x3=0,y3=1,z3=0;
    double x4=0,y4=0,z4=1;

    double E = 200e9;
    double nu = 0.3;
    double GG = E / (2.0*(1.0+nu));

    Matrix4d H;
    H << 1, x1, y1, z1,
         1, x2, y2, z2,
         1, x3, y3, z3,
         1, x4, y4, z4;

    double HDETER = H.determinant();
    Matrix4d Hinv = H.inverse();
    HDETER = std::fabs(HDETER);

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0+nu)*(1.0-2.0*nu);
    double g1 = E*(1.0-nu)/TEMP1;
    double g2 = E*nu/TEMP1;
    G(0,0)=g1; G(1,1)=g1; G(2,2)=g1;
    G(0,1)=g2; G(0,2)=g2;
    G(1,0)=g2; G(1,2)=g2;
    G(2,0)=g2; G(2,1)=g2;
    G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

    double Hv[17];
    Hv[1]=1; Hv[2]=x1; Hv[3]=y1; Hv[4]=z1;
    Hv[5]=1; Hv[6]=x2; Hv[7]=y2; Hv[8]=z2;
    Hv[9]=1; Hv[10]=x3; Hv[11]=y3; Hv[12]=z3;
    Hv[13]=1; Hv[14]=x4; Hv[15]=y4; Hv[16]=z4;

    Matrix<double,6,3> Cblk[4];
    for (int I = 1; I <= 4; I++) {
        double bx = Hinv(I-1,1);
        double by = Hinv(I-1,2);
        double bz = Hinv(I-1,3);
        Matrix<double,6,3> C = Matrix<double,6,3>::Zero();
        C(0,0) = bx;
        C(4,0) = by;
        C(2,1) = bz;
        C(4,1) = bz;
        C(0,2) = by;
        C(1,2) = bz;
        C(2,2) = bx;
        C(3,2) = by;
        C(4,2) = bx;
        Cblk[I-1] = C;
    }

    HDETER = HDETER / 6.0;

    int pivotNode = 1;
    int CpivIdx = pivotNode - 1;

    Matrix<double,3,6> GCT = Cblk[CpivIdx].transpose() * G * HDETER;

    Matrix<double,12,12> K = Matrix<double,12,12>::Zero();

    for (int I = 1; I <= 4; I++) {
        Matrix<double,3,3> T = GCT * Cblk[I-1];
        int r0 = (pivotNode-1)*3;
        int c0 = (I-1)*3;
        K.block<3,3>(r0, c0) = T;
    }

    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 12; i++) {
        std::cout << "[";
        for (int j = 0; j < 12; j++) {
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if (j < 11) std::cout << ",";
        }
        std::cout << "]";
        if (i < 11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}