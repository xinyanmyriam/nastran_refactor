#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    std::cout << "HELLO_TEST_MARKER" << std::endl;
    double x1=0,y1=0,z1=0;
    double x2=1,y2=0,z2=0;
    double x3=0,y3=1,z3=0;
    double x4=0,y4=0,z4=1;

    double E  = 200e9;
    double nu = 0.3;
    double GG = E / (2.0*(1.0+nu));

    Matrix4d M;
    M << 1.0, x1, y1, z1,
         1.0, x2, y2, z2,
         1.0, x3, y3, z3,
         1.0, x4, y4, z4;

    double detM = M.determinant();
    Matrix4d Minv = M.inverse();

    Matrix<double,6,6> G = Matrix<double,6,6>::Zero();
    double TEMP1 = (1.0+nu)*(1.0-2.0*nu);
    double c1 = E*(1.0-nu)/TEMP1;
    double c2 = E*nu/TEMP1;
    G(0,0)=c1; G(1,1)=c1; G(2,2)=c1;
    G(1,0)=c2; G(2,0)=c2; G(0,1)=c2; G(2,1)=c2; G(0,2)=c2; G(1,2)=c2;
    G(3,3)=GG; G(4,4)=GG; G(5,5)=GG;

    Matrix<double,6,3> Cmat[4];
    for (int i=0;i<4;i++){
        double b = Minv(i,1);
        double c = Minv(i,2);
        double d = Minv(i,3);
        Matrix<double,6,3> C = Matrix<double,6,3>::Zero();
        C(0,0)=b;
        C(1,1)=c;
        C(2,2)=d;
        C(3,1)=d; C(3,2)=c;
        C(4,0)=d; C(4,2)=b;
        C(5,0)=c; C(5,1)=b;
        Cmat[i]=C;
    }

    double V = std::fabs(detM)/6.0;

    Matrix<double,12,12> K = Matrix<double,12,12>::Zero();
    for (int i=0;i<4;i++){
        for (int j=0;j<4;j++){
            Matrix3d Kij = Cmat[i].transpose() * G * Cmat[j] * V;
            K.block<3,3>(3*i,3*j) = Kij;
        }
    }

    std::cout << std::setprecision(10);
    std::cout << "{\"stiffness_matrix\":[";
    for (int r=0;r<12;r++){
        std::cout << "[";
        for (int c=0;c<12;c++){
            std::cout << std::scientific << K(r,c);
            if (c<11) std::cout << ",";
        }
        std::cout << "]";
        if (r<11) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}