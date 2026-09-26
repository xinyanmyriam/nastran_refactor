#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    Matrix<double,4,2> nodes;
    nodes << 0,0, 2,0, 2,1.5, 0,1.5;

    Matrix3d D;
    double c = E/(1.0-nu*nu);
    D << c, c*nu, 0,
         c*nu, c, 0,
         0, 0, c*(1.0-nu)/2.0;

    int M[12] = {1,2,4, 2,3,1, 3,4,2, 4,1,3};

    Matrix<double,8,8> K = Matrix<double,8,8>::Zero();

    for (int j=0;j<4;j++) {
        int n[3] = {M[3*j]-1, M[3*j+1]-1, M[3*j+2]-1};
        double x1=nodes(n[0],0), y1=nodes(n[0],1);
        double x2=nodes(n[1],0), y2=nodes(n[1],1);
        double x3=nodes(n[2],0), y3=nodes(n[2],1);
        double A2 = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
        double A = 0.5*A2;
        double b1=y2-y3, c1=x3-x2;
        double b2=y3-y1, c2=x1-x3;
        double b3=y1-y2, c3=x2-x1;
        Matrix<double,3,6> B;
        B << b1,0,b2,0,b3,0,
             0,c1,0,c2,0,c3,
             c1,b1,c2,b2,c3,b3;
        B /= (2.0*A);
        Matrix<double,6,6> Kt = t * A * B.transpose() * D * B;
        for (int i=0;i<3;i++)
            for (int k=0;k<3;k++)
                K.block<2,2>(2*n[i],2*n[k]) += Kt.block<2,2>(2*i,2*k);
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<8;i++) {
        std::cout << "[";
        for (int j=0;j<8;j++) {
            std::cout << K(i,j);
            if (j<7) std::cout << ",";
        }
        std::cout << "]";
        if (i<7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}