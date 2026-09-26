#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace Eigen;

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    MatrixXd coords(4,2);
    coords << 0,0, 2,0, 2,1.5, 0,1.5;
    
    Matrix3d D;
    double c = E/(1.0-nu*nu);
    D << c, c*nu, 0,
         c*nu, c, 0,
         0, 0, c*(1.0-nu)/2.0;
    
    int M[12] = {0,1,3, 1,2,0, 2,3,1, 3,0,2};
    
    MatrixXd K8 = MatrixXd::Zero(8,8);
    
    for (int j=0;j<4;j++) {
        int nd[3] = {M[3*j], M[3*j+1], M[3*j+2]};
        double x1=coords(nd[0],0), y1=coords(nd[0],1);
        double x2=coords(nd[1],0), y2=coords(nd[1],1);
        double x3=coords(nd[2],0), y3=coords(nd[2],1);
        double A2 = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
        double A = 0.5*std::abs(A2);
        double b1=y2-y3, b2=y3-y1, b3=y1-y2;
        double c1=x3-x2, c2=x1-x3, c3=x2-x1;
        MatrixXd B(3,6);
        B << b1,0,b2,0,b3,0,
             0,c1,0,c2,0,c3,
             c1,b1,c2,b2,c3,b3;
        B /= A2;
        MatrixXd K = t*A*B.transpose()*D*B;
        for (int a=0;a<3;a++)
            for (int b=0;b<3;b++)
                for (int r=0;r<2;r++)
                    for (int s=0;s<2;s++)
                        K8(2*nd[a]+r, 2*nd[b]+s) += K(2*a+r, 2*b+s);
    }
    
    std::cout << "{\"stiffness_matrix\":[";
    for (int i=0;i<8;i++) {
        std::cout << "[";
        for (int j=0;j<8;j++) {
            std::cout << std::scientific << std::setprecision(6) << K8(i,j);
            if (j<7) std::cout << ",";
        }
        std::cout << "]";
        if (i<7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}