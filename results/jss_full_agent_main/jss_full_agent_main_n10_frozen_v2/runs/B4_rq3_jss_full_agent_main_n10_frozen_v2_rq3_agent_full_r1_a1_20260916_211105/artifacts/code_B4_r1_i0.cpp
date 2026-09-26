#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

int main() {
    double E = 200e9, nu = 0.3, t = 0.01;
    double x[4] = {0, 2, 2, 0};
    double y[4] = {0, 0, 1.5, 1.5};

    double c = E/(1.0-nu*nu);
    Matrix3d D;
    D << c, c*nu, 0,
         c*nu, c, 0,
         0, 0, c*(1.0-nu)/2.0;

    MatrixXd K = MatrixXd::Zero(8,8);

    // Map M (1-indexed): triangles
    int tri[4][3] = {{0,1,3},{1,2,0},{2,3,1},{3,0,2}};
    double th = t/2.0;

    for (int tr = 0; tr < 4; tr++) {
        int n1 = tri[tr][0], n2 = tri[tr][1], n3 = tri[tr][2];
        double x1=x[n1],y1=y[n1],x2=x[n2],y2=y[n2],x3=x[n3],y3=y[n3];
        double A2 = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
        double A = std::abs(A2)/2.0;
        MatrixXd B(3,6);
        double b1=y2-y3, b2=y3-y1, b3=y1-y2;
        double cc1=x3-x2, cc2=x1-x3, cc3=x2-x1;
        B << b1,0,b2,0,b3,0,
             0,cc1,0,cc2,0,cc3,
             cc1,b1,cc2,b2,cc3,b3;
        B /= A2;
        MatrixXd Ke = th*A*B.transpose()*D*B;
        int map[3] = {n1,n2,n3};
        for(int i=0;i<3;i++) for(int j=0;j<3;j++)
            for(int a=0;a<2;a++) for(int b=0;b<2;b++)
                K(2*map[i]+a, 2*map[j]+b) += Ke(2*i+a, 2*j+b);
    }

    std::cout << "{\"stiffness_matrix\":[";
    for(int i=0;i<8;i++){
        std::cout << "[";
        for(int j=0;j<8;j++){
            std::cout << std::scientific << std::setprecision(6) << K(i,j);
            if(j<7) std::cout << ",";
        }
        std::cout << "]";
        if(i<7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}