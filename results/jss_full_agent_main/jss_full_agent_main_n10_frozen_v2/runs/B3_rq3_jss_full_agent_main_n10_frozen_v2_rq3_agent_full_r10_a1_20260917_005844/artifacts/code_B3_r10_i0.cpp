#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
using namespace Eigen;
int main() {
    double x1=0,y1=0,x2=2,y2=0,x3=1,y3=1.5;
    double E_mod=2.1e11, nu=0.3, t=0.01;
    double A = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
    MatrixXd B(3,6);
    B << (y2-y3),0,(y3-y1),0,(y1-y2),0,
         0,(x3-x2),0,(x1-x3),0,(x2-x1),
         (x3-x2),(y2-y3),(x1-x3),(y3-y1),(x2-x1),(y1-y2);
    B /= (2*A);
    Matrix3d D;
    D << 1,nu,0, nu,1,0, 0,0,(1-nu)/2;
    D *= E_mod/(1-nu*nu);
    MatrixXd K = t*A*B.transpose()*D*B;
    std::cout << "MARKER_XYZ\n";
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int r=0;r<6;r++){
        std::cout << "[";
        for (int c=0;c<6;c++){ std::cout << K(r,c); if(c<5) std::cout << ","; }
        std::cout << "]";
        if(r<5) std::cout << ",";
    }
    std::cout << "]}\n";
    return 0;
}