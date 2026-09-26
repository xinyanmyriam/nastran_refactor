#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace Eigen;

Matrix<double,6,6> tri_stiff(const Matrix<double,3,2>& xy, double E, double nu, double t) {
    Matrix3d D;
    double c = E/(1.0 - nu*nu);
    D << c, c*nu, 0, c*nu, c, 0, 0, 0, c*(1.0-nu)/2.0;
    double x1=xy(0,0),y1=xy(0,1),x2=xy(1,0),y2=xy(1,1),x3=xy(2,0),y3=xy(2,1);
    double A2 = (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
    double A = 0.5*std::abs(A2);
    double b1=y2-y3,b2=y3-y1,b3=y1-y2;
    double c1=x3-x2,c2=x1-x3,c3=x2-x1;
    Matrix<double,3,6> B;
    B.setZero();
    B(0,0)=b1;B(0,2)=b2;B(0,4)=b3;
    B(1,1)=c1;B(1,3)=c2;B(1,5)=c3;
    B(2,0)=c1;B(2,1)=b1;B(2,2)=c2;B(2,3)=b2;B(2,4)=c3;B(2,5)=b3;
    B/=(2.0*A);
    return t*A*B.transpose()*D*B;
}

int main(){
    Matrix<double,4,2> nodes;
    nodes<<0,0, 2,0, 2,1.5, 0,1.5;
    double E=200e9,nu=0.3,t=0.01;
    int M[12]={1,2,4, 2,3,1, 3,4,2, 4,1,3};
    Matrix<double,8,8> K=Matrix<double,8,8>::Zero();
    for(int j=0;j<4;++j){
        int n[3]={M[3*j]-1,M[3*j+1]-1,M[3*j+2]-1};
        Matrix<double,3,2> xy;
        for(int a=0;a<3;++a){xy(a,0)=nodes(n[a],0);xy(a,1)=nodes(n[a],1);}
        Matrix<double,6,6> Kt=tri_stiff(xy,E,nu,t/2.0);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            for(int r=0;r<2;++r)for(int c=0;c<2;++c)
                K(2*n[a]+r,2*n[b]+c)+=Kt(2*a+r,2*b+c);
    }
    std::cout<<std::setprecision(6)<<std::scientific;
    std::cout<<"{\"stiffness_matrix\":[";
    for(int i=0;i<8;++i){
        std::cout<<"[";
        for(int j=0;j<8;++j){std::cout<<K(i,j);if(j<7)std::cout<<",";}
        std::cout<<"]";if(i<7)std::cout<<",";
    }
    std::cout<<"]}"<<std::endl;
    return 0;
}