#include <iostream>
#include <cmath>
#include <iomanip>
int main(){
    double nodes[4][2] = {{0,0},{2,0},{2,1.5},{0,1.5}};
    double E=200e9, nu=0.3, t=0.01;
    int M[12] = {1,2,4, 2,3,1, 3,4,2, 4,1,3};
    double K[64];
    for(int i=0;i<64;i++) K[i]=0;
    for(int j=0;j<4;j++){
        int n0=M[3*j]-1, n1=M[3*j+1]-1, n2=M[3*j+2]-1;
        double x1=nodes[n0][0], y1=nodes[n0][1];
        double x2=nodes[n1][0], y2=nodes[n1][1];
        double x3=nodes[n2][0], y3=nodes[n2][1];
        double A = 0.5*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));
        double b1=y2-y3, b2=y3-y1, b3=y1-y2;
        double c1=x3-x2, c2=x1-x3, c3=x2-x1;
        double B[3][6] = {{b1,0,b2,0,b3,0},{0,c1,0,c2,0,c3},{c1,b1,c2,b2,c3,b3}};
        double inv2A = 1.0/(2.0*A);
        for(int i=0;i<3;i++)for(int k=0;k<6;k++) B[i][k]*=inv2A;
        double f = E/(1.0-nu*nu);
        double D[3][3] = {{f,f*nu,0},{f*nu,f,0},{0,0,f*(1.0-nu)/2.0}};
        double DB[3][6];
        for(int i=0;i<3;i++)for(int k=0;k<6;k++){DB[i][k]=0;for(int m=0;m<3;m++)DB[i][k]+=D[i][m]*B[m][k];}
        double kt[6][6];
        for(int i=0;i<6;i++)for(int k=0;k<6;k++){kt[i][k]=0;for(int m=0;m<3;m++)kt[i][k]+=B[m][i]*DB[m][k];}
        double scale = (t/2.0)*std::abs(A);
        for(int i=0;i<6;i++)for(int k=0;k<6;k++)kt[i][k]*=scale;
        int idx[3]={n0,n1,n2};
        for(int a=0;a<3;a++)for(int b=0;b<3;b++)
            for(int i=0;i<2;i++)for(int jj=0;jj<2;jj++)
                K[(2*idx[a]+i)*8 + (2*idx[b]+jj)] += kt[2*a+i][2*b+jj];
    }
    std::cout << "K66=" << K[6*8+6] << " K77=" << K[7*8+7] << std::endl;
    std::cout << "{\"stiffness_matrix\":[";
    for(int i=0;i<8;i++){
        std::cout << "[";
        for(int j=0;j<8;j++){
            std::cout << std::scientific << std::setprecision(6) << K[i*8+j];
            if(j<7) std::cout << ",";
        }
        std::cout << "]";
        if(i<7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}