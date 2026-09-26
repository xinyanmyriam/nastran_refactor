#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

using namespace Eigen;

int main() {
    double x1=0.0,y1=0.0,z1=0.0;
    double x2=2.0,y2=0.0,z2=0.0;
    double x3=1.0,y3=1.5,z3=0.0;
    double Emod=2.1e11, nu=0.3, t=0.01;

    Vector3d Ivec(x2-x1, y2-y1, z2-z1);
    double XSUBB = Ivec.norm();
    Ivec /= XSUBB;

    Vector3d RC(x3-x1, y3-y1, z3-z1);
    double XSUBC = Ivec.dot(RC);

    Vector3d Kvec = Ivec.cross(RC);
    double YSUBC = Kvec.norm();
    Kvec /= YSUBC;

    Vector3d Jvec = Kvec.cross(Ivec);
    double TEMP = Jvec.norm();
    if (TEMP != 0.0) Jvec /= TEMP;

    double VOL = XSUBB * YSUBC * t / 2.0;
    double REELMU = 1.0 / XSUBB;
    double FLAMDA = 1.0 / YSUBC;
    double DELTA = XSUBC / XSUBB - 1.0;

    MatrixXd C(3,6);
    C(0,0)=-REELMU;      C(0,1)=0.0;
    C(1,0)=0.0;          C(1,1)=FLAMDA*DELTA;
    C(2,0)=FLAMDA*DELTA; C(2,1)=-REELMU;
    C(0,2)=REELMU;       C(0,3)=0.0;
    C(1,2)=0.0;          C(1,3)=-FLAMDA*REELMU*XSUBC;
    C(2,2)=-FLAMDA*REELMU*XSUBC; C(2,3)=REELMU;
    C(0,4)=0.0;          C(0,5)=0.0;
    C(1,4)=0.0;          C(1,5)=FLAMDA;
    C(2,4)=FLAMDA;       C(2,5)=0.0;

    double G11 = Emod/(1.0-n