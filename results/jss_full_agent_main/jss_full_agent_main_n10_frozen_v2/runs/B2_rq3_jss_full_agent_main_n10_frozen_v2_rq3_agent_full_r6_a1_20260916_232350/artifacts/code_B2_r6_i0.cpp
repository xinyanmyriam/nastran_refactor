#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace Eigen;

int main() {
    double E = 200e9;
    double G = 76.923e9;
    double A = 0.01;
    double I1 = 8.333e-6;
    double I2 = 8.333e-6;
    double FJ = 1.667e-5;
    double K1 = 0.0, K2 = 0.0, I12 = 0.0;

    Vector3d PA(0,0,0), PB(2,0,0);
    Vector3d VECI = PA - PB;
    VECI = -VECI;
    double FL = VECI.norm();
    VECI = VECI / FL;

    Vector3d SMALV0(0,0,1);
    SMALV0 = SMALV0 / SMALV0.norm();

    Vector3d VECK = VECI.cross(SMALV0);
    VECK.normalize();
    Vector3d VECJ = VECK.cross(VECI);
    VECJ.normalize();

    double L = FL;
    double LSQ = L*L;