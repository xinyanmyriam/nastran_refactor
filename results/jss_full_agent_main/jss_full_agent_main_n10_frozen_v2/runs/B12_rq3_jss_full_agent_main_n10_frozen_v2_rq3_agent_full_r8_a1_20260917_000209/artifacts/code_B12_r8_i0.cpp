#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace Eigen;

int main() {
    Matrix3d K;
    K << 2, -1, 0,
        -1, 2, -1,
         0, -1, 2;
    Matrix3d M = Matrix3d::Identity();
    Matrix3d C = Matrix3d::Zero();

    Vector3d u0;
    u0 << 1, 0, 0;
    Vector3d v0 = Vector3d::Zero();

    double dt = 0.01;
    int nsteps = 100;

    Matrix3d Keff = M / (dt * dt);
    Matrix3d KeffInv = Keff.inverse();

    Vector3d um1 = u0 - dt * v0;
    Vector3d un = u0;

    std::vector<Vector3d> history;
    history.push_back(u0);

    Matrix3d A = 2.0 * M / (dt * dt) - K;
    Matrix3d B = M / (dt * dt);

    Vector3d up1;
    for (int n = 0; n < nsteps; ++n) {
        up1 = KeffInv * (A * un - B * um1);
        um1 = un;
        un = up1;
        if ((n + 1) % 10 == 0) {
            history.push_back(un);
        }
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    for (size_t i = 0; i < history.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[" << history[i](0) << "," << history[i](1) << "," << history[i](2) << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}