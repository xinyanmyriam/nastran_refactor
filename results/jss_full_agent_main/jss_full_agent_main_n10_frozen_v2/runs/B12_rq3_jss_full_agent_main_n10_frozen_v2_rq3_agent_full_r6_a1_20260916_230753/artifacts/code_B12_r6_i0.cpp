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

    std::vector<Vector3d> history;
    history.push_back(u0);

    Vector3d un = u0;
    Vector3d unm1 = um1;

    for (int n = 0; n < nsteps; ++n) {
        Vector3d rhs = (2.0 * M / (dt * dt) - K) * un - (M / (dt * dt)) * unm1;
        Vector3d unp1 = KeffInv * rhs;

        unm1 = un;
        un = unp1;

        if ((n + 1) % 10 == 0) {
            history.push_back(un);
        }
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    for (size_t i = 0; i < history.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[" << history[i](0) << "," << history[i](1) << "," << history[i](2) << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}