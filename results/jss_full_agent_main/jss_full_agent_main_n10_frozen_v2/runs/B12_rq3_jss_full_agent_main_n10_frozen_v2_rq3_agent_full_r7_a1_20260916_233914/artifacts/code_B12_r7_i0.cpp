#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace Eigen;

int main() {
    MatrixXd K(3,3);
    K << 2,-1,0,
        -1,2,-1,
         0,-1,2;
    MatrixXd M = MatrixXd::Identity(3,3);

    VectorXd u0(3); u0 << 1,0,0;
    VectorXd v0(3); v0 << 0,0,0;

    double dt = 0.01;
    int nsteps = 100;

    MatrixXd Keff = M / (dt*dt);
    MatrixXd Keff_inv = Keff.inverse();

    VectorXd um1 = u0 - dt*v0;

    std::vector<VectorXd> hist;
    hist.push_back(u0);

    VectorXd u_prev = um1;
    VectorXd u_cur  = u0;

    MatrixXd A = 2.0*M/(dt*dt) - K;
    MatrixXd B = M/(dt*dt);

    for (int n = 1; n <= nsteps; ++n) {
        VectorXd rhs = A*u_cur - B*u_prev;
        VectorXd u_next = Keff_inv * rhs;
        u_prev = u_cur;
        u_cur = u_next;
        hist.push_back(u_cur);
    }

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    for (int s = 0; s <= 100; s += 10) {
        if (s > 0) std::cout << ",";
        std::cout << "[";
        for (int i = 0; i < 3; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << hist[s](i);
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}