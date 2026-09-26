#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace Eigen;

int main() {
    int n = 3;
    MatrixXd K(n, n);
    K << 2, -1, 0,
        -1, 2, -1,
         0, -1, 2;
    MatrixXd M = MatrixXd::Identity(n, n);

    VectorXd u0(n);
    u0 << 1, 0, 0;
    VectorXd v0 = VectorXd::Zero(n);

    double dt = 0.01;
    int nsteps = 100;

    // Central difference (explicit):
    //   K_eff   = M/dt^2
    //   u_{n+1} = K_eff^-1 * [ (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} ]
    MatrixXd Keff = M / (dt * dt);
    MatrixXd Keff_inv = Keff.inverse();
    MatrixXd A = 2.0 * M / (dt * dt) - K;
    MatrixXd B = M / (dt * dt);

    // Second-order starting: u_{-1} = u0 - dt*v0 + (dt^2/2)*a0
    // a0 = M^{-1}(P0 - K*u0 - C*v0) = -K*u0  (M=I, P0=0, C=0)
    VectorXd a0 = -K * u0;
    VectorXd u_prev = u0 - dt * v0 + 0.5 * dt * dt * a0;  // u_{-1}
    VectorXd u_curr = u0;                                  // u_0

    std::vector<VectorXd> history;
    history.push_back(u_curr);  // step 0

    for (int step = 1; step <= nsteps; ++step) {
        VectorXd u_next = Keff_inv * (A * u_curr - B * u_prev);
        u_prev = u_curr;
        u_curr = u_next;
        if (step % 10 == 0) {
            history.push_back(u_curr);
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