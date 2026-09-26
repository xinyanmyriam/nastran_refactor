#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    const int n = 3;
    const double dt = 0.01;
    const int nsteps = 100;

    Eigen::MatrixXd M = Eigen::MatrixXd::Identity(n, n);

    Eigen::MatrixXd K(n, n);
    K << 2, -1,  0,
        -1,  2, -1,
         0, -1,  2;

    Eigen::VectorXd u0(n);
    u0 << 1, 0, 0;
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(n);

    Eigen::MatrixXd Keff = M / (dt * dt);
    Eigen::MatrixXd Keff_inv = Keff.inverse();

    Eigen::VectorXd um1 = u0 - dt * v0;

    std::vector<Eigen::VectorXd> history;
    history.reserve(nsteps + 1);
    history.push_back(u0);

    Eigen::VectorXd u_prev = um1;
    Eigen::VectorXd u_curr = u0;

    Eigen::MatrixXd A = 2.0 * M / (dt * dt) - K;
    Eigen::MatrixXd B = M / (dt * dt);

    for (int step = 1; step <= nsteps; ++step) {
        Eigen::VectorXd rhs = A * u_curr - B * u_prev;
        Eigen::VectorXd u_next = Keff_inv * rhs;
        history.push_back(u_next);
        u_prev = u_curr;
        u_curr = u_next;
    }

    std::cout << std::setprecision(10);
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    bool first = true;
    for (int s = 0; s <= nsteps; s += 10) {
        if (!first) std::cout << ",";
        first = false;
        const Eigen::VectorXd& u = history[s];
        std::cout << "[" << std::scientific << u(0) << "," << u(1) << "," << u(2) << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}