#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    const int n = 3;
    const double dt = 0.01;
    const int nsteps = 100;

    Eigen::MatrixXd K(n, n);
    K << 2, -1, 0,
        -1, 2, -1,
         0, -1, 2;
    Eigen::MatrixXd M = Eigen::MatrixXd::Identity(n, n);

    Eigen::VectorXd u0(n);
    u0 << 1, 0, 0;
    Eigen::VectorXd v0 = Eigen::VectorXd::Zero(n);

    Eigen::MatrixXd Keff = M / (dt * dt);
    Eigen::MatrixXd Keff_inv = Keff.inverse();

    Eigen::VectorXd um1 = u0 - dt * v0;

    std::vector<Eigen::VectorXd> history;
    history.push_back(u0);

    Eigen::VectorXd u_prev = um1;
    Eigen::VectorXd u_cur = u0;

    Eigen::MatrixXd A = 2.0 * M / (dt * dt) - K;
    Eigen::MatrixXd B = M / (dt * dt);

    for (int step = 1; step <= nsteps; ++step) {
        Eigen::VectorXd rhs = A * u_cur - B * u_prev;
        Eigen::VectorXd u_next = Keff_inv * rhs;
        u_prev = u_cur;
        u_cur = u_next;
        if (step % 10 == 0) {
            history.push_back(u_cur);
        }
    }

    std::cout << std::setprecision(10);
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    for (size_t i = 0; i < history.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << history[i](j);
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}