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

    Eigen::VectorXd u0(n), v0(n);
    u0 << 1, 0, 0;
    v0 << 0, 0, 0;

    double d2 = dt*dt;
    Eigen::MatrixXd Keff = M/d2 + K;
    Eigen::MatrixXd KeffInv = Keff.inverse();
    Eigen::MatrixXd A = 2.0*M/d2;
    Eigen::MatrixXd B = M/d2;

    Eigen::VectorXd um1 = u0 - dt * v0;
    Eigen::VectorXd un = u0;

    std::vector<Eigen::VectorXd> history;
    history.push_back(u0);

    for (int step = 1; step <= nsteps; ++step) {
        Eigen::VectorXd unp1 = KeffInv * (A * un - B * um1);
        um1 = un;
        un = unp1;
        if (step % 10 == 0) history.push_back(un);
    }

    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    for (size_t i = 0; i < history.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << history[i](j);
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
    return 0;
}