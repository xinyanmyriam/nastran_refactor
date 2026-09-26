#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>

int main() {
    Eigen::Matrix3d K;
    K << 2, -1, 0,
        -1, 2, -1,
         0, -1, 2;
    Eigen::Matrix3d M = Eigen::Matrix3d::Identity();
    double dt = 0.01;
    int nsteps = 100;
    Eigen::Vector3d u0(1.0, 0.0, 0.0);
    Eigen::Vector3d v0(0.0, 0.0, 0.0);
    Eigen::Matrix3d Keff = M / (dt * dt);
    Eigen::Matrix3d KeffInv = Keff.inverse();
    Eigen::Matrix3d A = 2.0 * M / (dt * dt) - K;
    Eigen::Matrix3d B = M / (dt * dt);
    Eigen::Vector3d a0 = M.inverse() * (-K * u0);
    Eigen::Vector3d u1 = u0 + dt * v0 + 0.5 * dt * dt * a0;
    std::vector<Eigen::Vector3d> history;
    history.push_back(u0);
    history.push_back(u1);
    Eigen::Vector3d u_prev = u0, u_curr = u1;
    for (int n = 1; n < nsteps; ++n) {
        Eigen::Vector3d u_next = KeffInv * (A * u_curr - B * u_prev);
        u_prev = u_curr; u_curr = u_next;
        history.push_back(u_curr);
    }
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    bool first = true;
    for (int s = 0; s <= nsteps; s += 10) {
        if (!first) std::cout << ",";
        first = false;
        const Eigen::Vector3d& u = history[s];
        std::cout << "[" << u(0) << "," << u(1) << "," << u(2) << "]";
    }
    std::cout << "]}" << std::endl;
    return 0;
}