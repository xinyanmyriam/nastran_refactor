#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <sstream>
#include <string>

// Use double precision throughout
using Real = double;
using Vector = Eigen::VectorXd;
using Matrix = Eigen::MatrixXd;

int main() {
    // Test case parameters
    const int n = 3;  // 3-DOF system
    const Real dt = 0.01;
    const int nsteps = 100;

    // Stiffness matrix K = [2,-1,0; -1,2,-1; 0,-1,2]
    Matrix K(n, n);
    K << 2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0;

    // Mass matrix M = I (identity)
    Matrix M = Matrix::Identity(n, n);

    // Damping matrix C = 0
    Matrix C = Matrix::Zero(n, n);

    // Initial conditions
    Vector u0(n);
    u0 << 1.0, 0.0, 0.0;
    
    Vector v0(n);
    v0 << 0.0, 0.0, 0.0;

    // Central difference method:
    // u_{n+1} = 2*u_n - u_{n-1} + dt^2 * a_n
    // where a_n = M^{-1}*(P_n - C*v_n - K*u_n) = -K*u_n (since M=I, C=0, P=0)
    // So: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n

    // Compute u_{-1} using u0 and v0:
    // u_{-1} = u0 - dt*v0 - 0.5*dt^2*a0
    // a0 = -K*u0
    Vector a0 = -K * u0;
    Vector um1 = u0 - dt * v0 - 0.5 * dt * dt * a0;

    // Storage for displacement history (every 10th step including 0 and 100)
    std::vector<std::vector<Real>> timeseries;
    timeseries.reserve(11);

    // Initialize state vectors
    Vector u_prev = um1;   // u_{n-1}
    Vector u_curr = u0;    // u_n
    Vector u_next(n);      // u_{n+1}

    // Store step 0
    std::vector<Real> step0;
    step0.reserve(n);
    for (int i = 0; i < n; ++i) {
        step0.push_back(u_curr(i));
    }
    timeseries.push_back(step0);

    // Perform time integration
    for (int step = 1; step <= nsteps; ++step) {
        // Central difference: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n
        u_next = 2.0 * u_curr - u_prev - dt * dt * K * u_curr;

        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;

        // Store every 10th step (0, 10, 20, ..., 100)
        if (step % 10 == 0) {
            std::vector<Real> step_data;
            step_data.reserve(n);
            for (int i = 0; i < n; ++i) {
                step_data.push_back(u_curr(i));
            }
            timeseries.push_back(step_data);
        }
    }

    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    for (size_t i = 0; i < timeseries.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (size_t j = 0; j < timeseries[i].size(); ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(16) << timeseries[i][j];
        }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;

    return 0;
}