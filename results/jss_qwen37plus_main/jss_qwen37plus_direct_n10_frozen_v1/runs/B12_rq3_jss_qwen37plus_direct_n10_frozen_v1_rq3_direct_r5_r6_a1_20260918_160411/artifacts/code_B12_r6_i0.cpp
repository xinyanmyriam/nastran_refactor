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
    Vector u0(n), v0(n);
    u0 << 1.0, 0.0, 0.0;
    v0 << 0.0, 0.0, 0.0;

    // Precompute constants for central difference method
    // For central difference: u_{n+1} = (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}
    // Since M = I: u_{n+1} = (2/dt^2 - K) * u_n - (1/dt^2) * u_{n-1}
    Real inv_dt2 = 1.0 / (dt * dt);
    Matrix A = 2.0 * inv_dt2 * Matrix::Identity(n, n) - K;
    Real coeff = inv_dt2;

    // Storage for displacement history (we need steps 0,10,20,...,100 -> 11 points)
    std::vector<std::vector<Real>> timeseries;
    timeseries.reserve(11);

    // Initialize displacement vectors
    Vector u_prev = u0;  // u_{n-1}
    Vector u_curr = u0;  // u_n
    Vector u_next(n);    // u_{n+1}

    // Step 0: store initial state
    std::vector<Real> step0 = {u_curr(0), u_curr(1), u_curr(2)};
    timeseries.push_back(step0);

    // Compute u_1 using initial velocity: u_1 = u_0 + dt * v_0 + 0.5 * dt^2 * a_0
    // But for central difference, we need u_{-1} to start the recurrence.
    // From central difference: v_0 = (u_1 - u_{-1}) / (2*dt) => u_{-1} = u_1 - 2*dt*v_0
    // And u_1 = u_0 + dt*v_0 + 0.5*dt^2*a_0, where a_0 = M^{-1}*(P_0 - C*v_0 - K*u_0)
    // Here P_0 = 0 (no external load), C=0, so a_0 = -K*u_0
    Vector a0 = -K * u0;  // acceleration at t=0
    Vector u1 = u0 + dt * v0 + 0.5 * dt * dt * a0;
    
    // Now we have u0 and u1, so u_{-1} = 2*u0 - u1 (from central difference formula rearranged)
    // Actually: v0 = (u1 - u_{-1})/(2*dt) => u_{-1} = u1 - 2*dt*v0 = u1 since v0=0
    // But let's use the standard approach: for central difference, the first step is computed separately
    // Then the recurrence starts from n=1: u_{2} = A*u_{1} - coeff*u_{0}
    
    // So set:
    // u_prev = u0 (will be u_{n-1} for n=1)
    // u_curr = u1 (will be u_n for n=1)
    u_prev = u0;
    u_curr = u1;

    // Perform time integration for 100 steps
    for (int step = 1; step <= nsteps; ++step) {
        // Central difference recurrence:
        // u_{n+1} = A * u_n - coeff * u_{n-1}
        u_next = A * u_curr - coeff * u_prev;

        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;

        // Store every 10th step (0,10,20,...,100)
        if (step % 10 == 0) {
            std::vector<Real> step_data = {u_curr(0), u_curr(1), u_curr(2)};
            timeseries.push_back(step_data);
        }
    }

    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    for (size_t i = 0; i < timeseries.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 3; ++j) {
            std::cout << std::scientific << std::setprecision(16) << timeseries[i][j];
            if (j < 2) std::cout << ",";
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;

    return 0;
}