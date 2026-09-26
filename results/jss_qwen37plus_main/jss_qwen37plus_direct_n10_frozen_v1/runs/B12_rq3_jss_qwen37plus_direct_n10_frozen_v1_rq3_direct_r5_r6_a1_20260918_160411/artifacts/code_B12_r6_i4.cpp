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
    // For central difference: u_{n+1} = 2*u_n - dt^2*K*u_n - u_{n-1}
    // Since M = I: u_{n+1} = (2*I - dt^2*K) * u_n - u_{n-1}
    Real dt2 = dt * dt;
    Matrix A = 2.0 * Matrix::Identity(n, n) - dt2 * K;
    Real coeff = 1.0;

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
    // a_0 = M^{-1}*(P_0 - C*v_0 - K*u_0) = -K*u_0 (since P_0=0, C=0, M=I)
    Vector a0 = -K * u0;  // acceleration at t=0
    Vector u1 = u0 + dt * v0 + 0.5 * dt2 * a0;
    
    // Set up for central difference recurrence
    // u_prev should be u0 (u_{n-1} for n=1)
    // u_curr should be u1 (u_n for n=1)
    u_prev = u0;
    u_curr = u1;

    // Perform time integration for 100 steps
    // After 'step' iterations, u_curr contains u_{step+1}
    // We want to store u10, u20, ..., u100, which correspond to step=9,19,...,99
    for (int step = 1; step <= nsteps; ++step) {
        // Central difference recurrence:
        // u_{n+1} = A * u_n - coeff * u_{n-1}
        u_next = A * u_curr - coeff * u_prev;

        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;

        // Store every 10th step (0,10,20,...,100)
        // After step iterations, u_curr = u_{step+1}
        // So for u10, we need step+1 = 10 => step = 9
        // For u100, we need step+1 = 100 => step = 99
        if (step % 10 == 9) {
            // Store u10, u20, ..., u100 (step=9,19,...,99)
            std::vector<Real> step_data = {u_curr(0), u_curr(1), u_curr(2)};
            timeseries.push_back(step_data);
        }
    }

    // Ensure we have exactly 11 points: u0, u10, u20, ..., u100
    // Currently we have u0 + 10 points from the loop = 11 points
    // But verify size and truncate if needed to exactly 11
    if (timeseries.size() > 11) {
        timeseries.resize(11);
    } else if (timeseries.size() < 11) {
        // This shouldn't happen, but pad with zeros if necessary
        while (timeseries.size() < 11) {
            timeseries.push_back({0.0, 0.0, 0.0});
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