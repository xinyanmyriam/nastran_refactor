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
    // u_{n+1} = (2*M/dt^2 - K)^(-1) * [(M/dt^2)*u_{n-1} + (2*M/dt^2 - K)*u_n - (M/dt^2)*u_{n-1} + P_n]
    // But for undamped, unforced system: P_n = 0, so:
    // u_{n+1} = (2*M/dt^2 - K)^(-1) * [(2*M/dt^2 - K)*u_n - (M/dt^2)*u_{n-1}]
    //
    // However, we need u_{-1} to start the recurrence.
    // Using central difference approximation for velocity at t=0:
    // v0 = (u1 - u_{-1}) / (2*dt)  =>  u_{-1} = u1 - 2*dt*v0
    // But we don't know u1 yet. Instead, use the standard initialization:
    // u_{-1} = u0 - dt*v0 + 0.5*dt^2*a0, where a0 = M^{-1}*(P0 - K*u0 - C*v0)
    // Since P0 = 0, C = 0: a0 = -M^{-1}*K*u0 = -K*u0 (since M=I)
    Vector a0 = -K * u0;
    Vector um1 = u0 - dt * v0 + 0.5 * dt * dt * a0;

    // Precompute matrices for central difference
    Matrix M_dt2 = M / (dt * dt);
    Matrix A = 2.0 * M_dt2 - K;  // effective stiffness matrix

    // Solve for u1 using the equation:
    // M*(u1 - 2*u0 + um1)/dt^2 + K*u0 = 0
    // => M*u1 = 2*M*u0 - M*um1 - dt^2*K*u0
    // => u1 = 2*u0 - um1 - dt^2*M^{-1}*K*u0
    // Since M=I: u1 = 2*u0 - um1 - dt^2*K*u0
    Vector u1 = 2.0 * u0 - um1 - dt * dt * K * u0;

    // Store displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0
    displacements.push_back(u0);

    // Time integration loop
    Vector u_prev = um1;  // u_{n-1}
    Vector u_curr = u0;   // u_n
    Vector u_next(n);

    // We'll store every 10th step (0,10,20,...,100) -> 11 points
    std::vector<Vector> output_points;

    // Store step 0
    output_points.push_back(u0);

    // Integrate for nsteps steps
    for (int step = 1; step <= nsteps; ++step) {
        // Central difference formula: u_{n+1} = A^{-1} * [(2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}]
        // Since P=0 and C=0, and M=I: u_{n+1} = A^{-1} * [A * u_n - M_dt2 * u_{n-1}]
        // => u_{n+1} = u_n - A^{-1} * M_dt2 * u_{n-1}
        u_next = u_curr - A.ldlt().solve(M_dt2 * u_prev);

        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_next;

        // Store every 10th step (including step 100)
        if (step % 10 == 0) {
            output_points.push_back(u_curr);
        }
    }

    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    for (size_t i = 0; i < output_points.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(16) << output_points[i](j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;

    return 0;
}