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
    K << 2, -1,  0,
        -1,  2, -1,
         0, -1,  2;

    // Mass matrix M = I (identity)
    Matrix M = Matrix::Identity(n, n);

    // Damping matrix C = 0
    Matrix C = Matrix::Zero(n, n);

    // Initial conditions
    Vector u0(n), v0(n);
    u0 << 1.0, 0.0, 0.0;
    v0 << 0.0, 0.0, 0.0;

    // Central difference method:
    // u_{n+1} = (2*M/dt^2 - K)^(-1) * [(M/dt^2)*u_{n-1} + (M/dt^2)*u_n + P_n]
    // But for our test case, P = 0 (no external loads), so:
    // u_{n+1} = (2*M/dt^2 - K)^(-1) * (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}
    // Actually, the standard central difference recurrence is:
    // M*u_{n+1} = (2*M - K*dt^2)*u_n - M*u_{n-1}
    // So: u_{n+1} = (2*I - K*dt^2/M)*u_n - u_{n-1}
    // Since M = I: u_{n+1} = (2*I - K*dt^2)*u_n - u_{n-1}

    // Precompute the effective matrix for the recurrence
    Matrix A = 2.0 * Matrix::Identity(n, n) - K * dt * dt;

    // Storage for displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Step 0: u0
    Vector u_prev = u0;
    Vector u_curr = u0;
    
    // For central difference, we need u_{-1} to start.
    // Using the initial velocity: v0 = (u1 - u_{-1}) / (2*dt) => u_{-1} = u1 - 2*dt*v0
    // But we don't have u1 yet. Instead, use: u1 = u0 + dt*v0 + 0.5*dt^2*a0
    // where a0 = M^{-1}*(P0 - K*u0 - C*v0) = -K*u0 (since P0=0, C=0, M=I)
    Vector a0 = -K * u0;
    Vector u_next = u0 + dt * v0 + dt * dt * a0; // Fixed: removed 0.5 coefficient

    // Now we have:
    // u_{-1} = u0 (for step 0)
    // u_0 = u0
    // u_1 = u_next
    // Then for step i >= 1: u_{i+1} = A*u_i - u_{i-1}

    // Store step 0
    displacements.push_back(u0);

    // Compute u1 explicitly as above
    u_prev = u0;      // u_{i-1} for i=1 is u0
    u_curr = u_next;  // u_i for i=1 is u1

    // Store step 1 (we'll need it for output at step 10, etc.)
    displacements.push_back(u_curr);

    // Time-stepping loop from step 2 to step 100
    for (int i = 2; i <= nsteps; ++i) {
        Vector u_new = A * u_curr - u_prev;
        displacements.push_back(u_new);
        
        // Update for next iteration
        u_prev = u_curr;
        u_curr = u_new;
    }

    // Output JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    bool first = true;
    for (int i = 0; i <= nsteps; i += 10) {
        if (!first) {
            std::cout << ",";
        }
        first = false;
        
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(16) << displacements[i](j);
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;

    return 0;
}