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
    // But we need to derive u1 from initial conditions first
    
    // Central difference formula for acceleration: 
    // a0 = (u1 - 2*u0 + u_{-1}) / dt^2
    // and velocity: v0 = (u1 - u_{-1}) / (2*dt)
    // Solving: u1 = u0 + dt*v0 + 0.5*dt^2*a0
    // where a0 = M^{-1} * (P0 - C*v0 - K*u0)
    // Since P0 = 0 (no external load), C = 0, so a0 = -M^{-1}*K*u0
    
    // Compute initial acceleration
    Vector a0 = -M.ldlt().solve(K * u0);  // a0 = -M^{-1} K u0
    
    // Compute u1 using Taylor expansion (central difference starting step)
    Vector u1 = u0 + dt * v0 + 0.5 * dt * dt * a0;
    
    // Now we have u0 and u1, can iterate with central difference
    // u_{n+1} = (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1}
    Matrix A = 2.0 * M / (dt * dt) - K;
    Matrix B = M / (dt * dt);
    
    // Store displacement history
    std::vector<std::vector<Real>> timeseries;
    
    // Store step 0
    std::vector<Real> step0 = {u0(0), u0(1), u0(2)};
    timeseries.push_back(step0);
    
    // Current and previous displacements
    Vector un_minus_1 = u0;
    Vector un = u1;
    
    // Run integration for 100 steps
    for (int step = 1; step <= nsteps; ++step) {
        // u_{n+1} = A * u_n - B * u_{n-1}
        Vector un_plus_1 = A * un - B * un_minus_1;
        
        // Update for next iteration
        un_minus_1 = un;
        un = un_plus_1;
        
        // Store every 10th step (steps 0, 10, 20, ..., 100)
        if (step % 10 == 0) {
            std::vector<Real> step_data = {un(0), un(1), un(2)};
            timeseries.push_back(step_data);
        }
    }
    
    // Output as JSON
    std::cout << "{\"test\":\"TRD1C\",\"timeseries\":[";
    
    for (size_t i = 0; i < timeseries.size(); ++i) {
        if (i > 0) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < 3; ++j) {
            if (j > 0) std::cout << ",";
            std::cout << std::scientific << std::setprecision(16) << timeseries[i][j];
        }
        std::cout << "]";
    }
    
    std::cout << "]}" << std::endl;
    
    return 0;
}