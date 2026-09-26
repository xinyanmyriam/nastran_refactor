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
    // For central difference: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * M^{-1} * K * u_n
    // Since M = I, this simplifies to: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n
    
    // Compute initial acceleration
    Vector a0 = -M.ldlt().solve(K * u0);  // a0 = -M^{-1} K u0
    
    // Compute u1 using Taylor expansion (central difference starting step)
    Vector u1 = u0 + dt * v0 + 0.5 * dt * dt * a0;
    
    // Store displacement history - we need exactly 11 values: steps 0, 10, 20, ..., 100
    std::vector<std::vector<Real>> timeseries;
    
    // Store step 0
    std::vector<Real> step0 = {u0(0), u0(1), u0(2)};
    timeseries.push_back(step0);
    
    // Current and previous displacements
    Vector un_minus_1 = u0;
    Vector un = u1;
    
    // Run integration for 99 steps to get from u1 to u100 (since we need u100)
    // We already have step 0, so we need to compute u10, u20, ..., u100
    // After k iterations, un = u_{k+1}, so to get u10 we need k=9, u20 needs k=19, etc.
    for (int step = 1; step <= 99; ++step) {
        // Central difference formula: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n
        Vector un_plus_1 = 2.0 * un - un_minus_1 - dt * dt * K * un;
        
        // Update for next iteration
        un_minus_1 = un;
        un = un_plus_1;
        
        // Store this step if it's one of the required reporting steps (10,20,...,100)
        // After step iterations, un = u_{step+1}, so we want step+1 = 10,20,...,100
        // i.e., step = 9,19,29,...,99
        if (step % 10 == 9) {
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