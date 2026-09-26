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

    // Precompute dt^2 * K
    Matrix dt2K = dt * dt * K;

    // Storage for displacement history
    std::vector<Vector> displacements;
    displacements.reserve(nsteps + 1);

    // Initial state
    Vector u_prev = u0;
    Vector u_curr = u0;

    // Store step 0
    displacements.push_back(u_curr);

    // Compute u_1 using initial velocity v0 and acceleration a0
    // From equation of motion: M*a0 + C*v0 + K*u0 = 0 => a0 = -M^{-1}*K*u0 = -K*u0 (since M=I, C=0)
    Vector a0 = -K * u0;
    
    // Proper second-order initialization for central difference:
    // u_{-1} = u0 - dt*v0 - (dt^2/2)*a0
    Vector u_minus1 = u0 - dt * v0 - (dt * dt / 2.0) * a0;

    // Now run the integration for nsteps steps
    // We need to store u at steps: 0, 10, 20, ..., 100 (11 points total)
    std::vector<Vector> output_points;
    output_points.reserve(11);

    // Step 0
    output_points.push_back(u_curr);

    Vector u_n_minus1 = u_minus1;
    Vector u_n = u_curr;

    // Perform 100 time steps
    for (int step = 1; step <= nsteps; ++step) {
        // u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n
        Vector u_n_plus1 = 2.0 * u_n - u_n_minus1 - dt2K * u_n;

        // Update for next iteration
        u_n_minus1 = u_n;
        u_n = u_n_plus1;

        // Store every 10th step (including step 100)
        if (step % 10 == 0) {
            output_points.push_back(u_n);
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