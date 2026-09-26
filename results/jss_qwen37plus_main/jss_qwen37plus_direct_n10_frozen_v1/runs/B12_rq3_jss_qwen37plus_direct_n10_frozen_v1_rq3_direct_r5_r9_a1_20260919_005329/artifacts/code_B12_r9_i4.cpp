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
    // u_{n+1} = M^{-1} * [ (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} + P_n ]
    // But since C=0 and no external loads (P=0), we have:
    // u_{n+1} = (2*M/dt^2 - K)^{-1} * (M/dt^2) * (2*u_n - u_{n-1})
    // Actually, the standard central difference for M*u'' + K*u = 0 gives:
    // M*(u_{n+1} - 2*u_n + u_{n-1})/dt^2 + K*u_n = 0
    // => u_{n+1} = 2*u_n - u_{n-1} - dt^2 * M^{-1} * K * u_n
    // Since M = I: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n
    // => u_{n+1} = (2*I - dt^2*K) * u_n - u_{n-1}

    // However, the problem description states:
    //   K_eff = M/dt^2
    //   u_{n+1} = K_eff^-1 * [ (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} ]
    // With M=I: K_eff = I/dt^2, so K_eff^-1 = dt^2 * I
    // Then: u_{n+1} = dt^2 * [ (2*I/dt^2 - K) * u_n - (I/dt^2) * u_{n-1} ]
    //        = dt^2 * (2*I/dt^2 * u_n) - dt^2 * K * u_n - dt^2 * (I/dt^2) * u_{n-1}
    //        = 2*u_n - dt^2 * K * u_n - u_{n-1}
    // Which matches the standard form.

    // So we'll use: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n

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

    // Compute u_1 using initial velocity v0
    // From central difference: v0 = (u_1 - u_{-1}) / (2*dt)
    // and u_0 = (u_1 + u_{-1}) / 2 => u_{-1} = 2*u_0 - u_1
    // So: v0 = (u_1 - (2*u_0 - u_1)) / (2*dt) = (2*u_1 - 2*u_0) / (2*dt) = (u_1 - u_0) / dt
    // Therefore: u_1 = u_0 + dt * v_0
    Vector u_next = u_curr + dt * v0;

    // For the central difference scheme, we need u_{n-1}, u_n, then compute u_{n+1}
    // So at step 1, we have:
    //   u_{-1} = u_prev
    //   u_0   = u_curr
    //   u_1   = u_next (computed above)
    // Then for step 2 onwards: u_{n+1} = 2*u_n - u_{n-1} - dt^2 * K * u_n

    // But wait - the problem says "Derive the first step from the given initial conditions"
    // And the reference was produced by running TRD1C. Let's check what TRD1C does.

    // Looking at FORM1 subroutine: 
    //   U1(I) = U0(I) - DELTT*UDOT0(I)   [so u_{-1} = u0 - dt*v0]
    //   Then the integration starts with u_{-1}, u0, and computes u1, u2, etc.

    // So the state variables in TRD1C are:
    //   IU1 -> u_{n-1}
    //   IU2 -> u_n
    //   IU3 -> u_{n+1} (to be computed)
    // And initially: 
    //   IU1 = u_{-1} = u0 - dt*v0
    //   IU2 = u0
    //   Then u1 is computed as: u1 = 2*u0 - u_{-1} - dt^2*K*u0 = 2*u0 - (u0 - dt*v0) - dt^2*K*u0 = u0 + dt*v0 - dt^2*K*u0

    // But the problem says "Initial: u0 = [1,0,0], v0 = [0,0,0]", so u_{-1} = u0.
    // Then u1 = 2*u0 - u0 - dt^2*K*u0 = u0 - dt^2*K*u0.

    // However, the problem also says "the reference was produced by running it", and we must match that.
    // Let's follow the exact algorithm described in the problem:
    //   K_eff = M/dt^2 = I/dt^2
    //   u_{n+1} = K_eff^-1 * [ (2*M/dt^2 - K) * u_n - (M/dt^2) * u_{n-1} ]
    // With M=I: u_{n+1} = dt^2 * [ (2*I/dt^2 - K) * u_n - (I/dt^2) * u_{n-1} ]
    //          = 2*u_n - dt^2*K*u_n - u_{n-1}

    // So the recurrence is: u_{n+1} = 2*u_n - u_{n-1} - dt^2*K*u_n

    // Initial conditions give us u0 and v0.
    // We need u_{-1} to start the recurrence.
    // From central difference: v0 = (u1 - u_{-1})/(2*dt), but we don't know u1 yet.
    // Standard approach: u1 = u0 + dt*v0 + (dt^2/2)*a0, where a0 = -M^{-1}*K*u0 = -K*u0
    // So u1 = u0 + dt*v0 - (dt^2/2)*K*u0

    // But the problem states the exact recurrence used by TRD1C, so let's use that.
    // And it says "Derive the first step from the given initial conditions", meaning
    // we should compute u1 using v0 and the equation of motion.

    // From the equation of motion at t=0: M*a0 + C*v0 + K*u0 = 0 => a0 = -K*u0 (since M=I, C=0)
    // Then using Taylor: u1 = u0 + dt*v0 + (dt^2/2)*a0 = u0 + dt*v0 - (dt^2/2)*K*u0

    // However, looking more carefully at the Fortran code in FORM1:
    //   Line: U1(I) = U0(I)-DELTT*UDOT0(I)   [so u_{-1} = u0 - dt*v0]
    //   Then later in STEP, the recurrence u_{n+1} = 2*u_n - u_{n-1} - dt^2*K*u_n is used.
    // So with u_{-1} = u0 - dt*v0, then:
    //   u1 = 2*u0 - u_{-1} - dt^2*K*u0 = 2*u0 - (u0 - dt*v0) - dt^2*K*u0 = u0 + dt*v0 - dt^2*K*u0

    // This is different from the Taylor expansion (which has dt^2/2).
    // The central difference method uses a different approximation.

    // Actually, the standard central difference discretization of u'' is:
    //   u''_n ≈ (u_{n+1} - 2*u_n + u_{n-1}) / dt^2
    // So the equation M*u'' + K*u = 0 becomes:
    //   M*(u_{n+1} - 2*u_n + u_{n-1})/dt^2 + K*u_n = 0
    //   => u_{n+1} = 2*u_n - u_{n-1} - dt^2*M^{-1}*K*u_n

    // To start the iteration, we need u0 and u_{-1}.
    // Given u0 and v0, we can approximate u_{-1} using:
    //   v0 = (u1 - u_{-1})/(2*dt)  [central difference for first derivative]
    // But we don't know u1 yet.
    // Instead, use: v0 = (u0 - u_{-1})/dt  [backward difference] => u_{-1} = u0 - dt*v0

    // This is exactly what FORM1 does.

    // So with v0 = [0,0,0], we have u_{-1} = u0 = [1,0,0].

    // Therefore:
    //   u0 = [1,0,0]
    //   u_{-1} = [1,0,0]
    //   u1 = 2*u0 - u_{-1} - dt^2*K*u0 = u0 - dt^2*K*u0

    Vector u_minus1 = u0 - dt * v0;  // From FORM1: U1(I) = U0(I)-DELTT*UDOT0(I)

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