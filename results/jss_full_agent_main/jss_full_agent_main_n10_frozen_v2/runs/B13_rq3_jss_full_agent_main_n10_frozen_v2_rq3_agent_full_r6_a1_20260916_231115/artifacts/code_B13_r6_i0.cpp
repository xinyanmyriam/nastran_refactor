#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

// SROD1 - Rod element stress recovery (phase I)
// Translated from NASTRAN-95 Fortran 77.

int main() {
    // ---- Test case inputs ----
    const double L      = 2.0;
    const double A      = 0.01;
    const double E      = 2.1e11;
    const double nu     = 0.3;
    const double J      = 5e-6;
    const double C      = 0.005;   // J / max_radius

    // Node displacements (6 DOF each): pure axial
    // node_a = (0,0,0,0,0,0), node_b = (0.001,0,0,0,0,0)
    Eigen::Matrix<double,6,1> ua, ub;
    ua.setZero();
    ub.setZero();
    ub(0) = 0.001;

    // ---- Geometry: rod along x-axis ----
    // Node A at origin, Node B at (L,0,0)
    Eigen::Vector3d pa(0.0, 0.0, 0.0);
    Eigen::Vector3d pb(L,   0.0, 0.0);

    // Vector along rod, normalized
    Eigen::Vector3d xn = pb - pa;
    double XL = xn.norm();
    xn /= XL;

    // EOVERL = E / XL
    double EOVERL = E / XL;

    // G = E / (2*(1+nu))
    double G = E / (2.0 * (1.0 + nu));

    // GCOVRL = G * J / XL  (ECPT(6) is the polar moment J)
    double GCOVRL = G * J / XL;

    // ---- Axial stiffness terms ----
    // SAT = xn * EOVERL  (force per unit displacement at A)
    Eigen::Vector3d SAT = xn * EOVERL;
    // SAR = xn * GCOVRL  (torsional terms)
    Eigen::Vector3d SAR = xn * GCOVRL;

    // SBT = -xn * EOVERL
    Eigen::Vector3d SBT = -xn * EOVERL;
    // SBR = -xn * GCOVRL
    Eigen::Vector3d SBR = -xn * GCOVRL;

    // ---- Stress recovery ----
    // Axial force: F = SAT . ua + SBT . ub  (using axial components)
    // In the rod formulation, axial force = (E*A/L) * (u_b - u_a) along axis
    double axial_force = (E * A / XL) * (ub(0) - ua(0));

    // Axial stress = force / area
    double axial_stress = axial_force / A;

    // Torsional stress: tau = G * C * (theta_b - theta_a) / L
    // theta = rotation about rod axis (DOF 4 for x-axis rod)
    double theta_a = ua(3);
    double theta_b = ub(3);
    double torsional_stress = G * C * (theta_b - theta_a) / XL;

    // ---- Output JSON ----
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":" << axial_force << ","
              << "\"torsional_stress\":" << torsional_stress << "}"
              << std::endl;

    return 0;
}