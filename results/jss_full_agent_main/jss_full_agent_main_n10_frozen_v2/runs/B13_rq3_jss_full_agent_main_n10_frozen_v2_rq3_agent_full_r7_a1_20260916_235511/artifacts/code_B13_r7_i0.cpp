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
    // SAT = xn * EOVERL  (force per unit axial displacement at A)
    Eigen::Vector3d SAT = xn * EOVERL;
    // SAR = xn * GCOVRL  (torsional terms)
    Eigen::Vector3d SAR = xn * GCOVRL;

    // ---- Compute axial force from displacement ----
    // Axial displacement along rod axis: delta = xn . (ub - ua)
    Eigen::Vector3d da = ua.head<3>();
    Eigen::Vector3d db = ub.head<3>();
    double delta = xn.dot(db - da);

    // Axial force = (E*A/L) * delta
    double axial_force = (E * A / XL) * delta;

    // Axial stress = force / area
    double axial_stress = axial_force / A;

    // ---- Torsional stress ----
    // Rotation about rod axis: relative twist
    Eigen::Vector3d ra = ua.tail<3>();
    Eigen::Vector3d rb = ub.tail<3>();
    double twist = xn.dot(rb - ra);   // rotation about rod axis

    // Torsional shear stress = G * twist * C / L
    double torsional_stress = G * twist * C / XL;

    // ---- Output JSON ----
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{"
              << "\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":" << axial_force << ","
              << "\"torsional_stress\":" << torsional_stress
              << "}" << std::endl;

    return 0;
}