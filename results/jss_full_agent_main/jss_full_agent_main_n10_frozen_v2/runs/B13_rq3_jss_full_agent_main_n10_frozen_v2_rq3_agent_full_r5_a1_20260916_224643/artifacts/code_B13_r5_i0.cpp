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
    Eigen::Matrix<double,6,1> ua, ub;
    ua.setZero();
    ub.setZero();
    ub(0) = 0.001;

    // ---- Geometry: rod along x-axis ----
    Eigen::Vector3d pa(0.0, 0.0, 0.0);
    Eigen::Vector3d pb(L,   0.0, 0.0);

    // Vector along rod, normalized
    Eigen::Vector3d xn = pb - pa;
    double XL = xn.norm();
    xn /= XL;

    double EOVERL = E / XL;
    double G = E / (2.0 * (1.0 + nu));
    double GCOVRL = G * J / XL;

    // ---- Axial force and stress ----
    Eigen::Vector3d da = ua.head<3>();
    Eigen::Vector3d db = ub.head<3>();
    double delta_axial = (db - da).dot(xn);
    double axial_strain = delta_axial / XL;
    double axial_stress = E * axial_strain;
    double axial_force  = axial_stress * A;

    // ---- Torsional stress ----
    Eigen::Vector3d rot_a(ua(3), ua(4), ua(5));
    Eigen::Vector3d rot_b(ub(3), ub(4), ub(5));
    double twist = (rot_b - rot_a).dot(xn);
    double torsional_stress = G * C * twist / XL;

    // ---- Output JSON ----
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":"  << axial_force  << ","
              << "\"torsional_stress\":" << torsional_stress
              << "}" << std::endl;

    return 0;
}