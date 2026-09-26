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

    // Vector along rod, length, normalize
    Eigen::Vector3d XN = pb - pa;
    double XL = XN.norm();
    XN /= XL;

    double EOVERL = E / XL;

    // Shear modulus
    double G = E / (2.0 * (1.0 + nu));

    // GCOVRL = G * C / XL
    double GCOVRL = G * C / XL;

    // ---- Axial force and stress ----
    // Relative axial displacement along rod axis (translational DOFs 0..2)
    Eigen::Vector3d ta(ua(0), ua(1), ua(2));
    Eigen::Vector3d tb(ub(0), ub(1), ub(2));
    double delta_axial = (tb - ta).dot(XN);   // = 0.001
    double axial_strain = delta_axial / XL;
    double axial_stress = E * axial_strain;
    double axial_force  = axial_stress * A;

    // ---- Torsional stress ----
    // Relative rotation about rod axis (rotational DOFs 3..5)
    Eigen::Vector3d ra(ua(3), ua(4), ua(5));
    Eigen::Vector3d rb(ub(3), ub(4), ub(5));
    double dtheta = (rb - ra).dot(XN);   // = 0
    double torsional_stress = G * (dtheta / XL) * C;

    // ---- Output JSON ----
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":"  << axial_force  << ","
              << "\"torsional_stress\":" << torsional_stress
              << "}" << std::endl;

    return 0;
}