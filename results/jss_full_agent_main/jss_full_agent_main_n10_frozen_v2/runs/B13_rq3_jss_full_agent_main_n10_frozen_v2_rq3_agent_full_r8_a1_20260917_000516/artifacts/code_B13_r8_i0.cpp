#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

// SROD1 - Rod element stress recovery (phase I)
// Translated from NASTRAN-95 Fortran 77.

int main() {
    // ---- Test case inputs ----
    const double L    = 2.0;
    const double A    = 0.01;
    const double E    = 2.1e11;
    const double nu   = 0.3;
    const double J    = 5e-6;
    const double C    = 0.005;   // J / max_radius

    // Node displacements (6 DOF each): pure axial
    // node_a = (0,0,0,0,0,0), node_b = (0.001,0,0,0,0,0)
    Eigen::Matrix<double,6,1> ua, ub;
    ua.setZero();
    ub.setZero();
    ub(0) = 0.001;

    // ---- Geometry: rod along x-axis ----
    // Node coordinates
    Eigen::Vector3d pa(0.0, 0.0, 0.0);
    Eigen::Vector3d pb(L,   0.0, 0.0);

    // Vector along rod, length, normalize
    Eigen::Vector3d XN = pb - pa;
    double XL = XN.norm();
    XN /= XL;

    double EOVERL = E / XL;

    // Shear modulus
    double G = E / (2.0 * (1.0 + nu));

    // GCOVRL = G * J / (XL * C)  -- from Fortran: G * ECPT(6) / XL
    // ECPT(6) is the polar moment / radius ratio (C). Here ECPT(6)=C.
    double GCOVRL = G * C / XL;

    // ---- Axial stiffness terms ----
    // SAT = XN * EOVERL, SBT = -XN * EOVERL
    Eigen::Vector3d SAT = XN * EOVERL;
    Eigen::Vector3d SBT = -XN * EOVERL;

    // ---- Axial force and stress ----
    // Axial strain from relative axial displacement along rod axis
    // delta = (ub - ua) . XN
    Eigen::Vector3d da = ua.head<3>();
    Eigen::Vector3d db = ub.head<3>();
    double delta = (db - da).dot(XN);

    double axial_strain = delta / XL;
    double axial_stress = E * axial_strain;
    double axial_force  = axial_stress * A;

    // ---- Torsional stress ----
    // Rotation about rod axis: relative twist
    // theta = (rot_b - rot_a) . XN
    Eigen::Vector3d ra = ua.tail<3>();
    Eigen::Vector3d rb = ub.tail<3>();
    double theta = (rb - ra).dot(XN);

    // Torsional shear stress = G * theta * C / XL
    double torsional_stress = G * theta * C / XL;

    // ---- Output JSON ----
    std::cout << std::setprecision(10);
    std::cout << std::scientific;
    std::cout << "{"
              << "\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":"  << axial_force  << ","
              << "\"torsional_stress\":" << torsional_stress
              << "}" << std::endl;

    return 0;
}