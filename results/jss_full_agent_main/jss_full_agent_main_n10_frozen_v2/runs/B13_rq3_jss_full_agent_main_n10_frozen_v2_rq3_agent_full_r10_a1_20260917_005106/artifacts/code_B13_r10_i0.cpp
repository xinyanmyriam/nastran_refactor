#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

// SROD1 - Rod element stress recovery (phase I)
// Translated from NASTRAN-95 Fortran 77.

int main() {
    const double L      = 2.0;
    const double A      = 0.01;
    const double E      = 2.1e11;
    const double nu     = 0.3;
    const double J      = 5e-6;
    const double C      = 0.005;

    Eigen::Matrix<double,6,1> ua, ub;
    ua.setZero();
    ub.setZero();
    ub(0) = 0.001;

    Eigen::Vector3d pa(0.0, 0.0, 0.0);
    Eigen::Vector3d pb(L,   0.0, 0.0);

    Eigen::Vector3d XN = pb - pa;
    double XL = XN.norm();
    XN /= XL;

    double EOVERL = E / XL;
    double G = E / (2.0 * (1.0 + nu));
    double GCOVRL = G * J / XL;

    Eigen::Vector3d SAT = XN * EOVERL;
    Eigen::Vector3d SBT = -XN * EOVERL;

    Eigen::Vector3d ua_t = ua.head<3>();
    Eigen::Vector3d ub_t = ub.head<3>();

    // Axial force includes area factor; sign convention gives positive for extension
    double axial_force = -A * (SAT.dot(ua_t) + SBT.dot(ub_t));
    double axial_stress = axial_force / A;

    double dtheta = ub(3) - ua(3);
    double torsional_stress = G * C * (dtheta / XL);

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{"
              << "\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":" << axial_force << ","
              << "\"torsional_stress\":" << torsional_stress
              << "}" << std::endl;

    return 0;
}