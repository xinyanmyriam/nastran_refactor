#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>

// SROD1 - Rod element stress recovery (phase I)
// Translates the NASTRAN-95 Fortran routine.

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
    // Node coordinates
    Eigen::Vector3d pa(0.0, 0.0, 0.0);
    Eigen::Vector3d pb(L,   0.0, 0.0);

    // Vector along rod, length, normalize
    Eigen::Vector3d XN = pb - pa;
    double XL = XN.norm();
    XN /= XL;

    double EOVERL = E / XL;
    double G = E / (2.0 * (1.0 + nu));
    double GCOVRL = G * J / XL;   // ECPT(6) is the polar moment J

    // SAT/SAR (node A), SBT/SBR (node B)
    Eigen::Vector3d SAT = XN * EOVERL;
    Eigen::Vector3d SAR = XN * GCOVRL;
    Eigen::Vector3d SBT = -XN * EOVERL;
    Eigen::Vector3d SBR = -XN * GCOVRL;

    // ---- Stress recovery ----
    // Axial force: F = (E*A/L) * (u_b - u_a) projected on rod axis
    Eigen::Vector3d du = ub.head<3>() - ua.head<3>();
    double axial_elong = XN.dot(du);
    double axial_force  = (E * A / XL) * axial_elong;
    double axial_stress = axial_force / A;

    // Torsional stress: from rotational DOF difference about rod axis
    Eigen::Vector3d rot = ub.tail<3>() - ua.tail<3>();
    double twist = XN.dot(rot);          // rotation about rod axis
    double torsional_stress = G * C * twist / XL;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":"  << axial_force  << ","
              << "\"torsional_stress\":" << torsional_stress
              << "}" << std::endl;

    return 0;
}