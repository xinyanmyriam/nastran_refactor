#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <iomanip>

int main() {
    // Test case parameters
    const double L = 2.0;
    const double A = 0.01;
    const double E = 2.1e11;
    const double nu = 0.3;
    const double J = 5e-6;
    const double C = 0.005; // J / max_radius

    // Displacements: node_a = (0,0,0,0,0,0), node_b = (0.001,0,0,0,0,0)
    // Pure axial displacement along x
    double delta_x = 0.001;

    // Axial strain = delta_x / L
    double strain = delta_x / L;

    // Axial stress = E * strain
    double axial_stress = E * strain;

    // Axial force = stress * area
    double axial_force = axial_stress * A;

    // Torsional stress: no rotation applied (pure axial), so zero
    double torsional_stress = 0.0;

    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"test\":\"SROD1\","
              << "\"axial_stress\":" << axial_stress << ","
              << "\"axial_force\":" << axial_force << ","
              << "\"torsional_stress\":" << torsional_stress << "}"
              << std::endl;

    return 0;
}