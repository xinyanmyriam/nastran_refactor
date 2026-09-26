#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

// Constant Strain Triangle (CST) membrane stiffness, 6x6, in global coords.
// xy: 3x2 matrix, one row per node (x,y). Plane-stress constitutive.
Matrix<double,6,6> cst_stiffness(const Matrix<double,3,2>& xy,
                                 double E, double nu, double t) {
    double x1 = xy(0,0), y1 = xy(0,1);
    double x2 = xy(1,0), y2 = xy(1,1);
    double x3 = xy(2,0), y3 = xy(2,1);

    // Signed area (twice)
    double A = 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));

    // Strain-displacement matrix B (3x6)
    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

    Matrix<double,3,6> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    B /= (2.0 * A);

    // Plane-stress constitutive matrix D (3x3)
    double cc = E / (1.0 - nu * nu);
    Matrix3d D;
    D << cc,      cc * nu, 0.0,
         cc * nu, cc,      0.0,
         0.0,     0.0,     cc * (1.0 - nu) / 2.0;

    return t * std::abs(A) * B.transpose() * D * B;
}

int main() {
    // Test case geometry: rectangle 2.0 x 1.5
    Matrix<double,4,2> nodes;
    nodes << 0.0, 0.0,
             2.0, 0.0,
             2.0, 1.5,
             0.0, 1.5;

    double E  = 200e9;   // Pa
    double nu = 0.3;
    double t  = 0.01;    // m

    // NASTRAN CQUAD4 crossed-diagonal decomposition (map M):
    //   T1: 1,2,4   T2: 2,3,1   T3: 3,4,2   T4: 4,1,3
    int tris[4][3] = {{0,1,3}, {1,2,0}, {2,3,1}, {3,0,2}};

    Matrix<double,8,8> K = Matrix<double,8,8>::Zero();

    for (int tr = 0; tr < 4; ++tr) {
        Matrix<double,3,2> xy;
        for (int i = 0; i < 3; ++i)
            xy.row(i) = nodes.row(tris[tr][i]);

        Matrix<double,6,6> kt = cst_stiffness(xy, E, nu, t);

        // Assemble triangle DOFs (2 per node) into the 8x8 quad matrix
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                int gi = 2 * tris[tr][i];
                int gj = 2 * tris[tr][j];
                K(gi,     gj)     += kt(2*i,   2*j);
                K(gi,     gj + 1) += kt(2*i,   2*j + 1);
                K(gi + 1, gj)     += kt(2*i+1, 2*j);
                K(gi + 1, gj + 1) += kt(2*i+1, 2*j + 1);
            }
        }
    }

    // Emit as JSON
    std::cout << std::setprecision(6) << std::scientific;
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 8; ++i) {
        std::cout << "[";
        for (int j = 0; j < 8; ++j) {
            std::cout << K(i, j);
            if (j < 7) std::cout << ",";
        }
        std::cout << "]";
        if (i < 7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}