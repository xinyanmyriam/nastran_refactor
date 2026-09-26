#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace Eigen;

// Constant Strain Triangle (CST) membrane stiffness, 6x6 (2 DOF per node)
// nodes: 3x2 matrix of (x,y), E, nu, t
Matrix<double,6,6> cst_stiffness(const Matrix<double,3,2>& xy,
                                  double E, double nu, double t) {
    // Plane stress D matrix
    Matrix3d D;
    double c = E/(1.0 - nu*nu);
    D << c,      c*nu,   0.0,
         c*nu,   c,      0.0,
         0.0,    0.0,    c*(1.0-nu)/2.0;

    // B matrix (3x6)
    double x1=xy(0,0), y1=xy(0,1);
    double x2=xy(1,0), y2=xy(1,1);
    double x3=xy(2,0), y3=xy(2,1);

    double b1 = y2 - y3;
    double b2 = y3 - y1;
    double b3 = y1 - y2;
    double c1 = x3 - x2;
    double c2 = x1 - x3;
    double c3 = x2 - x1;

    // 2*Area
    double detJ = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1);
    double A = 0.5*detJ;

    Matrix<double,3,6> B;
    B << b1, 0.0, b2, 0.0, b3, 0.0,
         0.0, c1, 0.0, c2, 0.0, c3,
         c1, b1, c2, b2, c3, b3;
    B /= detJ;

    Matrix<double,6,6> K = t * A * (B.transpose() * D * B);
    return K;
}

int main() {
    // Test case
    double E = 200e9;
    double nu = 0.3;
    double t = 0.01;

    // Nodes (0,0),(2,0),(2,1.5),(0,1.5)
    Matrix<double,4,2> nodes;
    nodes << 0.0, 0.0,
             2.0, 0.0,
             2.0, 1.5,
             0.0, 1.5;

    // NASTRAN KQDMEM: pivot at node 1 (NPVT = first grid). Decompose into
    // 4 triangles using map M = {1,2,4, 2,3,1, 3,4,2, 4,1,3}
    // (1-based). Skip the triangle opposite the pivot (JNOT).
    // For pivot=1, JNOT = 3 (the triangle not containing node 1).
    int M[12] = {1,2,4, 2,3,1, 3,4,2, 4,1,3};

    Matrix<double,8,8> K = Matrix<double,8,8>::Zero();

    // Pivot = node 1 -> NPIVOT=1 -> JNOT = NPIVOT+2 = 3
    int NPIVOT = 1;
    int JNOT = NPIVOT + 2; // = 3

    for (int J = 1; J <= 4; ++J) {
        if (J == JNOT) continue;

        // Triangle J uses nodes M[3J-3 .. 3J-1] (1-based)
        int n0 = M[3*J-3] - 1;
        int n1 = M[3*J-2] - 1;
        int n2 = M[3*J-1] - 1;

        Matrix<double,3,2> xy;
        xy << nodes(n0,0), nodes(n0,1),
              nodes(n1,0), nodes(n1,1),
              nodes(n2,0), nodes(n2,1);

        Matrix<double,6,6> kt = cst_stiffness(xy, E, nu, t);

        // Assemble into global 8x8 (2 DOF per node)
        int gidx[3] = {n0, n1, n2};
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        K(2*gidx[a]+i, 2*gidx[b]+j) += kt(2*a+i, 2*b+j);
                    }
                }
            }
        }
    }

    // Print as JSON
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "{\"stiffness_matrix\":[";
    for (int i = 0; i < 8; ++i) {
        std::cout << "[";
        for (int j = 0; j < 8; ++j) {
            std::cout << K(i,j);
            if (j < 7) std::cout << ",";
        }
        std::cout << "]";
        if (i < 7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}