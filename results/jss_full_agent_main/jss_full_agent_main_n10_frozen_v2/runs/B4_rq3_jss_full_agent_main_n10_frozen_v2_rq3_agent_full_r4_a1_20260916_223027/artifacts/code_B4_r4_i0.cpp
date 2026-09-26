#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;

// Triangular membrane (CST) stiffness, 2x2 blocks per node pair.
static void ktrmem(const Vector3d& p0, const Vector3d& p1, const Vector3d& p2,
                   double E, double nu, double t,
                   Matrix2d K[3][3])
{
    Vector3d e1 = p1 - p0;
    double L1 = e1.norm();
    e1 /= L1;
    Vector3d n = (p1 - p0).cross(p2 - p0);
    n.normalize();
    Vector3d e2 = n.cross(e1);

    auto to2d = [&](const Vector3d& p) -> Vector2d {
        Vector3d d = p - p0;
        return Vector2d(d.dot(e1), d.dot(e2));
    };

    Vector2d a = to2d(p0);
    Vector2d b = to2d(p1);
    Vector2d c = to2d(p2);

    double area = 0.5 * ((b.x()-a.x())*(c.y()-a.y()) - (c.x()-a.x())*(b.y()-a.y()));
    double A = std::abs(area);

    double x1=a.x(), y1=a.y(), x2=b.x(), y2=b.y(), x3=c.x(), y3=c.y();
    double b1 = y2 - y3, b2 = y3 - y1, b3 = y1 - y2;
    double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;
    double twoA = 2.0 * area;

    Matrix<double,3,6> B;
    B << b1, 0, b2, 0, b3, 0,
         0, c1, 0, c2, 0, c3,
         c1, b1, c2, b2, c3, b3;
    B /= twoA;

    double fac = E / (1.0 - nu*nu);
    Matrix3d D;
    D << fac,      fac*nu, 0,
         fac*nu,   fac,    0,
         0,        0,      fac*(1.0-nu)/2.0;

    Matrix<double,6,6> Ke = t * A * (B.transpose() * D * B);

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            K[i][j] = Ke.block<2,2>(2*i, 2*j);
}

int main()
{
    Vector3d nodes[4] = {
        Vector3d(0.0, 0.0, 0.0),
        Vector3d(2.0, 0.0, 0.0),
        Vector3d(2.0, 1.5, 0.0),
        Vector3d(0.0, 1.5, 0.0)
    };
    double E = 200e9;
    double nu = 0.3;
    double t = 0.01;

    // NASTRAN KQDMEM mapping table (1-based node indices per triangle)
    int M[12] = {1,2,4, 2,3,1, 3,4,2, 4,1,3};

    MatrixXd K8 = MatrixXd::Zero(8,8);

    for (int j = 0; j < 4; ++j) {
        int mp = 3*j;
        int tri[3] = {M[mp+0]-1, M[mp+1]-1, M[mp+2]-1};

        Matrix2d K[3][3];
        // NASTRAN convention: each sub-triangle uses thickness t/2
        ktrmem(nodes[tri[0]], nodes[tri[1]], nodes[tri[2]], E, nu, t/2.0, K);

        for (int i = 0; i < 3; ++i)
            for (int k = 0; k < 3; ++k)
                K8.block<2,2>(2*tri[i], 2*tri[k]) += K[i][k];
    }

    std::cout << std::setprecision(10);
    std::cout << "{\"stiffness_matrix\":[";
    for (int r = 0; r < 8; ++r) {
        std::cout << "[";
        for (int c = 0; c < 8; ++c) {
            std::cout << std::scientific << K8(r,c);
            if (c < 7) std::cout << ",";
        }
        std::cout << "]";
        if (r < 7) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}