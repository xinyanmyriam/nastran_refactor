#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace Eigen;

// Translation of NASTRAN INVERD: in-place Gauss-Jordan inverse/solver.
// A: n x n matrix (will be replaced by its inverse)
// B: n x m right-hand sides (will be replaced by solutions)
// Returns determinant; sets ising (1 = ok, 2 = singular).
static double inverd(MatrixXd& A, MatrixXd& B, int n, int m, int& ising, bool computeDet)
{
    const double EPSI = 1.0e-36;
    double determ = computeDet ? 1.0 : 0.0;

    // INDEX stored as (n,3): col0 = irow, col1 = icolum, col2 = flag
    MatrixXi INDEX = MatrixXi::Zero(n, 3);

    for (int i = 0; i < n; ++i)
    {
        // Search for pivot
        double amax = 0.0;
        int irow = -1, icolum = -1;
        for (int j = 0; j < n; ++j)
        {
            if (INDEX(j, 2) == 1) continue;
            for (int k = 0; k < n; ++k)
            {
                if (INDEX(k, 2) - 1 > 0) { ising = 2; return determ; }
                if (INDEX(k, 2) - 1 == 0) continue; // ==1 handled above; here flag==0
                if (std::fabs(A(j, k)) <= amax) continue;
                irow = j;
                icolum = k;
                amax = std::fabs(A(j, k));
            }
        }
        if (irow < 0 || icolum < 0) { ising = 2; return determ; }

        INDEX(icolum, 2) += 1;
        INDEX(i, 0) = irow;
        INDEX(i, 1) = icolum;

        // Interchange rows to put pivot on diagonal
        if (irow != icolum)
        {
            determ = -determ;
            for (int l = 0; l < n; ++l)
                std::swap(A(irow, l), A(icolum, l));
            if (m > 0)
                for (int l = 0; l < m; ++l)
                    std::swap(B(irow, l), B(icolum, l));
        }

        double pivot = A(icolum, icolum);
        determ *= pivot;

        if (std::fabs(pivot) < EPSI) { ising = 2; return determ; }

        A(icolum, icolum) = 1.0;
        for (int l = 0; l < n; ++l)
            A(icolum, l) /= pivot;
        if (m > 0)
            for (int l = 0; l < m; ++l)
                B(icolum, l) /= pivot;

        // Reduce non-pivot rows
        for (int l1 = 0; l1 < n; ++l1)
        {
            if (l1 == icolum) continue;
            double t = A(l1, icolum);
            A(l1, icolum) = 0.0;
            if (std::fabs(t) < EPSI) continue;
            for (int l = 0; l < n; ++l)
                A(l1, l) -= A(icolum, l) * t;
            if (m > 0)
                for (int l = 0; l < m; ++l)
                    B(l1, l) -= B(icolum, l) * t;
        }
    }

    // Interchange columns
    for (int i = 0; i < n; ++i)
    {
        int l = n - 1 - i;
        if (INDEX(l, 0) == INDEX(l, 1)) continue;
        int jrow = INDEX(l, 0);
        int jcolum = INDEX(l, 1);
        for (int k = 0; k < n; ++k)
            std::swap(A(k, jrow), A(k, jcolum));
    }

    // Check singularity
    for (int k = 0; k < n; ++k)
    {
        if (INDEX(k, 2) != 1) { ising = 2; return determ; }
    }
    ising = 1;
    return determ;
}

int main()
{
    const int n = 4;
    MatrixXd K(n, n);
    K << 200, -100,    0,    0,
        -100,  200, -100,    0,
           0, -100,  200, -100,
           0,    0, -100,  100;

    VectorXd xexact(n);
    xexact << 1, 2, 3, 4;
    VectorXd b = K * xexact;

    // Solve
    MatrixXd A = K;
    MatrixXd B(n, 1);
    B.col(0) = b;
    int ising = 0;
    double det = inverd(A, B, n, 1, ising, true);
    VectorXd x = B.col(0);

    // Inverse
    MatrixXd Ainv = K;
    MatrixXd Bzero(n, 0);
    int ising2 = 0;
    inverd(Ainv, Bzero, n, 0, ising2, false);

    std::cout << std::setprecision(17);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < n; ++i)
    {
        if (i) std::cout << ",";
        std::cout << x(i);
    }
    std::cout << "],\"determinant\":" << det << ",\"inverse\":[";
    for (int i = 0; i < n; ++i)
    {
        std::cout << "[";
        for (int j = 0; j < n; ++j)
        {
            if (j) std::cout << ",";
            std::cout << Ainv(i, j);
        }
        std::cout << "]";
        if (i + 1 < n) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}