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
static double inverd(MatrixXd& A, MatrixXd& B, int n, int m, int& ising, bool computeDeterm)
{
    const double EPSI = 1.0e-36;
    double determ = computeDeterm ? 1.0 : 0.0;

    // INDEX stored as three integer columns (1-based in Fortran -> 0-based here)
    std::vector<int> idx1(n, 0), idx2(n, 0), idx3(n, 0);

    for (int i = 0; i < n; ++i) {
        // Search for pivot
        double amax = 0.0;
        int irow = 0, icolum = 0;
        for (int j = 0; j < n; ++j) {
            if (idx3[j] == 1) continue;
            for (int k = 0; k < n; ++k) {
                if (idx3[k] - 1 < 0) {
                    // idx3[k] == 0
                    if (std::fabs(A(j, k)) <= amax) continue;
                    irow = j;
                    icolum = k;
                    amax = std::fabs(A(j, k));
                } else if (idx3[k] - 1 > 0) {
                    // singular
                    ising = 2;
                    return determ;
                }
            }
        }
        idx3[icolum] += 1;
        idx1[i] = irow;
        idx2[i] = icolum;

        // Interchange rows to put pivot on diagonal
        if (irow != icolum) {
            determ = -determ;
            for (int l = 0; l < n; ++l) {
                std::swap(A(irow, l), A(icolum, l));
            }
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    std::swap(B(irow, l), B(icolum, l));
                }
            }
        }

        double pivot = A(icolum, icolum);
        determ *= pivot;

        if (std::fabs(pivot) < EPSI) {
            ising = 2;
            return determ;
        }

        A(icolum, icolum) = 1.0;
        for (int l = 0; l < n; ++l) {
            A(icolum, l) /= pivot;
        }
        if (m > 0) {
            for (int l = 0; l < m; ++l) {
                B(icolum, l) /= pivot;
            }
        }

        // Reduce non-pivot rows
        for (int l1 = 0; l1 < n; ++l1) {
            if (l1 == icolum) continue;
            double t = A(l1, icolum);
            A(l1, icolum) = 0.0;
            if (std::fabs(t) < EPSI) continue;
            for (int l = 0; l < n; ++l) {
                A(l1, l) -= A(icolum, l) * t;
            }
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    B(l1, l) -= B(icolum, l) * t;
                }
            }
        }
    }

    // Interchange columns
    for (int i = 0; i < n; ++i) {
        int l = n - 1 - i;
        if (idx1[l] == idx2[l]) continue;
        int jrow = idx1[l];
        int jcolum = idx2[l];
        for (int k = 0; k < n; ++k) {
            std::swap(A(k, jrow), A(k, jcolum));
        }
    }

    // Check singularity
    for (int k = 0; k < n; ++k) {
        if (idx3[k] != 1) {
            ising = 2;
            return determ;
        }
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
    MatrixXd Bempty(n, 0);
    int ising2 = 0;
    double det2 = inverd(Ainv, Bempty, n, 0, ising2, true);

    std::cout << std::setprecision(17);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < n; ++i) {
        if (i) std::cout << ",";
        std::cout << x(i);
    }
    std::cout << "],\"determinant\":" << det << ",\"inverse\":[";
    for (int i = 0; i < n; ++i) {
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j) std::cout << ",";
            std::cout << Ainv(i, j);
        }
        std::cout << "]";
        if (i + 1 < n) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}