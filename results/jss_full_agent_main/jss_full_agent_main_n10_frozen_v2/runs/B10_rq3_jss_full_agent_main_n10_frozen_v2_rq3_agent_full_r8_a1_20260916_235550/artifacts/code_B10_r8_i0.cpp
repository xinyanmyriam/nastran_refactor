#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

void INVERD(int ndim, Eigen::MatrixXd& A, int n,
            Eigen::MatrixXd& B, int m,
            double& determ, int& ising,
            Eigen::MatrixXi& INDEX)
{
    const double EPSI_DEFAULT = 1.0e-36;
    double epsi = EPSI_DEFAULT;

    determ = 1.0;
    if (ising < 0) determ = 0.0;

    for (int j = 0; j < n; ++j) INDEX(j, 2) = 0;

    for (int i = 0; i < n; ++i) {
        double amax = 0.0;
        int irow = 0, icolum = 0;
        for (int j = 0; j < n; ++j) {
            if (INDEX(j, 2) == 1) continue;
            for (int k = 0; k < n; ++k) {
                int ik3 = INDEX(k, 2);
                if (ik3 - 1 < 0) {
                    if (std::fabs(A(j, k)) <= amax) continue;
                    irow = j;
                    icolum = k;
                    amax = std::fabs(A(j, k));
                } else if (ik3 - 1 == 0) {
                    continue;
                } else {
                    ising = 2;
                    return;
                }
            }
        }
        INDEX(icolum, 2) = INDEX(icolum, 2) + 1;
        INDEX(i, 0) = irow;
        INDEX(i, 1) = icolum;

        if (irow != icolum) {
            determ = -determ;
            for (int l = 0; l < n; ++l) {
                double swap = A(irow, l);
                A(irow, l) = A(icolum, l);
                A(icolum, l) = swap;
            }
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    double swap = B(irow, l);
                    B(irow, l) = B(icolum, l);
                    B(icolum, l) = swap;
                }
            }
        }

        double pivot = A(icolum, icolum);
        determ = determ * pivot;

        if (std::fabs(pivot) < epsi) {
            ising = 2;
            return;
        }
        A(icolum, icolum) = 1.0;
        for (int l = 0; l < n; ++l)
            A(icolum, l) = A(icolum, l) / pivot;
        if (m > 0) {
            for (int l = 0; l < m; ++l)
                B(icolum, l) = B(icolum, l) / pivot;
        }

        for (int l1 = 0; l1 < n; ++l1) {
            if (l1 == icolum) continue;
            double t = A(l1, icolum);
            A(l1, icolum) = 0.0;
            if (std::fabs(t) < epsi) continue;
            for (int l = 0; l < n; ++l)
                A(l1, l) = A(l1, l) - A(icolum, l) * t;
            if (m > 0) {
                for (int l = 0; l < m; ++l)
                    B(l1, l) = B(l1, l) - B(icolum, l) * t;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        int l = n - 1 - i;
        if (INDEX(l, 0) == INDEX(l, 1)) continue;
        int jrow = INDEX(l, 0);
        int jcolum = INDEX(l, 1);
        for (int k = 0; k < n; ++k) {
            double swap = A(k, jrow);
            A(k, jrow) = A(k, jcolum);
            A(k, jcolum) = swap;
        }
    }

    for (int k = 0; k < n; ++k) {
        if (INDEX(k, 2) == 1) continue;
        ising = 2;
        return;
    }
    ising = 1;
}

int main() {
    const int n = 4;
    Eigen::MatrixXd K(n, n);
    K << 200, -100,   0,   0,
        -100,  200, -100,   0,
           0, -100,  200, -100,
           0,    0, -100,  100;

    Eigen::VectorXd xexact(n);
    xexact << 1, 2, 3, 4;
    Eigen::VectorXd b = K * xexact;

    Eigen::MatrixXd A = K;
    Eigen::MatrixXd B(n, 1);
    B.col(0) = b;
    double determ = 0.0;
    int ising = 1;
    Eigen::MatrixXi INDEX(n, 3);
    INDEX.setZero();

    INVERD(n, A, n, B, 1, determ, ising, INDEX);

    Eigen::VectorXd x = B.col(0);
    Eigen::MatrixXd Kinv = A;

    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < n; ++i) {
        if (i) std::cout << ",";
        std::cout << x(i);
    }
    std::cout << "],\"determinant\":" << determ << ",\"inverse\":[";
    for (int i = 0; i < n; ++i) {
        std::cout << "[";
        for (int j = 0; j < n; ++j) {
            if (j) std::cout << ",";
            std::cout << Kinv(i, j);
        }
        std::cout << "]";
        if (i + 1 < n) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}