#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

void INVERD(int ndim, Eigen::MatrixXd& A, int n, Eigen::MatrixXd& B, int m,
            double& determ, int& ising, Eigen::MatrixXi& INDEX)
{
    const double EPSI_DEFAULT = 1.0e-36;
    double EPSI = EPSI_DEFAULT;

    determ = 1.0;
    if (ising < 0) determ = 0.0;

    for (int j = 0; j < n; ++j) INDEX(j, 2) = 0;

    for (int i = 0; i < n; ++i) {
        double AMAX = 0.0;
        int IROW = 0, ICOLUM = 0;
        for (int j = 0; j < n; ++j) {
            if (INDEX(j, 2) == 1) continue;
            for (int k = 0; k < n; ++k) {
                int ik3 = INDEX(k, 2);
                if (ik3 - 1 < 0) {
                    if (std::fabs(A(j, k)) <= AMAX) continue;
                    IROW = j;
                    ICOLUM = k;
                    AMAX = std::fabs(A(j, k));
                } else if (ik3 - 1 > 0) {
                    ising = 2;
                    return;
                }
            }
        }
        INDEX(ICOLUM, 2) = INDEX(ICOLUM, 2) + 1;
        INDEX(i, 0) = IROW;
        INDEX(i, 1) = ICOLUM;

        if (IROW != ICOLUM) {
            determ = -determ;
            for (int L = 0; L < n; ++L) {
                double SWAP = A(IROW, L);
                A(IROW, L) = A(ICOLUM, L);
                A(ICOLUM, L) = SWAP;
            }
            if (m > 0) {
                for (int L = 0; L < m; ++L) {
                    double SWAP = B(IROW, L);
                    B(IROW, L) = B(ICOLUM, L);
                    B(ICOLUM, L) = SWAP;
                }
            }
        }

        double PIVOT = A(ICOLUM, ICOLUM);
        determ = determ * PIVOT;

        if (std::fabs(PIVOT) < EPSI) {
            ising = 2;
            return;
        }
        A(ICOLUM, ICOLUM) = 1.0;
        for (int L = 0; L < n; ++L)
            A(ICOLUM, L) = A(ICOLUM, L) / PIVOT;
        if (m > 0) {
            for (int L = 0; L < m; ++L)
                B(ICOLUM, L) = B(ICOLUM, L) / PIVOT;
        }

        for (int L1 = 0; L1 < n; ++L1) {
            if (L1 == ICOLUM) continue;
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            if (std::fabs(T) < EPSI) continue;
            for (int L = 0; L < n; ++L)
                A(L1, L) = A(L1, L) - A(ICOLUM, L) * T;
            if (m > 0) {
                for (int L = 0; L < m; ++L)
                    B(L1, L) = B(L1, L) - B(ICOLUM, L) * T;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        int L = n - 1 - i;
        if (INDEX(L, 0) == INDEX(L, 1)) continue;
        int JROW = INDEX(L, 0);
        int JCOLUM = INDEX(L, 1);
        for (int K = 0; K < n; ++K) {
            double SWAP = A(K, JROW);
            A(K, JROW) = A(K, JCOLUM);
            A(K, JCOLUM) = SWAP;
        }
    }

    for (int K = 0; K < n; ++K) {
        if (INDEX(K, 2) != 1) {
            ising = 2;
            return;
        }
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

    Eigen::VectorXd sol = B.col(0);
    Eigen::MatrixXd Kinv = A;

    std::cout << std::setprecision(17);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < n; ++i) {
        if (i) std::cout << ",";
        std::cout << sol(i);
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