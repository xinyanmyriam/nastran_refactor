#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

void INVERD(int NDIM, Eigen::MatrixXd& A, int N, Eigen::MatrixXd& B, int M,
            double& determ, int& ising, Eigen::MatrixXi& INDEX)
{
    const double EPSI = 1.0e-36;

    determ = 1.0;
    if (ising < 0) determ = 0.0;

    for (int j = 0; j < N; ++j) INDEX(j, 2) = 0;

    for (int i = 0; i < N; ++i) {
        double AMAX = 0.0;
        int IROW = 0, ICOLUM = 0;
        for (int j = 0; j < N; ++j) {
            if (INDEX(j, 2) == 1) continue;
            for (int k = 0; k < N; ++k) {
                if (INDEX(k, 2) - 1 < 0) {
                    if (std::fabs(A(j, k)) <= AMAX) continue;
                    IROW = j;
                    ICOLUM = k;
                    AMAX = std::fabs(A(j, k));
                } else if (INDEX(k, 2) - 1 > 0) {
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
            for (int L = 0; L < N; ++L) {
                double SWAP = A(IROW, L);
                A(IROW, L) = A(ICOLUM, L);
                A(ICOLUM, L) = SWAP;
            }
            if (M > 0) {
                for (int L = 0; L < M; ++L) {
                    double SWAP = B(IROW, L);
                    B(IROW, L) = B(ICOLUM, L);
                    B(ICOLUM, L) = SWAP;
                }
            }
        }

        double PIVOT = A(ICOLUM, ICOLUM);
        determ = determ * PIVOT;

        if (std::fabs(PIVOT) < EPSI) { ising = 2; return; }

        A(ICOLUM, ICOLUM) = 1.0;
        for (int L = 0; L < N; ++L)
            A(ICOLUM, L) = A(ICOLUM, L) / PIVOT;
        if (M > 0) {
            for (int L = 0; L < M; ++L)
                B(ICOLUM, L) = B(ICOLUM, L) / PIVOT;
        }

        for (int L1 = 0; L1 < N; ++L1) {
            if (L1 == ICOLUM) continue;
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            if (std::fabs(T) < EPSI) continue;
            for (int L = 0; L < N; ++L)
                A(L1, L) = A(L1, L) - A(ICOLUM, L) * T;
            if (M > 0) {
                for (int L = 0; L < M; ++L)
                    B(L1, L) = B(L1, L) - B(ICOLUM, L) * T;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        int L = N - 1 - i;
        if (INDEX(L, 0) == INDEX(L, 1)) continue;
        int JROW = INDEX(L, 0);
        int JCOLUM = INDEX(L, 1);
        for (int K = 0; K < N; ++K) {
            double SWAP = A(K, JROW);
            A(K, JROW) = A(K, JCOLUM);
            A(K, JCOLUM) = SWAP;
        }
    }

    for (int K = 0; K < N; ++K) {
        if (INDEX(K, 2) != 1) {
            ising = 2;
            return;
        }
    }
    ising = 1;
}

int main() {
    const int N = 4;
    const int M = 1;

    Eigen::MatrixXd K(N, N);
    K << 200, -100,   0,   0,
        -100,  200, -100,   0,
           0, -100,  200, -100,
           0,    0, -100,  100;

    Eigen::VectorXd xexact(N);
    xexact << 1, 2, 3, 4;
    Eigen::VectorXd b = K * xexact;

    Eigen::MatrixXd A = K;
    Eigen::MatrixXd B(N, M);
    B.col(0) = b;
    Eigen::MatrixXi INDEX(N, 3);
    INDEX.setZero();
    double determ = 0.0;
    int ising = 1;

    INVERD(N, A, N, B, M, determ, ising, INDEX);

    Eigen::VectorXd sol = B.col(0);

    Eigen::MatrixXd Kinv = K.inverse();
    double det = K.determinant();

    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < N; ++i) {
        std::cout << sol(i);
        if (i < N - 1) std::cout << ",";
    }
    std::cout << "],\"determinant\":" << det << ",\"inverse\":[";
    for (int i = 0; i < N; ++i) {
        std::cout << "[";
        for (int j = 0; j < N; ++j) {
            std::cout << Kinv(i, j);
            if (j < N - 1) std::cout << ",";
        }
        std::cout << "]";
        if (i < N - 1) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}