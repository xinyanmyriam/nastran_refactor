#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

static void INVERD(int NDIM, Eigen::MatrixXd& A, int N,
                   Eigen::MatrixXd& B, int M,
                   double& DETERM, int& ISING)
{
    Eigen::MatrixXi INDEX = Eigen::MatrixXi::Zero(N, 3);
    double EPSI = 1.0e-36;
    DETERM = 1.0;
    if (ISING < 0) DETERM = 0.0;
    for (int j = 0; j < N; ++j) INDEX(j, 2) = 0;
    for (int I = 0; I < N; ++I)
    {
        double AMAX = 0.0;
        int IROW = 0, ICOLUM = 0;
        for (int J = 0; J < N; ++J)
        {
            if (INDEX(J, 2) == 1) continue;
            for (int K = 0; K < N; ++K)
            {
                if (INDEX(K, 2) - 1 < 0)
                {
                    if (std::fabs(A(J, K)) <= AMAX) continue;
                    IROW = J; ICOLUM = K; AMAX = std::fabs(A(J, K));
                }
                else if (INDEX(K, 2) - 1 > 0) { ISING = 2; return; }
            }
        }
        INDEX(ICOLUM, 2) = INDEX(ICOLUM, 2) + 1;
        INDEX(I, 0) = IROW;
        INDEX(I, 1) = ICOLUM;
        if (IROW != ICOLUM)
        {
            DETERM = -DETERM;
            for (int L = 0; L < N; ++L)
            {
                double SWAP = A(IROW, L);
                A(IROW, L) = A(ICOLUM, L);
                A(ICOLUM, L) = SWAP;
            }
            if (M > 0)
                for (int L = 0; L < M; ++L)
                {
                    double SWAP = B(IROW, L);
                    B(IROW, L) = B(ICOLUM, L);
                    B(ICOLUM, L) = SWAP;
                }
        }
        double PIVOT = A(ICOLUM, ICOLUM);
        DETERM = DETERM * PIVOT;
        if (std::fabs(PIVOT) < EPSI) { ISING = 2; return; }
        A(ICOLUM, ICOLUM) = 1.0;
        for (int L = 0; L < N; ++L) A(ICOLUM, L) = A(ICOLUM, L) / PIVOT;
        if (M > 0)
            for (int L = 0; L < M; ++L) B(ICOLUM, L) = B(ICOLUM, L) / PIVOT;
        for (int L1 = 0; L1 < N; ++L1)
        {
            if (L1 == ICOLUM) continue;
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            if (std::fabs(T) < EPSI) continue;
            for (int L = 0; L < N; ++L) A(L1, L) = A(L1, L) - A(ICOLUM, L) * T;
            if (M > 0)
                for (int L = 0; L < M; ++L) B(L1, L) = B(L1, L) - B(ICOLUM, L) * T;
        }
    }
    for (int I = 0; I < N; ++I)
    {
        int L = N - 1 - I;
        if (INDEX(L, 0) == INDEX(L, 1)) continue;
        int JROW = INDEX(L, 0);
        int JCOLUM = INDEX(L, 1);
        for (int K = 0; K < N; ++K)
        {
            double SWAP = A(K, JROW);
            A(K, JROW) = A(K, JCOLUM);
            A(K, JCOLUM) = SWAP;
        }
    }
    for (int K = 0; K < N; ++K)
        if (INDEX(K, 2) != 1) { ISING = 2; return; }
    ISING = 1;
}

int main()
{
    Eigen::MatrixXd K(4, 4);
    K << 200, -100,   0,   0,
        -100,  200, -100,   0,
           0, -100,  200, -100,
           0,    0, -100,  100;
    Eigen::VectorXd xtrue(4);
    xtrue << 1, 2, 3, 4;
    Eigen::VectorXd b = K * xtrue;
    Eigen::MatrixXd A = K;
    Eigen::MatrixXd B(4, 1);
    B.col(0) = b;
    double DETERM = 0.0;
    int ISING = 1;
    INVERD(4, A, 4, B, 1, DETERM, ISING);
    Eigen::VectorXd x = B.col(0);
    Eigen::MatrixXd Kinv = A;
    std::cout << std::setprecision(17);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < 4; ++i) { std::cout << x(i); if (i < 3) std::cout << ","; }
    std::cout << "],\"determinant\":" << DETERM << ",\"inverse\":[";
    for (int i = 0; i < 4; ++i)
    {
        std::cout << "[";
        for (int j = 0; j < 4; ++j) { std::cout << Kinv(i, j); if (j < 3) std::cout << ","; }
        std::cout << "]";
        if (i < 3) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    return 0;
}