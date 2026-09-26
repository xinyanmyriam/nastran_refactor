#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

void INVERD(int NDIM, Eigen::MatrixXd& A, int N, Eigen::MatrixXd& B, int M,
            double& DETERM, int& ISING, Eigen::MatrixXi& INDEX)
{
    const double EPSI = 1.0e-36;
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
                int ik3 = INDEX(K, 2);
                if (ik3 - 1 < 0)
                {
                    if (std::fabs(A(J, K)) <= AMAX) continue;
                    IROW = J;
                    ICOLUM = K;
                    AMAX = std::fabs(A(J, K));
                }
                else if (ik3 - 1 > 0)
                {
                    ISING = 2;
                    return;
                }
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
            {
                for (int L = 0; L < M; ++L)
                {
                    double SWAP = B(IROW, L);
                    B(IROW, L) = B(ICOLUM, L);
                    B(ICOLUM, L) = SWAP;
                }
            }
        }

        double PIVOT = A(ICOLUM, ICOLUM);
        DETERM = DETERM * PIVOT;

        if (std::fabs(PIVOT) < EPSI)
        {
            ISING = 2;
            return;
        }
        A(ICOLUM, ICOLUM) = 1.0;
        for (int L = 0; L < N; ++L)
            A(ICOLUM, L) = A(ICOLUM, L) / PIVOT;
        if (M > 0)
        {
            for (int L = 0; L < M; ++L)
                B(ICOLUM, L) = B(ICOLUM, L) / PIVOT;
        }

        for (int L1 = 0; L1 < N; ++L1)
        {
            if (L1 == ICOLUM) continue;
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            if (std::fabs(T) < EPSI) continue;
            for (int L = 0; L < N; ++L)
                A(L1, L) = A(L1, L) - A(ICOLUM, L) * T;
            if (M > 0)
            {
                for (int L = 0; L < M; ++L)
                    B(L1, L) = B(L1, L) - B(ICOLUM, L) * T;
            }
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
    {
        if (INDEX(K, 2) != 1)
        {
            ISING = 2;
            return;
        }
    }
    ISING = 1;
}

int main()
{
    const int N = 4;
    Eigen::MatrixXd K(N, N);
    K << 200, -100,   0,   0,
        -100,  200, -100,   0,
           0, -100,  200, -100,
           0,    0, -100,  100;

    Eigen::VectorXd xexact(N);
    xexact << 1, 2, 3, 4;
    Eigen::VectorXd b = K * xexact;

    Eigen::MatrixXd A = K;
    Eigen::MatrixXd B(N, 1);
    B.col(0) = b;
    double DETERM = 0.0;
    int ISING = 1;
    Eigen::MatrixXi INDEX(N, 3);
    INDEX.setZero();

    INVERD(N, A, N, B, 1, DETERM, ISING, INDEX);

    Eigen::VectorXd solution = B.col(0);
    Eigen::MatrixXd inverse = A;

    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < N; ++i)
    {
        if (i) std::cout << ",";
        std::cout << solution(i);
    }
    std::cout << "],\"determinant\":" << DETERM << ",\"inverse\":[";
    for (int i = 0; i < N; ++i)
    {
        std::cout << "[";
        for (int j = 0; j < N; ++j)
        {
            if (j) std::cout << ",";
            std::cout << inverse(i, j);
        }
        std::cout << "]";
        if (i < N - 1) std::cout << ",";
    }
    std::cout << "]}" << std::endl;

    return 0;
}