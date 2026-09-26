#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

struct InverdResult {
    Eigen::MatrixXd inverse;
    Eigen::MatrixXd solution;
    double determ;
    int ising;
};

InverdResult inverd(const Eigen::MatrixXd& A_in, const Eigen::MatrixXd& B_in, int N, int M) {
    const double EPSI = 1.0e-36;
    Eigen::MatrixXd A = A_in;
    Eigen::MatrixXd B = B_in;
    double determ = 1.0;
    int ising = 1;
    std::vector<int> idx1(N + 1, 0);
    std::vector<int> idx2(N + 1, 0);
    std::vector<int> idx3(N + 1, 0);
    for (int j = 1; j <= N; ++j) idx3[j] = 0;
    for (int i = 1; i <= N; ++i) {
        double amax = 0.0;
        int irow = 0, icolum = 0;
        for (int j = 1; j <= N; ++j) {
            if (idx3[j] == 1) continue;
            for (int k = 1; k <= N; ++k) {
                if (idx3[k] - 1 < 0) {
                    if (std::fabs(A(j - 1, k - 1)) <= amax) continue;
                    irow = j; icolum = k;
                    amax = std::fabs(A(j - 1, k - 1));
                } else if (idx3[k] - 1 > 0) {
                    ising = 2;
                    InverdResult r; r.inverse = A; r.solution = B; r.determ = determ; r.ising = ising;
                    return r;
                }
            }
        }
        idx3[icolum] = idx3[icolum] + 1;
        idx1[i] = irow; idx2[i] = icolum;
        if (irow != icolum) {
            determ = -determ;
            for (int l = 1; l <= N; ++l) {
                double swap = A(irow - 1, l - 1);
                A(irow - 1, l - 1) = A(icolum - 1, l - 1);
                A(icolum - 1, l - 1) = swap;
            }
            if (M > 0) {
                for (int l = 1; l <= M; ++l) {
                    double swap = B(irow - 1, l - 1);
                    B(irow - 1, l - 1) = B(icolum - 1, l - 1);
                    B(icolum - 1, l - 1) = swap;
                }
            }
        }
        double pivot = A(icolum - 1, icolum - 1);
        determ = determ * pivot;
        if (std::fabs(pivot) < EPSI) {
            ising = 2;
            InverdResult r; r.inverse = A; r.solution = B; r.determ = determ; r.ising = ising;
            return r;
        }
        A(icolum - 1, icolum - 1) = 1.0;
        for (int l = 1; l <= N; ++l) A(icolum - 1, l - 1) = A(icolum - 1, l - 1) / pivot;
        if (M > 0) for (int l = 1; l <= M; ++l) B(icolum - 1, l - 1) = B(icolum - 1, l - 1) / pivot;
        for (int l1 = 1; l1 <= N; ++l1) {
            if (l1 == icolum) continue;
            double t = A(l1 - 1, icolum - 1);
            A(l1 - 1, icolum - 1) = 0.0;
            if (std::fabs(t) < EPSI) continue;
            for (int l = 1; l <= N; ++l) A(l1 - 1, l - 1) = A(l1 - 1, l - 1) - A(icolum - 1, l - 1) * t;
            if (M > 0) for (int l = 1; l <= M; ++l) B(l1 - 1, l - 1) = B(l1 - 1, l - 1) - B(icolum - 1, l - 1) * t;
        }
    }
    for (int i = 1; i <= N; ++i) {
        int l = N + 1 - i;
        if (idx1[l] == idx2[l]) continue;
        int jrow = idx1[l]; int jcolum = idx2[l];
        for (int k = 1; k <= N; ++k) {
            double swap = A(k - 1, jrow - 1);
            A(k - 1, jrow - 1) = A(k - 1, jcolum - 1);
            A(k - 1, jcolum - 1) = swap;
        }
    }
    for (int k = 1; k <= N; ++k) {
        if (idx3[k] != 1) {
            ising = 2;
            InverdResult r; r.inverse = A; r.solution = B; r.determ = determ; r.ising = ising;
            return r;
        }
    }
    ising = 1;
    InverdResult r; r.inverse = A; r.solution = B; r.determ = determ; r.ising = ising;
    return r;
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
    Eigen::MatrixXd B(N, M);
    B.col(0) = b;
    InverdResult res = inverd(K, B, N, M);
    Eigen::VectorXd x = res.solution.col(0);
    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"INVERD\",\"solution\":[";
    for (int i = 0; i < N; ++i) { if (i) std::cout << ","; std::cout << x(i); }
    std::cout << "],\"determinant\":" << res.determ << ",\"inverse\":[";
    for (int i = 0; i < N; ++i) {
        if (i) std::cout << ",";
        std::cout << "[";
        for (int j = 0; j < N; ++j) { if (j) std::cout << ","; std::cout << res.inverse(i, j); }
        std::cout << "]";
    }
    std::cout << "]}" << std::endl;
    return 0;
}