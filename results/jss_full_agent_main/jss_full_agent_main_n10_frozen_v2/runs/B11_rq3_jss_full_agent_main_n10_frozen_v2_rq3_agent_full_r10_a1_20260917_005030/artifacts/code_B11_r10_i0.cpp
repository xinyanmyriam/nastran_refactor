#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

// QRITER: Ortega-Kaiser QR iteration for eigenvalues of a tridiagonal matrix.
// VAL = diagonal terms (in/out: eigenvalues on return)
// O   = squares of off-diagonal terms
// N   = order of the problem
static const double EPSI = 1.0e-10;
static const double ZERO = 0.0;
static const double ONE  = 1.0;

void qriter(std::vector<double>& VAL, std::vector<double>& O, int N) {
    int MAX = 100 * N;
    int NEVER = 0;
    int M2 = N;
    int M2M1, M1 = 1, MM = 1;

    while (true) {
        // Search for a decoupled submatrix
        M2M1 = M2 - 1;
        bool allZero = true;
        for (int K = 1; K <= M2M1; ++K) {
            M1 = M2 - K;
            if (O[M1] != ZERO) { allZero = false; break; }
        }
        if (allZero) break;  // all off-diagonal zero -> done

        // Decoupled submatrix
        M2M1 = M1;
        M2   = M1 + 1;
        if (M2M1 != 1) {
            bool found = false;
            for (int K = 2; K <= M2M1; ++K) {
                M1 = M2 - K;
                if (O[M1] == ZERO) { found = true; break; }
            }
            if (found) M1 = M1 + 1;
        }
        MM = M1;

        // QR iteration for the decoupled submatrix
        bool converged = false;
        for (int ITER = 1; ITER <= MAX; ++ITER) {
            if (std::abs(VAL[M2]) + O[M2M1] == std::abs(VAL[M2])) {
                converged = true;
                break;
            }
            // Check if all diagonal terms are equal
            bool allEq = true;
            for (int K = M1; K <= M2M1; ++K) {
                if (VAL[K] != VAL[K+1]) { allEq = false; break; }
            }
            double SHIFT;
            if (allEq) {
                SHIFT = ZERO;
            } else {
                // Find the smallest diagonal term = SHIFT
                SHIFT = VAL[M2];
                for (int I = M1; I <= M2M1; ++I) {
                    if (std::abs(VAL[I]) < std::abs(SHIFT)) SHIFT = VAL[I];
                }
                // Reduce all terms by SHIFT
                for (int I = M1; I <= M2; ++I) {
                    VAL[I] = VAL[I] - SHIFT;
                }
            }

            // QR iteration
            double R = VAL[M1] * VAL[M1];
            double S = O[M1] / (R + O[M1]);
            double T = ZERO;
            double U = S * (VAL[M1] + VAL[M1+1]);
            VAL[M1] = VAL[M1] + U;
            if (M1 != M2M1) {
                int M1P1 = M1 + 1;
                for (int I = M1P1; I <= M2M1; ++I) {
                    double G = VAL[I] - U;
                    R = (ONE - T) * O[I-1];
                    double ONES = ONE - S;
                    if (std::abs(ONES) > EPSI) R = G * G / ONES;
                    R = R + O[I];
                    O[I-1] = S * R;
                    if (O[I-1] == ZERO) MM = I;
                    T = S;
                    S = O[I] / R;
                    U = S * (G + VAL[I+1]);
                    VAL[I] = U + G;
                }
            }
            VAL[M2] = VAL[M2] - U;
            R = (ONE - T) * O[M2M1];
            double ONES = ONE - S;
            if (std::abs(ONES) > EPSI) R = VAL[M2] * VAL[M2] / ONES;
            O[M2M1] = S * R;

            // Shift back
            if (SHIFT != ZERO) {
                for (int I = M1; I <= M2; ++I) {
                    VAL[I] = VAL[I] + SHIFT;
                }
            }
            M1 = MM;
        }
        if (!converged) NEVER = NEVER + 1;

        // Convergence achieved
        if (M1 != M2M1) {
            M2   = M2M1;
            M2M1 = M2 - 1;
            continue;
        } else {
            if (M1 <= 2) break;
            M2 = M1 - 1;
            continue;
        }
    }

    // Reorder eigenvalues algebraically in ascending order
    if (N == 1) return;
    std::vector<double> Ocopy(N + 1);
    for (int K = 1; K <= N; ++K) {
        int M = -1;
        for (int m = 1; m <= N; ++m) {
            if (VAL[m] != -10000.0) { M = m; break; }
        }
        if (M == -1) M = N;
        if (M != N) {
            int MP1 = M + 1;
            for (int I = MP1; I <= N; ++I) {
                if (VAL[I] == -10000.0) continue;
                if (VAL[M] > VAL[I]) M = I;
            }
        }
        Ocopy[K] = VAL[M];
        VAL[M] = -10000.0;
    }
    for (int I = 1; I <= N; ++I) VAL[I] = Ocopy[I];
}

int main() {
    int N = 5;
    std::vector<double> VAL(N + 1);
    std::vector<double> O(N + 1);
    // Diagonal = 2, off-diagonal = -1 (O stores squares of off-diagonals)
    for (int i = 1; i <= N; ++i) VAL[i] = 2.0;
    for (int i = 1; i <= N - 1; ++i) O[i] = 1.0;  // (-1)^2
    O[N] = 0.0;

    qriter(VAL, O, N);

    std::vector<double> eig;
    for (int i = 1; i <= N; ++i) eig.push_back(VAL[i]);
    std::sort(eig.begin(), eig.end());

    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (size_t i = 0; i < eig.size(); ++i) {
        if (i) std::cout << ",";
        std::cout << eig[i];
    }
    std::cout << "]}" << std::endl;
    return 0;
}