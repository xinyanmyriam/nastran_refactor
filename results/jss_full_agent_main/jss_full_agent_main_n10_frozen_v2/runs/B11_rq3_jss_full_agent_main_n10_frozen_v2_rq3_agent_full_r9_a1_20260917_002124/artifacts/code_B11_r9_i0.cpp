#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

static void qriter(std::vector<double>& VAL, std::vector<double>& O, int n) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE  = 1.0;
    int MAX = 100 * n;
    int M2 = n;
    while (true) {
        int M2M1 = M2 - 1;
        int M1 = 0;
        bool allZero = true;
        for (int K = 1; K <= M2M1; ++K) {
            M1 = M2 - K;
            if (O[M1] != ZERO) { allZero = false; break; }
        }
        if (allZero) break;
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
        int MM = M1;
        for (int ITER = 1; ITER <= MAX; ++ITER) {
            if (std::abs(VAL[M2]) + O[M2M1] == std::abs(VAL[M2])) break;
            bool allEq = true;
            for (int K = M1; K <= M2M1; ++K) {
                if (VAL[K] != VAL[K+1]) { allEq = false; break; }
            }
            double SHIFT = ZERO;
            if (!allEq) {
                SHIFT = VAL[M2];
                for (int I = M1; I <= M2M1; ++I)
                    if (std::abs(VAL[I]) < std::abs(SHIFT)) SHIFT = VAL[I];
                for (int I = M1; I <= M2; ++I) VAL[I] -= SHIFT;
            }
            double R = VAL[M1] * VAL[M1];
            double S = O[M1] / (R + O[M1]);
            double T = ZERO;
            double U = S * (VAL[M1] + VAL[M1+1]);
            VAL[M1] += U;
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
            VAL[M2] -= U;
            R = (ONE - T) * O[M2M1];
            double ONES = ONE - S;
            if (std::abs(ONES) > EPSI) R = VAL[M2] * VAL[M2] / ONES;
            O[M2M1] = S * R;
            if (SHIFT != ZERO)
                for (int I = M1; I <= M2; ++I) VAL[I] += SHIFT;
            M1 = MM;
        }
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
    std::vector<double> tmp(n + 1);
    for (int K = 1; K <= n; ++K) {
        int M = 0;
        for (int m = 1; m <= n; ++m)
            if (VAL[m] != -10000.0) { M = m; break; }
        if (M != n) {
            for (int I = M + 1; I <= n; ++I) {
                if (VAL[I] == -10000.0) continue;
                if (VAL[M] > VAL[I]) M = I;
            }
        }
        tmp[K] = VAL[M];
        VAL[M] = -10000.0;
    }
    for (int I = 1; I <= n; ++I) VAL[I] = tmp[I];
}

int main() {
    int n = 5;
    std::vector<double> VAL(n + 1), O(n + 1);
    for (int i = 1; i <= n; ++i) VAL[i] = 2.0;
    for (int i = 1; i <= n - 1; ++i) O[i] = 1.0;
    qriter(VAL, O, n);
    std::vector<double> eigs;
    for (int i = 1; i <= n; ++i) eigs.push_back(VAL[i]);
    std::sort(eigs.begin(), eigs.end());
    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (size_t i = 0; i < eigs.size(); ++i) {
        if (i) std::cout << ",";
        std::cout << eigs[i];
    }
    std::cout << "]}" << std::endl;
    return 0;
}