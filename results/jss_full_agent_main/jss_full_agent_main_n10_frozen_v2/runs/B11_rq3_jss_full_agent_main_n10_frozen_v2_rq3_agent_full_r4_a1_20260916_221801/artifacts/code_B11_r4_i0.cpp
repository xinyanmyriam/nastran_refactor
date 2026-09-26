#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

void qriter(std::vector<double>& val, std::vector<double>& o, int n) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE  = 1.0;
    int MAX = 100 * n;
    auto VAL = [&](int i) -> double& { return val[i - 1]; };
    auto O   = [&](int i) -> double& { return o[i - 1]; };
    int M2 = n;
    while (true) {
        int M2M1 = M2 - 1;
        int M1 = 0;
        bool allZero = true;
        for (int k = 1; k <= M2M1; ++k) {
            M1 = M2 - k;
            if (O(M1) != ZERO) { allZero = false; break; }
        }
        if (allZero) break;
        M2M1 = M1;
        M2   = M1 + 1;
        if (M2M1 != 1) {
            bool found = false;
            for (int k = 2; k <= M2M1; ++k) {
                M1 = M2 - k;
                if (O(M1) == ZERO) { found = true; break; }
            }
            if (found) M1 = M1 + 1;
        }
        int MM = M1;
        bool converged = false;
        for (int ITER = 1; ITER <= MAX; ++ITER) {
            if (std::fabs(VAL(M2)) + O(M2M1) == std::fabs(VAL(M2))) { converged = true; break; }
            bool allEq = true;
            for (int k = M1; k <= M2M1; ++k) {
                if (VAL(k) != VAL(k + 1)) { allEq = false; break; }
            }
            double SHIFT;
            if (allEq) { SHIFT = ZERO; }
            else {
                SHIFT = VAL(M2);
                for (int i = M1; i <= M2M1; ++i) {
                    if (std::fabs(VAL(i)) < std::fabs(SHIFT)) SHIFT = VAL(i);
                }
                for (int i = M1; i <= M2; ++i) VAL(i) = VAL(i) - SHIFT;
            }
            double R = VAL(M1) * VAL(M1);
            double S = O(M1) / (R + O(M1));
            double T = ZERO;
            double U = S * (VAL(M1) + VAL(M1 + 1));
            VAL(M1) = VAL(M1) + U;
            if (M1 != M2M1) {
                int M1P1 = M1 + 1;
                for (int i = M1P1; i <= M2M1; ++i) {
                    double G = VAL(i) - U;
                    R = (ONE - T) * O(i - 1);
                    double ONES = ONE - S;
                    if (std::fabs(ONES) > EPSI) R = G * G / ONES;
                    R = R + O(i);
                    O(i - 1) = S * R;
                    if (O(i - 1) == ZERO) MM = i;
                    T = S;
                    S = O(i) / R;
                    U = S * (G + VAL(i + 1));
                    VAL(i) = U + G;
                }
            }
            VAL(M2) = VAL(M2) - U;
            R = (ONE - T) * O(M2M1);
            double ONES = ONE - S;
            if (std::fabs(ONES) > EPSI) R = VAL(M2) * VAL(M2) / ONES;
            O(M2M1) = S * R;
            if (SHIFT != ZERO) {
                for (int i = M1; i <= M2; ++i) VAL(i) = VAL(i) + SHIFT;
            }
            M1 = MM;
        }
        if (M1 != M2M1) { M2 = M2M1; M2M1 = M2 - 1; continue; }
        else { if (M1 <= 2) break; M2 = M1 - 1; continue; }
    }
    if (n == 1) return;
    std::vector<double> tmp(n);
    for (int k = 1; k <= n; ++k) {
        int M = -1;
        for (int m = 1; m <= n; ++m) { if (VAL(m) != -10000.0) { M = m; break; } }
        if (M == -1) M = n;
        if (M != n) {
            int MP1 = M + 1;
            for (int i = MP1; i <= n; ++i) {
                if (VAL(i) == -10000.0) continue;
                if (VAL(M) > VAL(i)) M = i;
            }
        }
        tmp[k - 1] = VAL(M);
        VAL(M) = -10000.0;
    }
    for (int i = 1; i <= n; ++i) VAL(i) = tmp[i - 1];
}

int main() {
    int n = 5;
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> o   = {-1.0, -1.0, -1.0, -1.0};
    qriter(val, o, n);
    std::vector<double> eig(val.begin(), val.end());
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