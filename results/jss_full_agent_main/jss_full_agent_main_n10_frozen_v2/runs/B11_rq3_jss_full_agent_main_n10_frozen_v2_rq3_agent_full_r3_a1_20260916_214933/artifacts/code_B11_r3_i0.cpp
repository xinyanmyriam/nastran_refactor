#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

static void qriter(std::vector<double>& val, std::vector<double>& o, int n) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE  = 1.0;
    int MAX = 100 * n;
    int M2 = n;
    while (true) {
        int M2M1 = M2 - 1;
        int M1 = 0;
        bool found = false;
        for (int k = 1; k <= M2M1; ++k) {
            M1 = M2 - k;
            if (o[M1] != ZERO) { found = true; break; }
        }
        if (!found) break;
        M2M1 = M1;
        M2   = M1 + 1;
        if (M2M1 != 1) {
            bool found2 = false;
            for (int k = 2; k <= M2M1; ++k) {
                M1 = M2 - k;
                if (o[M1] == ZERO) { found2 = true; break; }
            }
            if (found2) M1 = M1 + 1;
        }
        int MM = M1;
        bool converged = false;
        for (int ITER = 1; ITER <= MAX; ++ITER) {
            if (std::abs(val[M2]) + o[M2M1] == std::abs(val[M2])) { converged = true; break; }
            bool allEqual = true;
            for (int k = M1; k <= M2M1; ++k) {
                if (val[k] != val[k+1]) { allEqual = false; break; }
            }
            double SHIFT = ZERO;
            if (!allEqual) {
                SHIFT = val[M2];
                for (int i = M1; i <= M2M1; ++i)
                    if (std::abs(val[i]) < std::abs(SHIFT)) SHIFT = val[i];
                for (int i = M1; i <= M2; ++i) val[i] = val[i] - SHIFT;
            }
            double R = val[M1] * val[M1];
            double S = o[M1] / (R + o[M1]);
            double T = ZERO;
            double U = S * (val[M1] + val[M1+1]);
            val[M1] = val[M1] + U;
            if (M1 != M2M1) {
                int M1P1 = M1 + 1;
                for (int i = M1P1; i <= M2M1; ++i) {
                    double G = val[i] - U;
                    R = (ONE - T) * o[i-1];
                    double ONES = ONE - S;
                    if (std::abs(ONES) > EPSI) R = G * G / ONES;
                    R = R + o[i];
                    o[i-1] = S * R;
                    if (o[i-1] == ZERO) MM = i;
                    T = S;
                    S = o[i] / R;
                    U = S * (G + val[i+1]);
                    val[i] = U + G;
                }
            }
            val[M2] = val[M2] - U;
            R = (ONE - T) * o[M2M1];
            double ONES = ONE - S;
            if (std::abs(ONES) > EPSI) R = val[M2] * val[M2] / ONES;
            o[M2M1] = S * R;
            if (SHIFT != ZERO) {
                for (int i = M1; i <= M2; ++i) val[i] = val[i] + SHIFT;
            }
            M1 = MM;
        }
        if (converged) {
            if (M1 != M2M1) { M2 = M2M1; M2M1 = M2 - 1; continue; }
            else { if (M1 <= 2) break; M2 = M1 - 1; continue; }
        } else break;
    }
    if (n == 1) return;
    std::vector<double> tmp(n + 1);
    for (int k = 1; k <= n; ++k) {
        int M = 0;
        for (int m = 1; m <= n; ++m) if (val[m] != -10000.0) { M = m; break; }
        if (M == 0) M = n;
        if (M != n) {
            for (int i = M + 1; i <= n; ++i) {
                if (val[i] == -10000.0) continue;
                if (val[M] > val[i]) M = i;
            }
        }
        tmp[k] = val[M];
        val[M] = -10000.0;
    }
    for (int i = 1; i <= n; ++i) val[i] = tmp[i];
}

int main() {
    const int n = 5;
    std::vector<double> val(n + 1, 0.0), o(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) val[i] = 2.0;
    for (int i = 1; i <= n - 1; ++i) o[i] = 1.0;
    qriter(val, o, n);
    std::vector<double> eig;
    for (int i = 1; i <= n; ++i) eig.push_back(val[i]);
    std::sort(eig.begin(), eig.end());
    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (size_t i = 0; i < eig.size(); ++i) { if (i) std::cout << ","; std::cout << eig[i]; }
    std::cout << "]}" << std::endl;
    return 0;
}