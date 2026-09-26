#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>

std::vector<double> qriter(std::vector<double> d, std::vector<double> e) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE  = 1.0;

    int n = (int)d.size();
    std::vector<double> O(n, 0.0);
    for (int i = 1; i <= n - 1; ++i) O[i] = e[i - 1] * e[i - 1];

    std::vector<double> VAL(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) VAL[i] = d[i - 1];

    int MAX = 100 * n;

    int M2 = n;
    bool done = false;

    while (!done) {
        int M2M1 = M2 - 1;
        int M1 = 0;
        bool allZero = true;
        for (int k = 1; k <= M2M1; ++k) {
            M1 = M2 - k;
            if (O[M1] != ZERO) { allZero = false; break; }
        }
        if (allZero) break;
        M2M1 = M1;
        M2 = M1 + 1;
        if (M2M1 != 1) {
            bool found = false;
            for (int k = 2; k <= M2M1; ++k) {
                M1 = M2 - k;
                if (O[M1] == ZERO) { found = true; break; }
            }
            if (found) M1 = M1 + 1;
        }
        int MM = M1;

        for (int ITER = 1; ITER <= MAX; ++ITER) {
            if (std::fabs(VAL[M2]) + O[M2M1] == std::fabs(VAL[M2])) break;
            bool allEqual = true;
            for (int k = M1; k <= M2M1; ++k) {
                if (VAL[k] != VAL[k + 1]) { allEqual = false; break; }
            }
            double SHIFT = ZERO;
            if (!allEqual) {
                SHIFT = VAL[M2];
                for (int i = M1; i <= M2M1; ++i) {
                    if (std::fabs(VAL[i]) < std::fabs(SHIFT)) SHIFT = VAL[i];
                }
                for (int i = M1; i <= M2; ++i) VAL[i] -= SHIFT;
            }
            double R = VAL[M1] * VAL[M1];
            double S = O[M1] / (R + O[M1]);
            double T = ZERO;
            double U = S * (VAL[M1] + VAL[M1 + 1]);
            VAL[M1] = VAL[M1] + U;
            if (M1 != M2M1) {
                int M1P1 = M1 + 1;
                for (int i = M1P1; i <= M2M1; ++i) {
                    double G = VAL[i] - U;
                    R = (ONE - T) * O[i - 1];
                    double ONES = ONE - S;
                    if (std::fabs(ONES) > EPSI) R = G * G / ONES;
                    R = R + O[i];
                    O[i - 1] = S * R;
                    if (O[i - 1] == ZERO) MM = i;
                    T = S;
                    S = O[i] / R;
                    U = S * (G + VAL[i + 1]);
                    VAL[i] = U + G;
                }
            }
            VAL[M2] = VAL[M2] - U;
            R = (ONE - T) * O[M2M1];
            double ONES = ONE - S;
            if (std::fabs(ONES) > EPSI) R = VAL[M2] * VAL[M2] / ONES;
            O[M2M1] = S * R;
            if (SHIFT != ZERO) {
                for (int i = M1; i <= M2; ++i) VAL[i] += SHIFT;
            }
            M1 = MM;
        }
        if (M1 == M2M1) {
            if (M1 <= 2) { done = true; break; }
            M2 = M1 - 1;
        } else {
            M2 = M2M1;
            M2M1 = M2 - 1;
        }
    }

    std::vector<double> ev(n);
    for (int i = 0; i < n; ++i) ev[i] = VAL[i + 1];
    std::sort(ev.begin(), ev.end());
    return ev;
}

int main() {
    int n = 5;
    std::vector<double> diag(n, 2.0);
    std::vector<double> off(n - 1, -1.0);
    std::vector<double> ev = qriter(diag, off);
    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (size_t i = 0; i < ev.size(); ++i) {
        if (i) std::cout << ",";
        std::cout << ev[i];
    }
    std::cout << "]}" << std::endl;
    return 0;
}