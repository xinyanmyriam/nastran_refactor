#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

static void qriter2(std::vector<double>& val, std::vector<double>& o, int n) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE  = 1.0;
    int MAX = 100 * n;

    int M2 = n;
    bool done = false;

    while (!done) {
        int M2M1 = M2 - 1;
        int M1 = 0;
        bool allZero = true;
        for (int k = 1; k <= M2M1; ++k) {
            M1 = M2 - k;
            if (o[M1] != ZERO) { allZero = false; break; }
        }
        if (allZero) break;

        M2M1 = M1;
        M2   = M1 + 1;
        if (M2M1 != 1) {
            bool found = false;
            for (int k = 2; k <= M2M1; ++k) {
                M1 = M2 - k;
                if (o[M1] == ZERO) { found = true; break; }
            }
            if (found) M1 = M1 + 1;
        }
        int MM = M1;

        bool innerDone = false;
        while (!innerDone) {
            for (int ITER = 1; ITER <= MAX; ++ITER) {
                if (std::fabs(val[M2]) + o[M2M1] == std::fabs(val[M2])) break;
                bool allEq = true;
                for (int k = M1; k <= M2M1; ++k) {
                    if (val[k] != val[k+1]) { allEq = false; break; }
                }
                double SHIFT;
                if (allEq) {
                    SHIFT = ZERO;
                } else {
                    SHIFT = val[M2];
                    for (int i = M1; i <= M2M1; ++i) {
                        if (std::fabs(val[i]) < std::fabs(SHIFT)) SHIFT = val[i];
                    }
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
                        if (std::fabs(ONES) > EPSI) R = G * G / ONES;
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
                if (std::fabs(ONES) > EPSI) R = val[M2] * val[M2] / ONES;
                o[M2M1] = S * R;

                if (SHIFT != ZERO) {
                    for (int i = M1; i <= M2; ++i) val[i] = val[i] + SHIFT;
                }
                M1 = MM;
            }
            if (M1 == M2M1) {
                if (M1 <= 2) { done = true; innerDone = true; }
                else { M2 = M1 - 1; innerDone = true; }
            } else {
                M2   = M2M1;
                M2M1 = M2 - 1;
            }
        }
    }

    std::sort(val.begin() + 1, val.begin() + n + 1);
}

int main() {
    const int n = 5;
    std::vector<double> val(n + 1, 0.0);
    std::vector<double> o(n + 1, 0.0);
    for (int i = 1; i <= n; ++i) val[i] = 2.0;
    for (int i = 1; i <= n - 1; ++i) o[i] = 1.0;

    qriter2(val, o, n);

    std::cout << std::setprecision(15);
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 1; i <= n; ++i) {
        if (i > 1) std::cout << ",";
        std::cout << val[i];
    }
    std::cout << "]}" << std::endl;
    return 0;
}