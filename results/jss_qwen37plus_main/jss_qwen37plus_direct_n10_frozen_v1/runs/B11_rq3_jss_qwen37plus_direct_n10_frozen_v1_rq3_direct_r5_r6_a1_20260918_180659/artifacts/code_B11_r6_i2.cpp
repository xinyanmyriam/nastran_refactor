#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Since we cannot use external JSON library, we'll create a simple JSON string builder
std::string to_json_string(const std::vector<double>& eigenvalues) {
    std::ostringstream oss;
    oss << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        if (i > 0) oss << ",";
        oss << std::setprecision(15) << eigenvalues[i];
    }
    oss << "]}";
    return oss.str();
}

// QR iteration for symmetric tridiagonal matrix eigenvalues (implicit shifted QR with Wilkinson shift)
void qriter(std::vector<double>& diag, std::vector<double>& off_diag) {
    const int n = static_cast<int>(diag.size());
    if (n <= 0) return;

    // Work on copies to avoid modifying inputs until final sort
    std::vector<double> d = diag;
    std::vector<double> e = off_diag; // off-diagonal elements (length n-1)

    const double EPS = 1e-12;
    const int MAX_ITER = 30 * n;

    int iter = 0;
    int m = n;

    while (m > 0 && iter < MAX_ITER) {
        // Find largest k < m such that e[k] == 0 (deflation)
        int k = m - 1;
        while (k > 0 && std::abs(e[k-1]) > EPS * (std::abs(d[k-1]) + std::abs(d[k]))) {
            k--;
        }

        if (k == 0) {
            // No deflation possible; entire matrix active
            k = m;
        } else {
            // Deflate at position k-1: set e[k-1] = 0
            e[k-1] = 0.0;
        }

        // If only one element left, done
        if (k == m) {
            m--;
            continue;
        }

        // Now work on submatrix [k, m-1] (0-indexed, inclusive)
        int l = k;
        while (l > 0 && std::abs(e[l-1]) > EPS * (std::abs(d[l-1]) + std::abs(d[l]))) {
            l--;
        }

        // Active submatrix is [l, m-1]
        if (l == m - 1) {
            // 1x1 block
            m--;
            continue;
        }

        // Wilkinson shift: eigenvalue of bottom-right 2x2 submatrix closest to d[m-1]
        double a = d[m-2], b = d[m-1], c = e[m-2];
        double mu;
        if (m - l == 2) {
            double discriminant = (a - b) * (a - b) + 4.0 * c * c;
            double sqrt_disc = std::sqrt(discriminant);
            double lambda1 = 0.5 * (a + b - sqrt_disc);
            double lambda2 = 0.5 * (a + b + sqrt_disc);
            mu = (std::abs(lambda1 - b) < std::abs(lambda2 - b)) ? lambda1 : lambda2;
        } else {
            mu = d[m-1];
        }

        // Apply shift to active submatrix [l, m-1]
        for (int i = l; i <= m-1; ++i) {
            d[i] -= mu;
        }

        // Bulge chase (Givens rotations) — implicit QR step
        double f = d[l];
        double g = e[l];
        for (int i = l; i <= m-2; ++i) {
            // Compute Givens rotation to zero out f
            double r = std::sqrt(f*f + g*g);
            if (r == 0.0) {
                r = EPS;
            }
            double c_val = f / r;
            double s_val = g / r;

            // Apply to (d[i], e[i])
            if (i > l) {
                e[i-1] = s_val * d[i-1];
            }
            d[i] = c_val * d[i-1] + s_val * e[i];
            e[i] = c_val * e[i] - s_val * d[i];

            if (i < m-2) {
                f = d[i];
                g = e[i+1];
                // Rotate next pair
                d[i] = c_val * d[i] + s_val * e[i+1];
                e[i+1] = c_val * e[i+1] - s_val * d[i];
                f = d[i];
                g = e[i+1];
            }
        }
        e[m-2] = s_val * d[m-1];
        d[m-1] = c_val * d[m-1];

        // Unshift
        for (int i = l; i <= m-1; ++i) {
            d[i] += mu;
        }

        // Check convergence of e[m-2]
        if (std::abs(e[m-2]) <= EPS * (std::abs(d[m-2]) + std::abs(d[m-1]))) {
            e[m-2] = 0.0;
            m--;
        }

        iter++;
    }

    // Copy computed eigenvalues into diag
    diag = d;

    // Sort eigenvalues in ascending order
    std::sort(diag.begin(), diag.end());
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    std::vector<double> diag = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> off_diag = {-1.0, -1.0, -1.0, -1.0}; // 4 elements for 5x5

    // Run QR iteration
    qriter(diag, off_diag);

    // Output JSON
    std::cout << to_json_string(diag) << std::endl;

    return 0;
}