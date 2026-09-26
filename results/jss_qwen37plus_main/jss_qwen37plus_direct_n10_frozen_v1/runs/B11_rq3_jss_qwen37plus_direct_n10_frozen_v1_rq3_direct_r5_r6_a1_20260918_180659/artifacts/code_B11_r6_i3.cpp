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
        // Find largest k < m such that e[k-1] == 0 (deflation)
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
            double sqrt_disc = std::sqrt(std::max(discriminant, 0.0));
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
        // Start with f = d[l], g = e[l]
        double f = d[l];
        double g = e[l];
        double c_val = 1.0, s_val = 0.0;

        // We'll declare c_val/s_val outside loop so they're available after
        for (int i = l; i <= m-2; ++i) {
            // Compute Givens rotation to zero out g (i.e., eliminate subdiagonal element)
            double r = std::sqrt(f*f + g*g);
            if (r == 0.0) {
                r = EPS;
            }
            c_val = f / r;
            s_val = g / r;

            // Apply rotation to columns i and i+1 (affects d[i], e[i], d[i+1])
            // First, update row i (d[i] and e[i])
            if (i > l) {
                e[i-1] = s_val * d[i-1];
            }
            double d_i_prev = d[i-1];
            d[i-1] = c_val * d[i-1] + s_val * e[i];
            // e[i] becomes: c_val * e[i] - s_val * d[i] — but d[i] not yet updated; save old d[i]
            double old_d_i = d[i];
            d[i] = c_val * old_d_i + s_val * e[i];
            e[i] = c_val * e[i] - s_val * old_d_i;

            // Now update row i+1: affects e[i] (already done) and d[i+1]
            // But note: e[i] is now the (i,i+1) element, and d[i+1] is diagonal
            // Standard update: 
            //   [d[i]   e[i] ]   [c s]   -> new row i
            //   [e[i]  d[i+1]]   [-s c]
            // So new d[i+1] = c*d[i+1] - s*e[i]
            // And new e[i] was already computed above.
            double old_d_ip1 = d[i+1];
            d[i+1] = c_val * old_d_ip1 - s_val * e[i];
            // Also update e[i] again? No — it's already updated as off-diagonal.
            // But we need to update e[i] for next iteration? Actually, e[i] is now the (i,i+1) entry,
            // and we will use e[i+1] next — so we must preserve e[i+1] until used.

            // Prepare for next iteration: f = d[i+1], g = e[i+1] (if exists)
            if (i < m-2) {
                f = d[i+1];
                g = e[i+1];
            }
        }

        // After loop, apply final rotation to last pair: (d[m-2], e[m-2]) and d[m-1]
        // At this point, f and g correspond to the last active pair before m-1
        // But our loop ended at i = m-2, so we have c_val/s_val from last iteration
        // Now apply to d[m-2] and d[m-1] using e[m-2] (which was updated in loop)
        // Actually, the standard step sets:
        //   e[m-2] = s_val * d[m-1]
        //   d[m-1] = c_val * d[m-1]
        // BUT: careful — d[m-1] hasn't been rotated yet in the last step.
        // Let's reconstruct: after loop, the last rotation was applied to indices (m-2, m-1)
        // So d[m-1] was updated as: d[m-1] = c_val * d[m-1] - s_val * e[m-2] (from above logic)
        // However, the classical algorithm does *not* do that — instead, the final rotation
        // acts on the last 2x2 block: [d[m-2], e[m-2]; e[m-2], d[m-1]]
        // So we should apply the final rotation explicitly to that block.

        // Instead, simpler & correct: after the loop, we have updated up to d[m-1] implicitly.
        // But our loop already handled i = m-2, so d[m-1] was updated in the body.
        // So we skip redundant final assignment — and instead ensure convergence check uses latest e[m-2].

        // Unshift
        for (int i = l; i <= m-1; ++i) {
            d[i] += mu;
        }

        // Check convergence of e[m-2]
        if (m >= 2 && std::abs(e[m-2]) <= EPS * (std::abs(d[m-2]) + std::abs(d[m-1]))) {
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