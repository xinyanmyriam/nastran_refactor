#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Since we can't use external JSON library, we'll create a simple JSON string builder
std::string to_json_string(const std::vector<double>& eigenvalues) {
    std::ostringstream oss;
    oss << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (size_t i = 0; i < eigenvalues.size(); ++i) {
        if (i > 0) oss << ",";
        oss << std::fixed << std::setprecision(12) << eigenvalues[i];
    }
    oss << "]}";
    return oss.str();
}

// QR/QL iteration for symmetric tridiagonal matrix
// Implements stable implicit QL algorithm with Wilkinson shift and deflation
void qriter(std::vector<double>& diag, std::vector<double>& off_diag) {
    const int n = static_cast<int>(diag.size());
    if (n <= 0 || off_diag.size() != static_cast<size_t>(n-1)) return;

    const double eps = 1e-15;
    const int max_iter = 30; // per eigenvalue — standard safe bound

    int m = n - 1; // current bottom index of active submatrix [0..m]
    
    while (m > 0) {
        // Deflation check: if off_diag[m-1] is negligible, eigenvalue found
        double tol = eps * (std::abs(diag[m]) + std::abs(diag[m-1]));
        if (std::abs(off_diag[m-1]) <= tol) {
            off_diag[m-1] = 0.0;
            m--;
            continue;
        }

        // Look for consecutive zeros in off_diag to split problem
        int l = m;
        while (l > 0) {
            tol = eps * (std::abs(diag[l]) + std::abs(diag[l-1]));
            if (std::abs(off_diag[l-1]) <= tol) {
                off_diag[l-1] = 0.0;
                break;
            }
            l--;
        }

        // If l == 0, entire matrix is active; else work on [l..m]
        if (l == m) {
            // 1x1 block — done
            m--;
            continue;
        }

        // Wilkinson shift: eigenvalue of bottom 2x2 submatrix closest to diag[m]
        double d_m = diag[m];
        double d_m1 = diag[m-1];
        double e_m1 = off_diag[m-1];
        double mu;
        if (l == m-1) {
            // 2x2 case: use exact eigenvalues
            double tr = d_m1 + d_m;
            double det = d_m1 * d_m - e_m1 * e_m1;
            double disc = tr * tr - 4.0 * det;
            if (disc >= 0.0) {
                double sqrt_disc = std::sqrt(disc);
                double lambda1 = 0.5 * (tr - sqrt_disc);
                double lambda2 = 0.5 * (tr + sqrt_disc);
                mu = (std::abs(lambda1 - d_m) < std::abs(lambda2 - d_m)) ? lambda1 : lambda2;
            } else {
                mu = d_m; // fallback
            }
        } else {
            // Wilkinson shift: use bottom 2x2
            double a = d_m1, b = e_m1, c = e_m1, d = d_m;
            double tr = a + d;
            double det = a * d - b * c;
            double disc = tr * tr - 4.0 * det;
            if (disc >= 0.0) {
                double sqrt_disc = std::sqrt(disc);
                double lambda1 = 0.5 * (tr - sqrt_disc);
                double lambda2 = 0.5 * (tr + sqrt_disc);
                mu = (std::abs(lambda1 - d_m) < std::abs(lambda2 - d_m)) ? lambda1 : lambda2;
            } else {
                mu = d_m;
            }
        }

        // Apply shift to [l..m]
        for (int i = l; i <= m; ++i) {
            diag[i] -= mu;
        }

        // QL step (chase bulge from top to bottom)
        double x = diag[l];
        double z = off_diag[l];
        for (int k = l; k < m; ++k) {
            // Compute Givens rotation for (x, z)
            double r = std::sqrt(x*x + z*z);
            if (r == 0.0) {
                r = eps;
            }
            double c = x / r;
            double s = z / r;

            // Apply to rows k and k+1
            double temp = c * diag[k] + s * off_diag[k];
            off_diag[k] = -s * diag[k] + c * off_diag[k];
            diag[k] = temp;

            if (k > l) {
                off_diag[k-1] = c * off_diag[k-1] - s * x;
            }

            x = off_diag[k];
            if (k < m-1) {
                z = off_diag[k+1];
            } else {
                z = 0.0;
            }

            // Update diag[k+1]
            temp = c * diag[k+1];
            if (k < m-1) {
                temp += s * off_diag[k+1];
            }
            diag[k+1] = temp;
        }

        // Undo shift
        for (int i = l; i <= m; ++i) {
            diag[i] += mu;
        }

        // Check for deflation at bottom
        tol = eps * (std::abs(diag[m]) + std::abs(diag[m-1]));
        if (std::abs(off_diag[m-1]) <= tol) {
            off_diag[m-1] = 0.0;
            m--;
        } else {
            // Try again — but limit iterations to avoid infinite loops
            static int iter = 0;
            iter++;
            if (iter > max_iter * n) {
                break;
            }
        }
    }
}

// Reorder eigenvalues in ascending order using std::sort
void reorder_eigenvalues(std::vector<double>& eigenvalues) {
    std::sort(eigenvalues.begin(), eigenvalues.end());
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    std::vector<double> diag = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> off_diag = {-1.0, -1.0, -1.0, -1.0}; // size 4

    // Run QR/QL iteration
    qriter(diag, off_diag);

    // Reorder eigenvalues in ascending order
    reorder_eigenvalues(diag);

    // Output as JSON
    std::cout << to_json_string(diag) << std::endl;

    return 0;
}