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
    
    while (m >= 0) {
        // Deflation check: if off_diag[m-1] is negligible, eigenvalue found
        if (m == 0) {
            // 1x1 block
            break;
        }
        
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

        // If l == m, then we have a 1x1 block at position m
        if (l == m) {
            m--;
            continue;
        }

        // Wilkinson shift: eigenvalue of bottom 2x2 submatrix closest to diag[m]
        double d_m = diag[m];
        double d_m1 = diag[m-1];
        double e_m1 = off_diag[m-1];
        double mu;
        
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

        // Apply shift to [l..m]
        for (int i = l; i <= m; ++i) {
            diag[i] -= mu;
        }

        // QL step (chase bulge from top to bottom)
        // Start with first rotation on rows l and l+1
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

            // Apply rotation to columns k and k+1 (affects rows k and k+1)
            // Update diagonal[k]
            double temp_diag_k = c * c * diag[k] + 2.0 * c * s * off_diag[k] + s * s * diag[k+1];
            // Update off-diagonal[k]
            double temp_off_k = c * s * (diag[k+1] - diag[k]) + (c*c - s*s) * off_diag[k];
            // Update diagonal[k+1]
            double temp_diag_k1 = s * s * diag[k] - 2.0 * c * s * off_diag[k] + c * c * diag[k+1];
            
            // But standard approach: apply rotation to eliminate subdiagonal
            // Better to use standard bulge chasing for tridiagonal
            
            // Standard approach: rotate rows k and k+1
            // First, save current values
            double d_k = diag[k];
            double d_k1 = diag[k+1];
            double e_k = (k < m-1) ? off_diag[k+1] : 0.0;
            
            // Apply rotation to eliminate off_diag[k] (make it zero)
            // New diag[k] = c*d_k + s*off_diag[k]
            // New off_diag[k] = -s*d_k + c*off_diag[k]
            // New diag[k+1] = c*d_k1 + s*e_k
            // New off_diag[k-1] gets updated if k > l
            
            // Actually, let's implement the standard QL step correctly:
            // We want to apply Givens rotation G(k,k+1) to eliminate off_diag[k]
            double old_off_k = off_diag[k];
            
            // Update diag[k] and off_diag[k]
            double new_diag_k = c * d_k + s * old_off_k;
            double new_off_k = -s * d_k + c * old_off_k;
            
            // Update off_diag[k-1] if applicable
            if (k > l) {
                off_diag[k-1] = c * off_diag[k-1] - s * x;
            }
            
            // Update diag[k+1] and off_diag[k+1]
            double new_diag_k1 = c * d_k1;
            if (k < m-1) {
                new_diag_k1 += s * off_diag[k+1];
            }
            
            // Store updates
            diag[k] = new_diag_k;
            off_diag[k] = new_off_k;
            diag[k+1] = new_diag_k1;
            
            // Prepare for next rotation
            x = new_off_k;
            if (k < m-1) {
                z = off_diag[k+1];
            } else {
                z = 0.0;
            }
        }

        // Undo shift
        for (int i = l; i <= m; ++i) {
            diag[i] += mu;
        }

        // Check for deflation at bottom
        if (m > 0) {
            tol = eps * (std::abs(diag[m]) + std::abs(diag[m-1]));
            if (std::abs(off_diag[m-1]) <= tol) {
                off_diag[m-1] = 0.0;
                m--;
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