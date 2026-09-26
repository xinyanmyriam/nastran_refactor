#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Forward declarations for Eigen types
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Constants
const double EPSI = 1.0e-10;
const double ZERO = 0.0;
const double ONE = 1.0;

// QR iteration for tridiagonal matrix eigenvalues
// Input: diagonal vector 'd' and off-diagonal vector 'e' (size n-1)
// Output: eigenvalues in sorted ascending order
std::vector<double> qriter(const std::vector<double>& d_in, const std::vector<double>& e_in) {
    int n = static_cast<int>(d_in.size());
    if (n == 0) return {};
    if (n == 1) return {d_in[0]};
    
    // Copy input to working arrays
    std::vector<double> d = d_in;
    std::vector<double> e = e_in; // e[i] is the off-diagonal between d[i] and d[i+1]
    
    // Main iteration loop
    int m = n - 1; // current size of active submatrix
    
    while (m > 0) {
        // Check for convergence: if e[m-1] is negligible, we have an eigenvalue
        if (std::abs(e[m-1]) <= EPSI * (std::abs(d[m-1]) + std::abs(d[m]))) {
            m--;
            continue;
        }
        
        // Find the largest unreduced submatrix
        int l = m;
        while (l > 0 && std::abs(e[l-1]) > EPSI * (std::abs(d[l-1]) + std::abs(d[l]))) {
            l--;
        }
        
        // Compute Wilkinson shift from bottom-right 2x2 submatrix
        double s = ZERO;
        if (l < m) {
            // 2x2 submatrix: [d[m-1], e[m-1]; e[m-1], d[m]]
            double a = d[m-1];
            double b = e[m-1];
            double c = d[m];
            
            double trace = a + c;
            double det = a * c - b * b;
            double sqrt_disc = std::sqrt(trace * trace - 4.0 * det);
            
            // Choose eigenvalue closer to d[m] as shift
            double lambda1 = 0.5 * (trace + sqrt_disc);
            double lambda2 = 0.5 * (trace - sqrt_disc);
            
            if (std::abs(lambda1 - d[m]) < std::abs(lambda2 - d[m])) {
                s = lambda1;
            } else {
                s = lambda2;
            }
        } else {
            // Single element, shift is just d[l]
            s = d[l];
        }
        
        // Apply implicit QR step with shift s
        // First, compute the first Householder reflection (or Givens) to introduce bulge
        double x = d[l] - s;
        double y = (l < m) ? e[l] : ZERO;
        
        // Bulge chasing: apply sequence of Givens rotations
        for (int k = l; k < m; ++k) {
            // Compute Givens rotation to zero out e[k] in the shifted matrix
            double c, s_rot;
            double r;
            
            if (k == l) {
                // First rotation: act on [x, y]^T where x = d[l]-s, y = e[l]
                r = std::sqrt(x*x + y*y);
                if (r == ZERO) {
                    c = ONE;
                    s_rot = ZERO;
                } else {
                    c = x / r;
                    s_rot = y / r;
                }
                
                // Apply to first row/column
                if (k < m) {
                    double temp = c * d[k] + s_rot * e[k];
                    e[k] = -s_rot * d[k] + c * e[k];
                    d[k] = temp;
                }
                
                if (k + 1 < n) {
                    double temp = c * e[k] + s_rot * d[k+1];
                    d[k+1] = -s_rot * e[k] + c * d[k+1];
                    e[k] = temp;
                }
            } else {
                // Subsequent rotations: act on [e[k-1], d[k+1]] or similar
                // Standard bulge chasing: apply Givens to eliminate bulge
                double a = e[k-1];
                double b = d[k+1];
                
                r = std::sqrt(a*a + b*b);
                if (r == ZERO) {
                    c = ONE;
                    s_rot = ZERO;
                } else {
                    c = a / r;
                    s_rot = b / r;
                }
                
                // Apply to submatrix
                if (k < m) {
                    double temp = c * d[k] + s_rot * e[k];
                    e[k] = -s_rot * d[k] + c * e[k];
                    d[k] = temp;
                }
                
                if (k + 1 < n) {
                    double temp = c * e[k] + s_rot * d[k+1];
                    d[k+1] = -s_rot * e[k] + c * d[k+1];
                    e[k] = temp;
                }
            }
        }
        
        // After bulge chasing, check if we can deflate
        if (std::abs(e[m-1]) <= EPSI * (std::abs(d[m-1]) + std::abs(d[m]))) {
            m--;
        }
    }
    
    // Simple convergence check and cleanup
    // In practice, we might need more iterations, but for this test case,
    // let's use a more robust approach with explicit QR steps
    
    // Reset and use standard implicit QR with better convergence
    d = d_in;
    e = e_in;
    m = n - 1;
    
    const int max_iter = 30;
    for (int iter = 0; iter < max_iter && m > 0; ++iter) {
        // Look for deflation from bottom
        int l = m;
        while (l > 0 && std::abs(e[l-1]) > EPSI * (std::abs(d[l-1]) + std::abs(d[l]))) {
            l--;
        }
        
        if (l == m) {
            // Isolated eigenvalue
            m--;
            continue;
        }
        
        // Compute Wilkinson shift from bottom 2x2
        double a = d[m-1];
        double b = e[m-1];
        double c = d[m];
        double trace = a + c;
        double det = a * c - b * b;
        double sqrt_disc = std::sqrt(std::max(0.0, trace * trace - 4.0 * det));
        double lambda1 = 0.5 * (trace + sqrt_disc);
        double lambda2 = 0.5 * (trace - sqrt_disc);
        double s = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
        
        // Implicit QR step: chase bulge from top to bottom
        double x = d[l] - s;
        double y = (l < n-1) ? e[l] : ZERO;
        
        for (int k = l; k < m; ++k) {
            // Compute Givens rotation to zero the subdiagonal element
            double r = std::sqrt(x*x + y*y);
            if (r == ZERO) {
                break;
            }
            double c_rot = x / r;
            double s_rot = y / r;
            
            // Apply to rows k and k+1
            if (k == l) {
                // First rotation affects d[l] and e[l]
                double temp_d = c_rot * d[k] + s_rot * e[k];
                double temp_e = -s_rot * d[k] + c_rot * e[k];
                d[k] = temp_d;
                e[k] = temp_e;
            }
            
            // Apply to column k+1
            if (k + 1 < n) {
                double temp_d_next = c_rot * e[k] + s_rot * d[k+1];
                double temp_e_next = -s_rot * e[k] + c_rot * d[k+1];
                if (k < n-1) e[k] = temp_d_next;
                if (k + 1 < n) d[k+1] = temp_e_next;
            }
            
            // Prepare for next rotation
            if (k + 1 < m) {
                x = e[k];
                y = (k + 1 < n-1) ? e[k+1] : ZERO;
            }
        }
        
        // Check for deflation
        if (std::abs(e[m-1]) <= EPSI * (std::abs(d[m-1]) + std::abs(d[m]))) {
            m--;
        }
    }
    
    // Return eigenvalues (they should be in d, but sort to be safe)
    std::vector<double> result(d.begin(), d.end());
    std::sort(result.begin(), result.end());
    
    return result;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    std::vector<double> diagonal = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> off_diagonal = {-1.0, -1.0, -1.0, -1.0}; // size n-1
    
    // Compute eigenvalues using QR iteration
    std::vector<double> eigenvalues = qriter(diagonal, off_diagonal);
    
    // Sort the eigenvalues (though our algorithm should return them sorted)
    std::sort(eigenvalues.begin(), eigenvalues.end());
    
    // Generate JSON output
    std::ostringstream json;
    json << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < static_cast<int>(eigenvalues.size()); ++i) {
        if (i > 0) json << ",";
        json << std::fixed << std::setprecision(10) << eigenvalues[i];
    }
    json << "]}";
    
    std::cout << json.str() << std::endl;
    
    return 0;
}