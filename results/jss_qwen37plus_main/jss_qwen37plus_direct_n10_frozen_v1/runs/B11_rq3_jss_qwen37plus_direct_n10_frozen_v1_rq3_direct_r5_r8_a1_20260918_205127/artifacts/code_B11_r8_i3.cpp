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
    
    const int max_iter = 30 * n; // More robust iteration limit
    
    for (int iter = 0; iter < max_iter && m > 0; ++iter) {
        // Look for deflation from bottom
        int l = m;
        while (l > 0 && std::abs(e[l-1]) > EPSI * (std::abs(d[l-1]) + std::abs(d[l]))) {
            l--;
        }
        
        if (l == m) {
            // Isolated eigenvalue at position m
            m--;
            continue;
        }
        
        if (l == 0) {
            // Entire matrix is unreduced, use Wilkinson shift from bottom 2x2
            double a = d[m-1];
            double b = e[m-1];
            double c = d[m];
            double trace = a + c;
            double det = a * c - b * b;
            double disc = trace * trace - 4.0 * det;
            double sqrt_disc = std::sqrt(std::max(0.0, disc));
            double lambda1 = 0.5 * (trace + sqrt_disc);
            double lambda2 = 0.5 * (trace - sqrt_disc);
            double s = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
            
            // Implicit QR step with shift s
            double x = d[0] - s;
            double y = e[0];
            
            // Chase the bulge from top to bottom
            for (int k = 0; k < m; ++k) {
                // Compute Givens rotation to zero y
                double r = std::sqrt(x*x + y*y);
                if (r == ZERO) {
                    break;
                }
                double c_rot = x / r;
                double s_rot = y / r;
                
                // Apply to rows k and k+1 (affecting d[k], e[k], d[k+1])
                // First, update d[k] and e[k]
                double temp_d = c_rot * d[k] + s_rot * e[k];
                double temp_e = -s_rot * d[k] + c_rot * e[k];
                d[k] = temp_d;
                e[k] = temp_e;
                
                // Then update e[k] and d[k+1] (but e[k] was just updated, so use original for next)
                // Actually, for the next step, we need the new e[k] and d[k+1]
                if (k + 1 < n) {
                    double temp_d_next = c_rot * e[k] + s_rot * d[k+1];
                    double temp_e_next = -s_rot * e[k] + c_rot * d[k+1];
                    if (k < n-1) e[k] = temp_d_next;
                    if (k + 1 < n) d[k+1] = temp_e_next;
                }
                
                // Prepare for next rotation: x becomes e[k] (after update), y becomes e[k+1] if exists
                if (k + 1 < m) {
                    x = e[k];
                    y = e[k+1];
                } else if (k == m-1) {
                    x = e[k];
                    y = ZERO;
                }
            }
        } else {
            // Use Wilkinson shift from bottom 2x2 submatrix
            double a = d[m-1];
            double b = e[m-1];
            double c = d[m];
            double trace = a + c;
            double det = a * c - b * b;
            double disc = trace * trace - 4.0 * det;
            double sqrt_disc = std::sqrt(std::max(0.0, disc));
            double lambda1 = 0.5 * (trace + sqrt_disc);
            double lambda2 = 0.5 * (trace - sqrt_disc);
            double s = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
            
            // Implicit QR step with shift s on submatrix [l:m]
            // Start bulge at position l
            double x = d[l] - s;
            double y = (l < n-1) ? e[l] : ZERO;
            
            // Chase bulge from l to m-1
            for (int k = l; k < m; ++k) {
                // Compute Givens rotation to zero the subdiagonal element
                double r = std::sqrt(x*x + y*y);
                if (r == ZERO) {
                    break;
                }
                double c_rot = x / r;
                double s_rot = y / r;
                
                // Apply rotation to columns k and k+1
                // Update d[k] and e[k]
                double temp_d = c_rot * d[k] + s_rot * e[k];
                double temp_e = -s_rot * d[k] + c_rot * e[k];
                d[k] = temp_d;
                if (k < n-1) e[k] = temp_e;
                
                // Update e[k] and d[k+1]
                if (k + 1 < n) {
                    double temp_d_next = c_rot * e[k] + s_rot * d[k+1];
                    double temp_e_next = -s_rot * e[k] + c_rot * d[k+1];
                    if (k < n-1) e[k] = temp_d_next;
                    if (k + 1 < n) d[k+1] = temp_e_next;
                }
                
                // Prepare for next rotation
                if (k + 1 < m) {
                    x = e[k];
                    y = e[k+1];
                } else if (k == m-1) {
                    x = e[k];
                    y = ZERO;
                }
            }
        }
        
        // Check for deflation at bottom
        if (m > 0 && std::abs(e[m-1]) <= EPSI * (std::abs(d[m-1]) + std::abs(d[m]))) {
            m--;
        }
    }
    
    // Return eigenvalues (they should be in d[0..n-1], but sort to be safe)
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