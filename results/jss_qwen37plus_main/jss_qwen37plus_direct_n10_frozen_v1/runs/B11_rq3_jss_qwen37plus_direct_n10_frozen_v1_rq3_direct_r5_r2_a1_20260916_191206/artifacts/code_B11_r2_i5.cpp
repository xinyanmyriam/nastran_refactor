#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Forward declaration for the main QR iteration function
void qriter(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr);

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    const int n = 5;
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> o = {-1.0, -1.0, -1.0, -1.0}; // ACTUAL off-diagonals, not squares
    std::vector<int> loc(n);
    
    // Initialize loc to 0-based indices
    for (int i = 0; i < n; ++i) {
        loc[i] = i;
    }
    
    // Call QR iteration (qr = 0 means perform iteration)
    qriter(val, o, loc, 0);
    
    // Sort eigenvalues in ascending order
    std::sort(val.begin(), val.end());
    
    // Output JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < n; ++i) {
        if (i > 0) std::cout << ",";
        std::cout << std::fixed << std::setprecision(12) << val[i];
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}

void qriter(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr) {
    const int n = static_cast<int>(val.size());
    if (n == 0) return;
    
    // Constants
    const double EPSI = 1.0e-12;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    
    // If qr != 0, just reorder (but we'll do full iteration for test case)
    if (qr != 0) {
        // Just sort the eigenvalues
        std::vector<double> temp = val;
        std::vector<int> temp_loc = loc;
        
        // Simple selection sort
        for (int k = 0; k < n; ++k) {
            int m = -1;
            for (int i = 0; i < n; ++i) {
                if (temp[i] != -10000.0) {
                    m = i;
                    break;
                }
            }
            if (m == -1) break;
            
            for (int i = m + 1; i < n; ++i) {
                if (temp[i] != -10000.0 && temp[m] > temp[i]) {
                    m = i;
                }
            }
            
            o[k] = temp[m];
            temp[m] = -10000.0;
            loc[k] = m;
        }
        
        for (int i = 0; i < n; ++i) {
            val[i] = o[i];
        }
        return;
    }
    
    // Main QR iteration loop - use standard implicit QR with Wilkinson shift
    const int max_iter = 30 * n;
    std::vector<double> d = val;  // diagonal
    std::vector<double> e = o;    // off-diagonal (subdiagonal)
    
    // Work on submatrices
    int l = 0;
    while (l < n) {
        // Find start of current submatrix
        int m = l;
        while (m < n - 1) {
            if (std::abs(e[m]) <= EPSI * (std::abs(d[m]) + std::abs(d[m+1]))) {
                e[m] = 0.0;
                break;
            }
            m++;
        }
        
        if (m == l) {
            // Single element
            l++;
            continue;
        }
        
        // QR iteration for submatrix from l to m (inclusive)
        int iter = 0;
        bool converged = false;
        
        while (!converged && iter < max_iter) {
            iter++;
            
            // Check convergence: if e[m-1] is small enough
            if (m > l && std::abs(e[m-1]) <= EPSI * (std::abs(d[m-1]) + std::abs(d[m]))) {
                e[m-1] = 0.0;
                converged = true;
                continue;
            }
            
            // Wilkinson shift: use eigenvalue of bottom 2x2 submatrix
            double s = 0.0;
            if (m > l) {
                double a = d[m-1];
                double b = e[m-1];
                double c = d[m];
                double disc = (a - c) * (a - c) + 4.0 * b * b;
                if (disc >= 0.0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (a + c + sqrt_disc);
                    double lambda2 = 0.5 * (a + c - sqrt_disc);
                    s = (std::abs(lambda1 - d[m]) < std::abs(lambda2 - d[m])) ? lambda1 : lambda2;
                } else {
                    s = d[m];
                }
            } else {
                s = d[l];
            }
            
            // Apply shift to submatrix: form matrix H - sI
            // Start the bulge chase from position l
            double f = d[l] - s;
            double g = e[l];
            
            // Perform QR step (Givens rotations) - corrected implementation
            for (int i = l; i < m; i++) {
                double r = std::sqrt(f * f + g * g);
                if (r == 0.0) {
                    f = 0.0;
                    g = 0.0;
                    break;
                }
                
                double c = f / r;
                double s_rot = g / r;
                
                // Update d[i] and e[i-1] (if i > l)
                if (i > l) {
                    e[i-1] = r;
                }
                
                // Update d[i] using rotation
                double temp_d = c * d[i] + s_rot * g;
                double temp_g = -s_rot * d[i] + c * g;
                
                // Update e[i] and prepare g for next iteration
                if (i < m) {
                    double temp_e = c * e[i] + s_rot * temp_g;
                    g = -s_rot * e[i] + c * temp_g;
                    e[i] = temp_e;
                } else {
                    g = temp_g;
                }
                
                d[i] = temp_d;
            }
            
            // Update the last diagonal element d[m]
            if (m > l) {
                e[m-1] = g;
                d[m] = d[m] + f; // This was the critical error - need to add the shift effect properly
            }
        }
        
        if (!converged) {
            // Failed to converge, use last value
        }
        
        l = m + 1;
    }
    
    // Copy results back
    val = d;
    o = e;
    
    // Reorder eigenvalues algebraically in ascending order
    std::vector<double> temp_val = val;
    std::vector<int> temp_loc = loc;
    
    for (int k = 0; k < n; ++k) {
        int m = -1;
        for (int i = 0; i < n; ++i) {
            if (temp_val[i] != -10000.0) {
                m = i;
                break;
            }
        }
        if (m == -1) break;
        
        for (int i = m + 1; i < n; ++i) {
            if (temp_val[i] != -10000.0 && temp_val[m] > temp_val[i]) {
                m = i;
            }
        }
        
        o[k] = temp_val[m];
        temp_val[m] = -10000.0;
        loc[k] = m;
    }
    
    // Copy back
    for (int i = 0; i < n; ++i) {
        val[i] = o[i];
    }
}