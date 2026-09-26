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

// Portable M_PI definition
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// QR iteration for tridiagonal matrix eigenvalues
// Input: diagonal vector 'd' and off-diagonal vector 'e' (size n-1)
// Output: eigenvalues in 'd' (sorted ascending)
void QRITER(VectorXd& d, const VectorXd& e) {
    int n = static_cast<int>(d.size());
    if (n <= 0) return;
    
    // Handle trivial case
    if (n == 1) return;
    
    // Copy off-diagonal squares: o[i] = e[i]^2, for i in [0, n-2]
    VectorXd o = e.array().square();
    
    // Simulate common blocks
    int MAX = 100 * n;
    
    // Work on submatrices [l, u], 0-indexed
    int l = 0;
    int u = n - 1;
    
    while (l < u) {
        // Deflation: find largest m such that [l, m] is decoupled
        // i.e., find rightmost nonzero off-diagonal in [l, u-1]
        int m = u;
        for (int k = u - 1; k >= l; --k) {
            // o[k] corresponds to off-diagonal between k and k+1
            if (k < o.size() && std::abs(o(k)) > EPSI) {
                m = k + 1; // submatrix goes to index k+1 (so [l, k+1])
                break;
            }
        }
        
        // If no nonzero off-diagonal found in [l, u-1], then [l, u] is decoupled
        if (m == u && (u - 1 < 0 || u - 1 >= o.size() || std::abs(o(u-1)) <= EPSI)) {
            // All off-diagonals zero => eigenvalues are diagonals
            // Just proceed to next block
            u = u - 1;
            continue;
        }
        
        // Now work on submatrix [l, m-1] (size = m - l), where m-1 >= l
        // But ensure m-1 <= u and m > l
        if (m <= l) {
            l = l + 1;
            continue;
        }
        u = m - 1;
        
        // QR iteration on [l, u]
        int iter_count = 0;
        bool converged = false;
        
        for (iter_count = 0; iter_count < MAX; ++iter_count) {
            // Check deflation at bottom: if |e[u-1]| is negligible relative to |d[u]| and |d[u-1]|
            if (u > l) {
                int idx = u - 1; // o[idx] is off-diagonal between u-1 and u
                if (idx >= 0 && idx < o.size() && std::abs(o(idx)) <= EPSI * (std::abs(d(u)) + std::abs(d(u-1)))) {
                    // Deflate: treat u as separate eigenvalue
                    u = u - 1;
                    if (u <= l) {
                        converged = true;
                        break;
                    }
                    continue;
                }
            }
            
            // Wilkinson shift: use eigenvalue of bottom 2x2 submatrix
            double shift = ZERO;
            if (u > l) {
                double a = d(u-1);
                double b = d(u);
                double c_val = o(u-1); // off-diagonal between u-1 and u
                double trace = a + b;
                double det = a * b - c_val;
                double disc = trace * trace - 4.0 * det;
                if (disc >= 0.0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (trace - sqrt_disc);
                    double lambda2 = 0.5 * (trace + sqrt_disc);
                    // Choose shift closest to d[u]
                    shift = (std::abs(lambda1 - d(u)) < std::abs(lambda2 - d(u))) ? lambda1 : lambda2;
                } else {
                    shift = 0.5 * trace;
                }
            } else {
                shift = d(l);
            }
            
            // Apply shift
            for (int i = l; i <= u; ++i) {
                d(i) -= shift;
            }
            
            // Perform implicit QR step (chasing bulge)
            double f = d(l);
            double g = o(l);
            double c = 0.0, s = 0.0; // Declare c and s here to be accessible later
            
            for (int k = l; k < u; ++k) {
                double r = std::sqrt(f*f + g*g);
                if (r == ZERO) {
                    r = EPSI;
                }
                c = f / r;
                s = g / r;
                
                // Update d[k]
                if (k > l) {
                    o(k-1) = s * f + c * g;
                }
                f = c * d(k) + s * o(k);
                d(k) = c * d(k) - s * o(k);
                
                if (k < u) {
                    g = s * d(k+1);
                    d(k+1) = c * d(k+1) - s * o(k);
                    o(k) = c * o(k) + s * d(k+1);
                }
            }
            
            // Update last element - now c and s are in scope
            o(u-1) = f * s;
            d(u) = f * c;
            
            // Undo shift
            for (int i = l; i <= u; ++i) {
                d(i) += shift;
            }
            
            // Check convergence: if any off-diagonal is negligible, deflate
            bool deflated = false;
            for (int k = l; k < u; ++k) {
                if (k < o.size() && std::abs(o(k)) <= EPSI * (std::abs(d(k)) + std::abs(d(k+1)))) {
                    o(k) = ZERO;
                    // Split here
                    u = k;
                    deflated = true;
                    break;
                }
            }
            if (deflated) {
                break;
            }
        }
        
        // Move to next block
        if (u <= l) {
            l = l + 1;
            u = n - 1;
        } else {
            // Continue with [l, u-1] if deflated at bottom, else adjust
            if (u > l) {
                u = u - 1;
            }
        }
    }
    
    // Final sort ascending
    std::vector<double> temp;
    for (int i = 0; i < n; ++i) {
        temp.push_back(d(i));
    }
    std::sort(temp.begin(), temp.end());
    for (int i = 0; i < n; ++i) {
        d(i) = temp[i];
    }
}

// Analytical solution for test case: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
std::vector<double> analytical_eigenvalues() {
    std::vector<double> result;
    for (int k = 1; k <= 5; ++k) {
        double lambda = 2.0 - 2.0 * std::cos(k * M_PI / 6.0);
        result.push_back(lambda);
    }
    return result;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    int n = 5;
    VectorXd d(n);
    VectorXd e(n-1);
    
    // Initialize diagonal
    for (int i = 0; i < n; ++i) {
        d(i) = 2.0;
    }
    
    // Initialize off-diagonal (subdiagonal, symmetric)
    for (int i = 0; i < n-1; ++i) {
        e(i) = -1.0;
    }
    
    // Call QRITER
    QRITER(d, e);
    
    // Extract and sort eigenvalues (QRITER does sort, but double-check)
    std::vector<double> eigenvalues;
    for (int i = 0; i < n; ++i) {
        eigenvalues.push_back(d(i));
    }
    std::sort(eigenvalues.begin(), eigenvalues.end());
    
    // Format JSON output
    std::ostringstream json;
    json << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (i > 0) json << ",";
        json << std::fixed << std::setprecision(12) << eigenvalues[i];
    }
    json << "]}";
    
    std::cout << json.str() << std::endl;
    
    return 0;
}