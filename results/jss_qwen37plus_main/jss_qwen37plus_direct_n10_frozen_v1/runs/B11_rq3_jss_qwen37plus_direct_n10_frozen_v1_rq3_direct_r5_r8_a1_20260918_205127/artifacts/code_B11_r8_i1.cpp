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
    int m = n;
    int max_iter = 30; // Standard max iterations per eigenvalue
    
    while (m > 1) {
        // Find the largest k such that e[k] is negligible
        // Search from bottom up for first non-negligible e[i]
        int k;
        for (k = m - 2; k >= 0; --k) {
            if (std::abs(e[k]) > EPSI * (std::abs(d[k]) + std::abs(d[k+1]))) {
                break;
            }
        }
        
        if (k < 0) {
            // All off-diagonals are zero, we're done
            break;
        }
        
        // If k == m-2, then the last 2x2 block is coupled
        // Otherwise, we have a submatrix from 0 to k+1
        if (k == m - 2) {
            // Apply QR step to the bottom 2x2 submatrix
            // Compute Wilkinson shift from bottom 2x2 submatrix
            double a = d[m-2], b = e[m-2], c = d[m-1];
            double mu;
            if (m == 2) {
                // For 2x2, use exact eigenvalues
                double disc = (a - c) * (a - c) + 4.0 * b * b;
                if (disc >= 0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (a + c - sqrt_disc);
                    double lambda2 = 0.5 * (a + c + sqrt_disc);
                    // Use smaller eigenvalue as shift
                    mu = (lambda1 < lambda2) ? lambda1 : lambda2;
                } else {
                    mu = a;
                }
            } else {
                // Wilkinson shift: eigenvalue of bottom 2x2 closest to d[m-1]
                double trace = a + c;
                double det = a * c - b * b;
                double disc = trace * trace - 4.0 * det;
                if (disc >= 0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (trace - sqrt_disc);
                    double lambda2 = 0.5 * (trace + sqrt_disc);
                    // Choose shift closest to d[m-1]
                    mu = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
                } else {
                    mu = c;
                }
            }
            
            // Apply shift
            d[m-2] -= mu;
            d[m-1] -= mu;
            
            // Perform one QR step on 2x2
            double f = d[m-2];
            double g = e[m-2];
            double h = d[m-1];
            
            double r = std::sqrt(f * f + g * g);
            double c_val = f / r;
            double s = g / r;
            
            // Update d[m-2] and d[m-1]
            d[m-2] = c_val * f - s * g;
            d[m-1] = s * f + c_val * h;
            
            // Update e[m-2]
            e[m-2] = s * h;
            
            // Add shift back
            d[m-2] += mu;
            d[m-1] += mu;
            
            // Check convergence
            if (std::abs(e[m-2]) <= EPSI * (std::abs(d[m-2]) + std::abs(d[m-1]))) {
                e[m-2] = 0.0;
                m--;
            }
        } else {
            // We have a larger submatrix, need to apply implicit QR step
            // Use Wilkinson shift from bottom 2x2 of current submatrix [k+1, m-1]
            int i0 = k + 1;
            if (i0 < m - 1) {
                double a = d[m-2], b = e[m-2], c = d[m-1];
                double trace = a + c;
                double det = a * c - b * b;
                double disc = trace * trace - 4.0 * det;
                double mu;
                if (disc >= 0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (trace - sqrt_disc);
                    double lambda2 = 0.5 * (trace + sqrt_disc);
                    mu = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
                } else {
                    mu = c;
                }
                
                // Apply shift to entire submatrix [i0, m-1]
                for (int i = i0; i < m; ++i) {
                    d[i] -= mu;
                }
                
                // Perform implicit QR step using Givens rotations
                double f = d[i0];
                double g = e[i0];
                
                for (int i = i0; i < m - 1; ++i) {
                    double r = std::sqrt(f * f + g * g);
                    double c_val = f / r;
                    double s = g / r;
                    
                    // Update d[i] and d[i+1]
                    double temp = c_val * d[i] + s * e[i];
                    e[i] = c_val * e[i] - s * d[i+1];
                    d[i] = temp;
                    
                    if (i > i0) {
                        e[i-1] = c_val * e[i-1] + s * d[i];
                        d[i] = c_val * d[i] - s * e[i-1];
                    }
                    
                    // Update next f,g
                    if (i < m - 2) {
                        f = d[i+1];
                        g = e[i+1];
                    }
                }
                
                // Add shift back
                for (int i = i0; i < m; ++i) {
                    d[i] += mu;
                }
            }
        }
        
        // Check for deflation
        bool deflated = false;
        for (int i = 0; i < m - 1; ++i) {
            if (std::abs(e[i]) <= EPSI * (std::abs(d[i]) + std::abs(d[i+1]))) {
                e[i] = 0.0;
                deflated = true;
                // Split matrix at position i+1
                m = i + 1;
                break;
            }
        }
        
        if (!deflated) {
            // No deflation, continue with same m
        }
    }
    
    // Now handle remaining submatrices by simple iteration
    // Reset m to n and do standard QR iteration
    m = n;
    int iter_count = 0;
    const int max_total_iter = 1000;
    
    while (m > 1 && iter_count < max_total_iter) {
        iter_count++;
        
        // Find deflation point from bottom
        int k = m - 2;
        while (k >= 0 && std::abs(e[k]) <= EPSI * (std::abs(d[k]) + std::abs(d[k+1]))) {
            k--;
        }
        
        if (k < 0) {
            break; // All deflated
        }
        
        // Determine shift from bottom 2x2 submatrix
        double a = d[k], b = e[k], c = d[k+1];
        double trace = a + c;
        double det = a * c - b * b;
        double disc = trace * trace - 4.0 * det;
        double mu;
        if (disc >= 0) {
            double sqrt_disc = std::sqrt(disc);
            double lambda1 = 0.5 * (trace - sqrt_disc);
            double lambda2 = 0.5 * (trace + sqrt_disc);
            mu = (std::abs(lambda1 - c) < std::abs(lambda2 - c)) ? lambda1 : lambda2;
        } else {
            mu = c;
        }
        
        // Apply shift
        for (int i = 0; i < m; ++i) {
            d[i] -= mu;
        }
        
        // Perform QR step: apply Givens rotations from top to bottom
        double f = d[0];
        double g = e[0];
        
        for (int i = 0; i < m - 1; ++i) {
            double r = std::sqrt(f * f + g * g);
            if (r == 0.0) {
                r = EPSI;
            }
            double c_val = f / r;
            double s = g / r;
            
            // Update d[i] and d[i+1]
            double temp = c_val * d[i] + s * e[i];
            e[i] = c_val * e[i] - s * d[i+1];
            d[i] = temp;
            
            if (i > 0) {
                e[i-1] = c_val * e[i-1] + s * d[i];
                d[i] = c_val * d[i] - s * e[i-1];
            }
            
            // Update for next iteration
            if (i < m - 2) {
                f = d[i+1];
                g = e[i+1];
            }
        }
        
        // Add shift back
        for (int i = 0; i < m; ++i) {
            d[i] += mu;
        }
        
        // Check for deflation
        for (int i = 0; i < m - 1; ++i) {
            if (std::abs(e[i]) <= EPSI * (std::abs(d[i]) + std::abs(d[i+1]))) {
                e[i] = 0.0;
                m = i + 1;
                break;
            }
        }
    }
    
    // Sort eigenvalues in ascending order
    std::vector<double> result = d;
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