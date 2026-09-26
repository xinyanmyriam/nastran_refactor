#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Portable M_PI definition (not standard C++)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Forward declarations for Eigen compatibility
using Eigen::VectorXd;
using Eigen::MatrixXd;

// Mock common block equivalents as global variables (simplified for test case)
int N = 5;                    // order of matrix
int MAX = 0;                   // max iterations (computed)
int NV = 0;                    // number of eigenvectors (not used in test)
int NE = 0;                    // number of eigenvalues to print (not used)
int NFOUND = 0;                // number of rigid modes (0 in test)
int IOPTN = 0;                 // option flag (not MGIV in test)
double DLMDAS = 0.0;           // shift parameter (not used in test)
int NEVER = 0;                 // failure counter

// Constants
const double EPSI = 1.0e-10;
const double ZERO = 0.0;
const double ONE = 1.0;

// QRITER implementation adapted from Fortran 77 - FIXED VERSION
void QRITER(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr) {
    // Initialize MAX = 100*N
    MAX = 100 * N;
    if (NV > N) NV = N;
    if (NE == 0) NE = N;
    if (NE < NV) NE = NV;

    NEVER = 0;
    
    // If QR != 0, just reorder (but in our test QR=0, so do full iteration)
    if (qr != 0) {
        // Reordering only - skip to sorting
        goto reorder;
    }

    // Start with full matrix [M1,M2] = [1,N]
    int M1 = 1;
    int M2 = N;
    
    // Main iteration loop
    while (M1 <= M2) {
        // Find decoupled submatrix: search from bottom up for zero off-diagonal
        // Look for the largest k such that o[k-1] == 0 (k from M1 to M2-1)
        int new_M2 = M2;
        for (int k = M2 - 1; k >= M1; --k) {
            if (std::abs(o[k-1]) > EPSI) {
                new_M2 = k + 1;
                break;
            }
        }
        
        // If no non-zero off-diagonal found, then M2 is isolated
        if (new_M2 == M1) {
            // Single element
            M1 = M2;
            M2 = M2 - 1;
            continue;
        }
        
        // Now find top of decoupled block by searching upward from new_M2-1
        int new_M1 = M1;
        for (int k = new_M2 - 1; k > M1; --k) {
            if (std::abs(o[k-2]) > EPSI) {
                new_M1 = k;
                break;
            }
        }
        
        M1 = new_M1;
        M2 = new_M2;
        
        // QR iteration for submatrix [M1, M2]
        bool converged = false;
        
        for (int iter = 1; iter <= MAX && !converged; ++iter) {
            // Check convergence: if bottom subdiagonal is small enough
            if (M1 < M2 && std::abs(o[M2-2]) <= EPSI * (std::abs(val[M2-2]) + std::abs(val[M2-1]))) {
                converged = true;
                continue;
            }
            
            // Compute Wilkinson shift for 2x2 bottom submatrix
            double SHIFT = ZERO;
            if (M2 - M1 > 1) {
                // 2x2 submatrix: [val[M2-2], o[M2-2]; o[M2-2], val[M2-1]]
                double a = val[M2-2];
                double b = o[M2-2];
                double c = val[M2-1];
                double trace = a + c;
                double det = a * c - b * b;
                double disc = trace * trace - 4.0 * det;
                
                if (disc >= 0.0) {
                    double sqrt_disc = std::sqrt(disc);
                    double lambda1 = 0.5 * (trace + sqrt_disc);
                    double lambda2 = 0.5 * (trace - sqrt_disc);
                    // Choose shift closest to val[M2-1]
                    if (std::abs(lambda1 - val[M2-1]) < std::abs(lambda2 - val[M2-1])) {
                        SHIFT = lambda1;
                    } else {
                        SHIFT = lambda2;
                    }
                } else {
                    // Complex eigenvalues, use trace/2
                    SHIFT = 0.5 * trace;
                }
            } else {
                // 1x1 case, no shift needed
                SHIFT = ZERO;
            }
            
            // Apply shift to diagonal elements
            for (int i = M1; i <= M2; ++i) {
                val[i-1] -= SHIFT;
            }
            
            // Perform QR step with implicit shift (bulge chasing)
            // Start with first element
            double d = val[M1-1];
            double e = o[M1-1];
            double r = std::sqrt(d*d + e*e);
            
            if (std::abs(r) > EPSI) {
                double c = d / r;
                double s = e / r;
                
                // Apply Givens rotation to first 2x2 block
                val[M1-1] = r;
                o[M1-1] = s * val[M1];
                val[M1] = c * val[M1] + s * o[M1-1];
                
                // Chase the bulge down the matrix
                for (int i = M1; i < M2 - 1; ++i) {
                    d = o[i-1];
                    e = val[i];
                    r = std::sqrt(d*d + e*e);
                    
                    if (std::abs(r) > EPSI) {
                        c = d / r;
                        s = e / r;
                        
                        o[i-1] = r;
                        val[i] = c * val[i+1] + s * o[i];
                        o[i] = s * val[i+1] - c * o[i];
                    } else {
                        o[i-1] = ZERO;
                        val[i] = val[i+1];
                        o[i] = ZERO;
                    }
                }
                
                // Last element
                if (M2 > M1 + 1) {
                    d = o[M2-2];
                    e = val[M2-1];
                    r = std::sqrt(d*d + e*e);
                    
                    if (std::abs(r) > EPSI) {
                        c = d / r;
                        s = e / r;
                        o[M2-2] = r;
                        val[M2-1] = c * val[M2-1];
                    } else {
                        o[M2-2] = ZERO;
                        val[M2-1] = ZERO;
                    }
                }
            }
            
            // Shift back
            if (SHIFT != ZERO) {
                for (int i = M1; i <= M2; ++i) {
                    val[i-1] += SHIFT;
                }
            }
            
            // Check for convergence again after QR step
            if (M1 < M2 && std::abs(o[M2-2]) <= EPSI * (std::abs(val[M2-2]) + std::abs(val[M2-1]))) {
                converged = true;
            }
        }
        
        if (!converged) {
            NEVER++;
        }
        
        // Deflate the matrix: if bottom subdiagonal is zero, we have an eigenvalue
        if (M1 < M2 && std::abs(o[M2-2]) <= EPSI) {
            // Found eigenvalue at val[M2-1]
            M2 = M2 - 1;
        } else if (M1 == M2) {
            // Single element done
            M1 = M1 + 1;
            M2 = M2 - 1;
        } else {
            // Continue with same bounds
        }
    }

reorder:
    // Reorder eigenvalues algebraically in ascending order
    // First, create temporary arrays for sorting
    std::vector<double> sorted_vals = val;
    std::vector<int> temp_loc(N);
    
    // Initialize loc to identity
    for (int i = 0; i < N; ++i) {
        loc[i] = i + 1; // 1-indexed positions
    }
    
    // Simple selection sort to match Fortran logic
    for (int k = 0; k < N; ++k) {
        int min_idx = k;
        for (int i = k; i < N; ++i) {
            if (sorted_vals[i] < sorted_vals[min_idx]) {
                min_idx = i;
            }
        }
        
        // Swap values
        double temp_val = sorted_vals[k];
        sorted_vals[k] = sorted_vals[min_idx];
        sorted_vals[min_idx] = temp_val;
        
        // Swap locations
        int temp_loc_val = loc[k];
        loc[k] = loc[min_idx];
        loc[min_idx] = temp_loc_val;
    }
    
    // Copy sorted values back
    val = sorted_vals;
    
    // If rigid modes were found, set first NFOUND to zero (not applicable here)
    if (NFOUND > 0) {
        for (int i = 0; i < NFOUND; ++i) {
            val[i] = ZERO;
        }
    }
}

// Analytical solution for test case: lambda_k = 2 - 2*cos(k*pi/6), k=1..5
std::vector<double> analytical_eigenvalues() {
    std::vector<double> exact(5);
    for (int k = 1; k <= 5; ++k) {
        exact[k-1] = 2.0 - 2.0 * std::cos(k * M_PI / 6.0);
    }
    return exact;
}

int main() {
    // Test case: 5x5 tridiagonal matrix
    // diagonal = [2,2,2,2,2], off-diagonal = [-1,-1,-1,-1]
    // Note: Fortran uses 1-based indexing, C++ uses 0-based
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0}; // diagonal
    std::vector<double> o = {-1.0, -1.0, -1.0, -1.0};   // off-diagonal (subdiagonal)
    std::vector<int> loc(5);                             // location array
    
    // Set global parameters for test
    N = 5;
    NV = 0;
    NE = 0;
    NFOUND = 0;
    IOPTN = 0;
    DLMDAS = 0.0;
    
    // Call QRITER with QR=0 (full iteration)
    QRITER(val, o, loc, 0);
    
    // Output as JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setprecision(15) << val[i];
        if (i < 4) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}