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

// QRITER implementation adapted from Fortran 77
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
        int M2M1 = M2 - 1;
        int new_M1 = M1;
        int new_M2 = M2;
        
        // Search upward from M2-1 to find first non-zero off-diagonal
        for (int k = 1; k <= M2M1; ++k) {
            int candidate = M2 - k;
            if (o[candidate-1] != ZERO) { // o is 0-indexed, Fortran is 1-indexed
                new_M2 = candidate + 1;
                break;
            }
        }
        
        // If all off-diagonals are zero, we're done for this block
        bool all_zero = true;
        for (int i = M1; i < M2; ++i) {
            if (o[i-1] != ZERO) {
                all_zero = false;
                break;
            }
        }
        if (all_zero) {
            M1 = M2;
            continue;
        }
        
        // Now find top of decoupled block
        int temp_M1 = new_M2;
        for (int k = 2; k <= new_M2 - 1; ++k) {
            int candidate = new_M2 - k;
            if (o[candidate-1] == ZERO) {
                temp_M1 = candidate + 1;
                break;
            }
        }
        M1 = temp_M1;
        M2 = new_M2;
        
        // QR iteration for submatrix [M1, M2]
        int MM = M1;
        bool converged = false;
        
        for (int iter = 1; iter <= MAX; ++iter) {
            // Check convergence: if |val[M2]| + o[M2-1] == |val[M2]| => o[M2-1] == 0
            if (std::abs(o[M2-1]) <= EPSI) {
                converged = true;
                break;
            }
            
            // Check if all diagonal elements equal
            bool all_equal = true;
            for (int k = M1; k < M2; ++k) {
                if (std::abs(val[k-1] - val[k]) > EPSI) {
                    all_equal = false;
                    break;
                }
            }
            
            double SHIFT = ZERO;
            if (!all_equal) {
                // Find smallest diagonal term (by absolute value)
                SHIFT = val[M2-1];
                for (int i = M1; i <= M2; ++i) {
                    if (std::abs(val[i-1]) < std::abs(SHIFT)) {
                        SHIFT = val[i-1];
                    }
                }
                
                // Shift the diagonal
                for (int i = M1; i <= M2; ++i) {
                    val[i-1] -= SHIFT;
                }
            }
            
            // QR iteration step
            double R = val[M1-1] * val[M1-1];
            double S = o[M1-1] / (R + o[M1-1]);
            double T = ZERO;
            double U = S * (val[M1-1] + val[M1]);
            val[M1-1] += U;
            
            if (M1 == M2M1) {
                // 2x2 case
                val[M2-1] -= U;
                R = (ONE - T) * o[M2M1-1];
                double ONES = ONE - S;
                if (std::abs(ONES) > EPSI) {
                    R = val[M2-1] * val[M2-1] / ONES;
                }
                o[M2M1-1] = S * R;
            } else {
                // General case
                int M1P1 = M1 + 1;
                for (int i = M1P1; i <= M2M1; ++i) {
                    double G = val[i-1] - U;
                    R = (ONE - T) * o[i-2];
                    double ONES = ONE - S;
                    if (std::abs(ONES) > EPSI) {
                        R = G * G / ONES;
                    }
                    R = R + o[i-1];
                    o[i-2] = S * R;
                    if (std::abs(o[i-2]) <= EPSI) {
                        MM = i;
                    }
                    T = S;
                    S = o[i-1] / R;
                    U = S * (G + val[i]);
                    val[i-1] = U + G;
                }
                
                val[M2-1] -= U;
                R = (ONE - T) * o[M2M1-1];
                double ONES = ONE - S;
                if (std::abs(ONES) > EPSI) {
                    R = val[M2-1] * val[M2-1] / ONES;
                }
                o[M2M1-1] = S * R;
            }
            
            // Shift back
            if (SHIFT != ZERO) {
                for (int i = M1; i <= M2; ++i) {
                    val[i-1] += SHIFT;
                }
            }
            
            M1 = MM;
        }
        
        if (!converged) {
            NEVER++;
        }
        
        // Update bounds for next submatrix
        if (M1 == M2M1) {
            M2 = M2M1;
            M2M1 = M2 - 1;
        } else {
            if (M1 <= 2) {
                break;
            }
            M2 = M1 - 1;
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
        int min_idx = -1;
        for (int m = 0; m < N; ++m) {
            if (sorted_vals[m] != -10000.0) {
                min_idx = m;
                break;
            }
        }
        if (min_idx == -1) break;
        
        // Find minimum element
        for (int i = min_idx + 1; i < N; ++i) {
            if (sorted_vals[i] != -10000.0 && sorted_vals[i] < sorted_vals[min_idx]) {
                min_idx = i;
            }
        }
        
        // Store in output
        val[k] = sorted_vals[min_idx];
        sorted_vals[min_idx] = -10000.0;
        loc[k] = min_idx + 1; // 1-indexed
    }
    
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
    
    // Sort the computed eigenvalues (QRITER should have done this, but ensure)
    std::sort(val.begin(), val.end());
    
    // Output as JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::setprecision(15) << val[i];
        if (i < 4) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}