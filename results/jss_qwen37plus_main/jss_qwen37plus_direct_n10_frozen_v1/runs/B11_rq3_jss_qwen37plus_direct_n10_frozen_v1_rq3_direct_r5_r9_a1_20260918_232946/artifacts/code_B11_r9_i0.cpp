#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <Eigen/Dense>

// Simulate the Fortran COMMON blocks and global parameters
struct GlobalState {
    int N = 0;           // order of the problem
    int MAX = 0;         // maximum iterations
    int NV = 0;          // number of eigenvectors
    int NE = 0;          // number of eigenvalues to print
    int NFOUND = 0;      // number of rigid modes found
    int IOPTN = 0;       // option flag
    double LFREQ = 0.0;  // low frequency bound
    double HFREQ = 0.0;  // high frequency bound
    double DLMDAS = 0.0; // shift parameter for MGIV method
    int NEVER = 0;       // iteration failure counter
};

GlobalState g_state;

// QR iteration for tridiagonal matrix eigenvalues
// Input: diagonal (val), off-diagonal (o) - both as vectors of doubles
// Output: eigenvalues in val (sorted ascending)
void QRITER(std::vector<double>& val, std::vector<double>& o, std::vector<int>& loc, int qr) {
    const double EPSI = 1.0e-10;
    const double ZERO = 0.0;
    const double ONE = 1.0;
    
    // Initialize loc array
    for (int i = 0; i < static_cast<int>(val.size()); ++i) {
        loc[i] = i + 1; // Fortran 1-based indexing
    }
    
    // Set up global state
    g_state.N = static_cast<int>(val.size());
    g_state.MAX = 100 * g_state.N;
    if (g_state.NV > g_state.N) g_state.NV = g_state.N;
    if (g_state.NE == 0) g_state.NE = g_state.N;
    if (g_state.NE < g_state.NV) g_state.NE = g_state.NV;
    
    // If QR == 0, perform actual QR iteration
    if (qr == 0) {
        int M1 = 1;
        int M2 = g_state.N;
        
        // Main iteration loop
        bool converged = false;
        for (int iter = 1; iter <= g_state.MAX && !converged; ++iter) {
            // Check for decoupled submatrix
            int M2M1 = M2 - 1;
            bool all_zeros = true;
            for (int k = 1; k <= M2M1; ++k) {
                int idx = k - 1; // convert to 0-based
                if (std::abs(o[idx]) > EPSI) {
                    all_zeros = false;
                    break;
                }
            }
            
            if (all_zeros) {
                converged = true;
                continue;
            }
            
            // Find decoupled submatrix boundaries
            M2M1 = M2 - 1;
            int new_M2 = M2;
            for (int k = 1; k <= M2M1; ++k) {
                int M1_temp = M2 - k;
                if (std::abs(o[M1_temp - 1]) > EPSI) {
                    new_M2 = M1_temp + 1;
                    break;
                }
            }
            M2 = new_M2;
            M2M1 = M2 - 1;
            
            if (M2M1 == 1) {
                // Only two elements, handle separately
                int MM = 1;
                // Continue with iteration
            } else {
                // Find M1 for decoupled submatrix
                int MM = M1;
                for (int k = 2; k <= M2M1; ++k) {
                    int M1_temp = M2 - k;
                    if (std::abs(o[M1_temp - 1]) > EPSI) {
                        MM = M1_temp + 1;
                        break;
                    }
                }
                M1 = MM;
            }
            
            // Check convergence condition
            if (std::abs(val[M2-1]) + std::abs(o[M2M1-1]) == std::abs(val[M2-1])) {
                converged = true;
                continue;
            }
            
            // Check if all diagonal elements are equal
            bool all_equal = true;
            for (int k = M1; k <= M2M1; ++k) {
                if (std::abs(val[k-1] - val[k]) > EPSI) {
                    all_equal = false;
                    break;
                }
            }
            
            double SHIFT = ZERO;
            if (!all_equal) {
                // Find smallest diagonal term
                SHIFT = val[M2-1];
                for (int i = M1; i <= M2M1; ++i) {
                    if (std::abs(val[i-1]) < std::abs(SHIFT)) {
                        SHIFT = val[i-1];
                    }
                }
                
                // Reduce all terms by shift
                for (int i = M1; i <= M2; ++i) {
                    val[i-1] -= SHIFT;
                }
            }
            
            // QR iteration
            double R = val[M1-1] * val[M1-1];
            double S = o[M1-1] / (R + o[M1-1]);
            double T = ZERO;
            double U = S * (val[M1-1] + val[M1]);
            val[M1-1] += U;
            
            if (M1 != M2M1) {
                for (int i = M1+1; i <= M2M1; ++i) {
                    double G = val[i-1] - U;
                    R = (ONE - T) * o[i-2];
                    double ONES = ONE - S;
                    if (std::abs(ONES) > EPSI) {
                        R = G * G / ONES;
                    }
                    R = R + o[i-1];
                    o[i-2] = S * R;
                    if (std::abs(o[i-2]) <= EPSI) {
                        M1 = i;
                    }
                    T = S;
                    S = o[i-1] / R;
                    U = S * (G + val[i]);
                    val[i-1] = U + G;
                }
            }
            
            val[M2-1] -= U;
            R = (ONE - T) * o[M2M1-1];
            double ONES = ONE - S;
            if (std::abs(ONES) > EPSI) {
                R = val[M2-1] * val[M2-1] / ONES;
            }
            o[M2M1-1] = S * R;
            
            // Shift back
            if (std::abs(SHIFT) > EPSI) {
                for (int i = M1; i <= M2; ++i) {
                    val[i-1] += SHIFT;
                }
            }
            
            // Update M1 for next iteration
            int MM = M1;
            for (int i = M1; i <= M2M1; ++i) {
                if (std::abs(o[i-1]) <= EPSI) {
                    MM = i + 1;
                    break;
                }
            }
            M1 = MM;
            
            // Check convergence
            if (M1 == M2M1) {
                if (M1 == M2) {
                    converged = true;
                } else {
                    M2 = M2M1;
                    M2M1 = M2 - 1;
                    iter = 0; // restart iteration count for this submatrix
                }
            } else if (M1 <= 2) {
                converged = true;
            } else {
                M2 = M1 - 1;
                M1 = 1;
                iter = 0; // restart iteration count
            }
        }
        
        // Handle iteration failure
        if (!converged) {
            g_state.NEVER++;
        }
    }
    
    // Reorder eigenvalues algebraically in ascending order
    // First, create temporary arrays for sorting
    std::vector<double> temp_val = val;
    std::vector<int> temp_loc(loc.size());
    
    // Simple selection sort (mimics Fortran logic)
    for (int k = 0; k < g_state.N; ++k) {
        int min_idx = -1;
        for (int m = 0; m < g_state.N; ++m) {
            if (temp_val[m] != -10000.0) {
                min_idx = m;
                break;
            }
        }
        if (min_idx == -1) break;
        
        // Find minimum element
        for (int i = min_idx + 1; i < g_state.N; ++i) {
            if (temp_val[i] != -10000.0 && temp_val[i] < temp_val[min_idx]) {
                min_idx = i;
            }
        }
        
        temp_loc[k] = min_idx + 1; // 1-based indexing
        val[k] = temp_val[min_idx];
        temp_val[min_idx] = -10000.0;
    }
    
    // Copy back
    for (int i = 0; i < g_state.N; ++i) {
        loc[i] = temp_loc[i];
    }
    
    // Handle rigid modes
    if (g_state.NFOUND > 0) {
        for (int i = 0; i < g_state.NFOUND; ++i) {
            val[i] = ZERO;
        }
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
    std::vector<double> val = {2.0, 2.0, 2.0, 2.0, 2.0};
    std::vector<double> o = {-1.0, -1.0, -1.0, -1.0}; // off-diagonals (subdiagonal only)
    std::vector<int> loc(5);
    
    // Set global parameters for test
    g_state.N = 5;
    g_state.NV = 0;
    g_state.NE = 0;
    g_state.NFOUND = 0;
    g_state.IOPTN = 0;
    g_state.LFREQ = 0.0;
    g_state.HFREQ = 0.0;
    g_state.DLMDAS = 0.0;
    
    // Perform QR iteration
    QRITER(val, o, loc, 0);
    
    // Sort eigenvalues (the algorithm should produce sorted, but ensure)
    std::sort(val.begin(), val.end());
    
    // Output JSON
    std::cout << "{\"test\":\"QRITER\",\"eigenvalues\":[";
    for (int i = 0; i < 5; ++i) {
        std::cout << std::fixed << std::setprecision(12) << val[i];
        if (i < 4) std::cout << ",";
    }
    std::cout << "]}" << std::endl;
    
    return 0;
}