#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <string>
#include <sstream>

// For portability: use std::abs for double
using std::abs;

// INVERD implementation in C++17 using Eigen
// This implements the same algorithm as the Fortran INVERD subroutine:
// - Full pivoting (row and column exchanges)
// - In-place matrix inversion
// - Solution of linear system Kx = b
// - Determinant computation
// - Singularity detection

void INVERD(
    Eigen::MatrixXd& A,           // Input/output: matrix to invert (modified in-place)
    int N,                        // Size of the square submatrix to invert
    Eigen::MatrixXd& B,           // Input/output: RHS vectors (M columns), solution vectors on output
    int M,                        // Number of RHS vectors (0 if none)
    double& DETERM,               // Output: determinant
    int& ISING,                   // Output: 1=non-singular, 2=singular
    Eigen::MatrixXi& INDEX        // Workspace: N x 3 matrix
) {
    // Initialize
    const double EPSI = 1.0e-36;
    
    // On Windows, use larger epsilon for double precision stability
    // Mimicking the Fortran logic: if MACH == 5 (VAX), use 1e-18
    // Since we don't have MACH, use a reasonable value for modern systems
    const double eps = 1.0e-18;
    
    DETERM = 1.0;
    if (ISING < 0) {
        DETERM = 0.0;
    }
    
    // Initialize INDEX(:,3) = 0
    INDEX.col(2).setZero();
    
    // Main loop over pivot steps
    for (int I = 0; I < N; ++I) {
        // Search for pivot: find element with maximum absolute value among unmarked rows/columns
        double AMAX = 0.0;
        int IROW = -1;
        int ICOLUM = -1;
        
        for (int J = 0; J < N; ++J) {
            if (INDEX(J, 2) == 1) continue;  // Skip marked columns
            
            for (int K = 0; K < N; ++K) {
                if (INDEX(K, 2) == 1) continue;  // Skip marked rows
                
                double abs_val = abs(A(J, K));
                if (abs_val > AMAX) {
                    AMAX = abs_val;
                    IROW = J;
                    ICOLUM = K;
                }
            }
        }
        
        if (IROW == -1 || ICOLUM == -1) {
            ISING = 2;
            return;
        }
        
        // Mark this column as used
        INDEX(ICOLUM, 2) += 1;
        // Store row and column indices for later unscrambling
        INDEX(I, 0) = IROW;   // INDEX(I,1) in Fortran (1-indexed)
        INDEX(I, 1) = ICOLUM; // INDEX(I,2) in Fortran (1-indexed)
        
        // Interchange rows to put pivot on diagonal
        if (IROW != ICOLUM) {
            DETERM = -DETERM;
            
            // Swap rows IROW and ICOLUM in A
            A.row(IROW).swap(A.row(ICOLUM));
            
            // Swap rows in B if present
            if (M > 0) {
                for (int L = 0; L < M; ++L) {
                    double swap_val = B(IROW, L);
                    B(IROW, L) = B(ICOLUM, L);
                    B(ICOLUM, L) = swap_val;
                }
            }
        }
        
        // Divide pivot row by pivot element
        double PIVOT = A(ICOLUM, ICOLUM);
        DETERM *= PIVOT;
        
        if (abs(PIVOT) < eps) {
            ISING = 2;
            return;
        }
        
        // Normalize pivot row
        A.row(ICOLUM) /= PIVOT;
        
        // Normalize corresponding RHS entries
        if (M > 0) {
            B.row(ICOLUM) /= PIVOT;
        }
        
        // Reduce other rows
        for (int L1 = 0; L1 < N; ++L1) {
            if (L1 == ICOLUM) continue;
            
            double T = A(L1, ICOLUM);
            A(L1, ICOLUM) = 0.0;
            
            if (abs(T) < eps) continue;
            
            // Update row L1: A[L1,:] -= T * A[ICOLUM,:]
            A.row(L1) -= T * A.row(ICOLUM);
            
            // Update RHS if present
            if (M > 0) {
                B.row(L1) -= T * B.row(ICOLUM);
            }
        }
    }
    
    // Interchange columns (unscramble)
    for (int I = 0; I < N; ++I) {
        int L = N - 1 - I;  // Fortran: N+1-I (1-indexed)
        if (INDEX(L, 0) == INDEX(L, 1)) continue;
        
        int JROW = INDEX(L, 0);
        int JCOLUM = INDEX(L, 1);
        
        // Swap columns JROW and JCOLUM in A
        A.col(JROW).swap(A.col(JCOLUM));
    }
    
    // Check singularity: all columns must be marked
    ISING = 1;
    for (int K = 0; K < N; ++K) {
        if (INDEX(K, 2) != 1) {
            ISING = 2;
            break;
        }
    }
}

int main() {
    // Test case: 4x4 SPD matrix K
    // K = [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
    Eigen::MatrixXd K(4, 4);
    K << 200.0, -100.0,   0.0,   0.0,
        -100.0,  200.0, -100.0,   0.0,
           0.0, -100.0,  200.0, -100.0,
           0.0,    0.0, -100.0,  100.0;
    
    // Compute b = K * [1,2,3,4]^T
    Eigen::VectorXd x_true(4);
    x_true << 1.0, 2.0, 3.0, 4.0;
    Eigen::VectorXd b = K * x_true;
    
    // Prepare inputs for INVERD
    Eigen::MatrixXd A = K;  // Copy for modification
    Eigen::MatrixXd B(4, 1);
    B.col(0) = b;
    
    int N = 4;
    int M = 1;
    double DETERM = 0.0;
    int ISING = 1;
    Eigen::MatrixXi INDEX(N, 3);  // N x 3 workspace
    
    // Call INVERD
    INVERD(A, N, B, M, DETERM, ISING, INDEX);
    
    // Extract solution
    Eigen::VectorXd solution = B.col(0);
    
    // Compute inverse: solve A * X = I
    Eigen::MatrixXd A_for_inverse = K;
    Eigen::MatrixXd B_inv(4, 4);
    B_inv.setIdentity();
    Eigen::MatrixXi INDEX_inv(4, 3);
    double DETERM_inv = 0.0;
    int ISING_inv = 1;
    
    INVERD(A_for_inverse, 4, B_inv, 4, DETERM_inv, ISING_inv, INDEX_inv);
    
    // Format output as JSON-like plain text
    std::cout << std::fixed << std::setprecision(6);

    // test
    std::cout << "\"test\": \"INVERD\"\n";

    // solution
    std::cout << "\"solution\": [";
    for (int i = 0; i < 4; ++i) {
        std::cout << solution(i);
        if (i < 3) std::cout << ", ";
    }
    std::cout << "]\n";

    // determinant
    std::cout << "\"determinant\": " << DETERM << "\n";

    // inverse matrix
    std::cout << "\"inverse\": [\n";
    for (int i = 0; i < 4; ++i) {
        std::cout << " [";
        for (int j = 0; j < 4; ++j) {
            std::cout << B_inv(i, j);
            if (j < 3) std::cout << ", ";
        }
        std::cout << "]";
        if (i < 3) std::cout << ",";
        std::cout << "\n";
    }
    std::cout << "]\n";

    return 0;
}