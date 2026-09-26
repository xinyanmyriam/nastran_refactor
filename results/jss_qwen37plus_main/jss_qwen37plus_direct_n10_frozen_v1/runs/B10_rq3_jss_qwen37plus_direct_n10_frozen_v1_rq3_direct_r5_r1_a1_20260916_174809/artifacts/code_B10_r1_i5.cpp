#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

// For compatibility with MSVC and Eigen, we'll use double precision throughout
using namespace Eigen;

// INVERD implementation in C++17
// This is a direct translation of the Fortran 77 subroutine to modern C++
// using Eigen for matrix operations, but preserving the exact algorithm logic
void INVERD(int ndim, MatrixXd& A, int n, MatrixXd& B, int m, 
            double& determinant, int& ising, MatrixXi& index) {
    // Initialize
    const double EPSI = 1.0e-18; // Use 1e-18 as in original Fortran comment for double precision
    determinant = 1.0;
    ising = 1;
    
    // Initialize index(1:n,3) = 0 → set third column (index 2) to zero
    index.block(0, 2, n, 1).setZero();
    
    // Initialize index for identity mapping
    for (int i = 0; i < n; ++i) {
        index(i, 0) = i; // row mapping
        index(i, 1) = i; // column mapping
    }
    
    // Main loop: i from 0 to n-1 (0-indexed equivalent of Fortran i=1..n)
    for (int i = 0; i < n; ++i) {
        // Search for pivot: find element with maximum absolute value in column i, rows i to n-1
        double amax = 0.0;
        int irow = -1;
        
        for (int j = i; j < n; ++j) {
            double abs_val = std::abs(A(j, i));
            if (abs_val > amax) {
                amax = abs_val;
                irow = j;
            }
        }
        
        if (irow == -1) {
            ising = 2;
            determinant = 0.0;
            return;
        }
        
        // Store row mapping: INDEX(I,1) = IROW
        index(i, 0) = irow;
        
        // Interchange rows to put pivot on diagonal
        if (irow != i) {
            determinant = -determinant;
            
            // Swap rows irow and i in A
            A.row(irow).swap(A.row(i));
            
            // Swap rows in B if present
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    double swap_val = B(irow, l);
                    B(irow, l) = B(i, l);
                    B(i, l) = swap_val;
                }
            }
        }
        
        // Divide pivot row by pivot element
        double pivot = A(i, i);
        determinant *= pivot;
        
        if (std::abs(pivot) < EPSI) {
            ising = 2;
            determinant = 0.0;
            return;
        }
        
        // Set pivot element to 1.0 and divide entire row by pivot
        A(i, i) = 1.0;
        A.row(i) /= pivot;
        
        // Divide corresponding elements in B
        if (m > 0) {
            for (int l = 0; l < m; ++l) {
                B(i, l) /= pivot;
            }
        }
        
        // Reduce non-pivot rows: eliminate column i in all other rows
        for (int l1 = 0; l1 < n; ++l1) {
            if (l1 == i) continue;
            
            double t = A(l1, i);
            A(l1, i) = 0.0;
            
            if (std::abs(t) < EPSI) continue;
            
            // Update row l1: A(l1,:) = A(l1,:) - t * A(i,:)
            A.row(l1) -= t * A.row(i);
            
            // Update B if present
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    B(l1, l) -= t * B(i, l);
                }
            }
        }
    }
    
    // Check if matrix is singular
    ising = 1;
    for (int k = 0; k < n; ++k) {
        if (std::abs(A(k, k)) < EPSI) {
            ising = 2;
            break;
        }
    }
    
    // Apply column permutations to solution if needed
    // For partial pivoting, only row permutations matter for solution extraction
    // So solution is already correct in B
}

int main() {
    // Test case setup
    const int n = 4;
    const int ndim = 4;
    
    // K matrix: [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
    MatrixXd K(n, n);
    K << 200.0, -100.0,   0.0,   0.0,
        -100.0,  200.0, -100.0,   0.0,
           0.0, -100.0,  200.0, -100.0,
           0.0,    0.0, -100.0,  100.0;
    
    // b = K * [1,2,3,4]^T
    VectorXd x_true(n);
    x_true << 1.0, 2.0, 3.0, 4.0;
    VectorXd b = K * x_true;
    
    // Create copies for INVERD
    MatrixXd A = K;
    MatrixXd B(ndim, 1);
    B.col(0) = b;
    
    // Allocate workspace
    MatrixXi index(n, 3);
    
    // Call INVERD
    double determinant = 0.0;
    int ising = 1;
    
    INVERD(ndim, A, n, B, 1, determinant, ising, index);
    
    // Extract solution vector
    VectorXd solution = B.col(0);
    
    // Compute inverse: A now contains the inverse after INVERD
    MatrixXd K_inv = A;
    
    // Output plain-text result (no JSON)
    std::cout << "INVERD\n";
    std::cout << "determinant: " << std::setprecision(16) << determinant << "\n";
    std::cout << "ising: " << ising << "\n";
    std::cout << "solution:";
    for (int i = 0; i < n; ++i) {
        std::cout << " " << std::setprecision(16) << solution(i);
    }
    std::cout << "\n";
    std::cout << "inverse:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << std::setprecision(16) << K_inv(i, j);
            if (j < n-1) std::cout << " ";
        }
        std::cout << "\n";
    }
    
    return 0;
}