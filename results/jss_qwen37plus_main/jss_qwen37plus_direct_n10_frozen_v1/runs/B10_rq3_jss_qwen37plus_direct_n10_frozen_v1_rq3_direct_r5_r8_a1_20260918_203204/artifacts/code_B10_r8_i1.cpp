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
    const double EPSI = 1.0e-36;
    double mach = 5.0; // Simulate MACH=5 from Fortran (VAX-like)
    if (mach == 5.0) {
        // Use more realistic epsilon for double precision
        // Original Fortran used 1.D-18 for VAX, but we'll use machine epsilon scaled
        // Since we're using double, use std::numeric_limits<double>::epsilon() * 1e10
        // But per spec, use 1e-18 as in original comment
    }
    
    determinant = 1.0;
    if (ising < 0) {
        determinant = 0.0;
    }
    
    // Initialize index(1:n,3) = 0
    index.block(0, 2, n, 1).setZero();
    
    // Main loop: i from 1 to n (1-indexed in Fortran, 0-indexed in C++)
    for (int i = 0; i < n; ++i) {
        // Search for pivot: find element with maximum absolute value among unmarked rows/columns
        double amax = 0.0;
        int irow = -1;
        int icolum = -1;
        
        for (int j = 0; j < n; ++j) {
            if (index(j, 2) == 1) continue; // INDEX(J,3) == 1
            
            for (int k = 0; k < n; ++k) {
                if (index(k, 2) == 1) continue; // INDEX(K,3) == 1
                
                double abs_val = std::abs(A(j, k));
                if (abs_val > amax) {
                    amax = abs_val;
                    irow = j;
                    icolum = k;
                }
            }
        }
        
        if (irow == -1 || icolum == -1) {
            ising = 2;
            return;
        }
        
        // Mark this column as used: INDEX(ICOLUM,3) = INDEX(ICOLUM,3) + 1
        index(icolum, 2) += 1;
        
        // Store row and column permutations: INDEX(I,1) = IROW, INDEX(I,2) = ICOLUM
        index(i, 0) = irow;
        index(i, 1) = icolum;
        
        // Interchange rows to put pivot on diagonal
        if (irow != icolum) {
            determinant = -determinant;
            
            // Swap rows irow and icolum in A
            A.row(irow).swap(A.row(icolum));
            
            // Swap rows in B if m > 0
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    double swap_val = B(irow, l);
                    B(irow, l) = B(icolum, l);
                    B(icolum, l) = swap_val;
                }
            }
        }
        
        // Divide pivot row by pivot element
        double pivot = A(icolum, icolum);
        determinant *= pivot;
        
        if (std::abs(pivot) < EPSI) {
            ising = 2;
            return;
        }
        
        // Set diagonal element to 1.0 and divide entire row by pivot
        A(icolum, icolum) = 1.0;
        A.row(icolum) /= pivot;
        
        // Divide corresponding rows in B
        if (m > 0) {
            for (int l = 0; l < m; ++l) {
                B(icolum, l) /= pivot;
            }
        }
        
        // Reduce non-pivot rows
        for (int l1 = 0; l1 < n; ++l1) {
            if (l1 == icolum) continue;
            
            double t = A(l1, icolum);
            A(l1, icolum) = 0.0;
            
            if (std::abs(t) < EPSI) continue;
            
            // Update row l1: A[l1] = A[l1] - t * A[icolum]
            A.row(l1) -= t * A.row(icolum);
            
            // Update B if m > 0
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    B(l1, l) -= t * B(icolum, l);
                }
            }
        }
    }
    
    // Interchange columns (reverse order)
    for (int i = 0; i < n; ++i) {
        int l = n - 1 - i; // L = N + 1 - I in Fortran (1-indexed)
        if (index(l, 0) == index(l, 1)) continue;
        
        int jrow = index(l, 0);
        int jcolum = index(l, 1);
        
        // Swap columns jrow and jcolum in A
        A.col(jrow).swap(A.col(jcolum));
    }
    
    // Check if matrix is singular
    ising = 1;
    for (int k = 0; k < n; ++k) {
        if (index(k, 2) != 1) {
            ising = 2;
            break;
        }
    }
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
    
    // Prepare matrices for INVERD
    MatrixXd A = K; // Will be modified in-place
    MatrixXd B(ndim, 1);
    B.col(0) = b;
    
    double determinant = 0.0;
    int ising = 1;
    MatrixXi index(n, 3); // INDEX(N,3)
    
    // Call INVERD
    INVERD(ndim, A, n, B, 1, determinant, ising, index);
    
    // Extract solution vector
    VectorXd solution = B.col(0);
    
    // Compute inverse matrix (A now contains the inverse after INVERD)
    MatrixXd inverse = A;
    
    // Output in simple structured plain-text format (replaces JSON)
    std::cout << "test: INVERD" << std::endl;
    
    std::cout << "solution:";
    for (int i = 0; i < n; ++i) {
        std::cout << " " << std::fixed << std::setprecision(6) << solution(i);
    }
    std::cout << std::endl;
    
    std::cout << "determinant: " << std::fixed << std::setprecision(6) << determinant << std::endl;
    
    std::cout << "inverse:" << std::endl;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << std::fixed << std::setprecision(6) << inverse(i, j);
            if (j < n-1) std::cout << " ";
        }
        std::cout << std::endl;
    }
    
    return 0;
}