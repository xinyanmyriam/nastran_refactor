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
        // Store row and column mappings: INDEX(I,1) = IROW, INDEX(I,2) = ICOLUM
        index(i, 0) = irow;
        index(i, 1) = icolum;
        
        // Interchange rows to put pivot on diagonal
        if (irow != icolum) {
            determinant = -determinant;
            
            // Swap rows irow and icolum in A
            A.row(irow).swap(A.row(icolum));
            
            // Swap rows in B if present
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
        
        // Set pivot element to 1.0
        A(icolum, icolum) = 1.0;
        
        // Divide entire row icolum by pivot
        A.row(icolum) /= pivot;
        
        // Divide corresponding elements in B
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
            
            // Update row l1: A[l1] -= t * A[icolum]
            A.row(l1) -= t * A.row(icolum);
            
            // Update B if present
            if (m > 0) {
                for (int l = 0; l < m; ++l) {
                    B(l1, l) -= t * B(icolum, l);
                }
            }
        }
    }
    
    // Interchange columns
    for (int i = 0; i < n; ++i) {
        int l = n - 1 - i; // l = N + 1 - I in Fortran (1-indexed), so l = n-1-i in 0-indexed
        
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
    // Test case: 4x4 SPD matrix K
    // K = [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
    MatrixXd K(4, 4);
    K << 200, -100,   0,   0,
         -100,  200, -100,   0,
            0, -100,  200, -100,
            0,    0, -100,  100;
    
    // b = K * [1,2,3,4]^T
    VectorXd x_true(4);
    x_true << 1, 2, 3, 4;
    VectorXd b = K * x_true;
    
    // Prepare matrices for INVERD
    // A will be modified in-place (copy of K)
    MatrixXd A = K;
    MatrixXd B(4, 1); // b as column vector
    B.col(0) = b;
    
    int ndim = 4;
    int n = 4;
    int m = 1;
    double determinant = 0.0;
    int ising = 0;
    
    // Index matrix: n x 3
    MatrixXi index(n, 3);
    
    // Call INVERD
    INVERD(ndim, A, n, B, m, determinant, ising, index);
    
    // Extract solution vector
    VectorXd solution = B.col(0);
    
    // Compute inverse: A now contains the inverse after INVERD
    MatrixXd inverse = A;
    
    // Output results in plain text (replacing JSON)
    std::cout << "{\n";
    std::cout << "  \"test\": \"INVERD\",\n";
    std::cout << "  \"solution\": [" << std::fixed << std::setprecision(12)
              << solution(0) << ", " << solution(1) << ", " << solution(2) << ", " << solution(3) << "],\n";
    std::cout << "  \"determinant\": " << std::scientific << std::setprecision(6) << determinant << ",\n";
    std::cout << "  \"inverse\": [\n";
    for (int i = 0; i < 4; ++i) {
        std::cout << "    [";
        for (int j = 0; j < 4; ++j) {
            std::cout << std::fixed << std::setprecision(12) << inverse(i, j);
            if (j < 3) std::cout << ", ";
        }
        std::cout << "]";
        if (i < 3) std::cout << ",";
        std::cout << "\n";
    }
    std::cout << "  ]\n";
    std::cout << "}\n";

    return 0;
}