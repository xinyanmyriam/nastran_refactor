#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <map>
#include <array>
#include <limits>
#include <Eigen/Dense>

// Minimal JSON serializer — no external deps
namespace json {
    template<typename T>
    std::string to_string(const T& v);

    std::string escape(const std::string& s) {
        std::string result;
        result.reserve(s.size() * 2);
        for (char c : s) {
            switch (c) {
                case '"': result += "\\\""; break;
                case '\\': result += "\\\\"; break;
                case '\b': result += "\\b"; break;
                case '\f': result += "\\f"; break;
                case '\n': result += "\\n"; break;
                case '\r': result += "\\r"; break;
                case '\t': result += "\\t"; break;
                default:
                    if (static_cast<unsigned char>(c) < 32) {
                        std::ostringstream oss;
                        oss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c);
                        result += oss.str();
                    } else {
                        result += c;
                    }
                    break;
            }
        }
        return result;
    }

    template<>
    std::string to_string<double>(const double& v) {
        std::ostringstream oss;
        oss << std::setprecision(15) << v;
        return oss.str();
    }

    template<>
    std::string to_string<int>(const int& v) {
        return std::to_string(v);
    }

    template<>
    std::string to_string<std::string>(const std::string& v) {
        return "\"" + escape(v) + "\"";
    }

    template<typename T>
    std::string to_string(const std::vector<T>& v) {
        if (v.empty()) return "[]";
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) oss << ",";
            oss << to_string(v[i]);
        }
        oss << "]";
        return oss.str();
    }

    template<>
    std::string to_string<std::vector<std::vector<double>>>(const std::vector<std::vector<double>>& v) {
        if (v.empty()) return "[]";
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) oss << ",";
            oss << to_string(v[i]);
        }
        oss << "]";
        return oss.str();
    }

    template<>
    std::string to_string<std::vector<double>>(const std::vector<double>& v) {
        if (v.empty()) return "[]";
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < v.size(); ++i) {
            if (i > 0) oss << ",";
            oss << to_string(v[i]);
        }
        oss << "]";
        return oss.str();
    }

    // Object: map<string, string> — we'll build key-value pairs manually
    struct object {
        std::map<std::string, std::string> data;
        object() = default;
        object& operator[](const std::string& key) {
            data[key] = "";
            return *this;
        }
        object& operator=(const std::string& value) {
            if (!data.empty()) {
                data.rbegin()->second = value;
            }
            return *this;
        }
        std::string dump(int indent = 0) const {
            if (data.empty()) return "{}";
            std::ostringstream oss;
            oss << "{";
            std::string pad(indent + 2, ' ');
            bool first = true;
            for (const auto& p : data) {
                if (!first) oss << ",";
                oss << "\n" << pad << "\"" << escape(p.first) << "\":" << (p.second.empty() ? "null" : " " + p.second);
                first = false;
            }
            oss << "\n" << std::string(indent, ' ') << "}";
            return oss.str();
        }
    };
}

// Using declarations
using MatrixXd = Eigen::MatrixXd;
using VectorXd = Eigen::VectorXd;

// INVERD equivalent in C++17
// Performs full matrix inversion and/or linear solve via Gauss-Jordan elimination with full pivoting
// Returns: determinant, singularity flag, modifies A in-place to hold inverse, and B to hold solution(s)
//
// Parameters:
//   A: input/output square matrix (N x N), overwritten with inverse on success
//   N: size of matrix
//   B: right-hand side(s); if M > 0, shape (N x M), overwritten with solution(s)
//   M: number of RHS vectors (0 means no solve)
//   DETERM: output determinant
//   ISING: output singularity flag (1=non-singular, 2=singular)
//
void INVERD(MatrixXd& A, int N, MatrixXd& B, int M, double& DETERM, int& ISING) {
    // Initialize
    const double EPSI = 1.0e-12; // safe for double-precision 4x4
    
    DETERM = 1.0;
    ISING = 1;

    // Store pivot row and column permutations
    std::vector<int> IROW(N);   // IROW[i] = row index used at step i
    std::vector<int> ICOLUM(N); // ICOLUM[i] = column index used at step i
    std::vector<bool> ROW_USED(N, false);
    std::vector<bool> COL_USED(N, false);

    // Main loop over pivot steps
    for (int i = 0; i < N; ++i) {
        // Search for pivot: max |A(j,k)| where j,k not yet used
        double AMAX = 0.0;
        int best_j = -1, best_k = -1;

        for (int j = 0; j < N; ++j) {
            if (ROW_USED[j]) continue;
            for (int k = 0; k < N; ++k) {
                if (COL_USED[k]) continue;
                double abs_val = std::abs(A(j, k));
                if (abs_val > AMAX) {
                    AMAX = abs_val;
                    best_j = j;
                    best_k = k;
                }
            }
        }

        if (best_j == -1 || best_k == -1 || AMAX < EPSI) {
            ISING = 2;
            return;
        }

        IROW[i] = best_j;
        ICOLUM[i] = best_k;
        ROW_USED[best_j] = true;
        COL_USED[best_k] = true;

        // Interchange rows: bring pivot row to row i
        if (best_j != i) {
            DETERM = -DETERM;
            A.row(i).swap(A.row(best_j));
            if (M > 0) {
                B.row(i).swap(B.row(best_j));
            }
        }

        // Pivot element
        double PIVOT = A(i, best_k);
        DETERM *= PIVOT;

        if (std::abs(PIVOT) < EPSI) {
            ISING = 2;
            return;
        }

        // Normalize pivot row: divide by pivot
        A.row(i) /= PIVOT;
        if (M > 0) {
            B.row(i) /= PIVOT;
        }

        // Eliminate other rows (all rows except i)
        for (int l = 0; l < N; ++l) {
            if (l == i) continue;
            double T = A(l, best_k);
            if (std::abs(T) < EPSI) continue;
            A.row(l) -= T * A.row(i);
            if (M > 0) {
                B.row(l) -= T * B.row(i);
            }
        }
    }

    // Now apply inverse column permutation to get true inverse
    // At step i, we used column ICOLUM[i]; so the i-th column of current A corresponds to original column ICOLUM[i]
    // To get inverse in natural order, we need to place column i of result into position ICOLUM[i]
    // i.e., construct inv such that inv.col(ICOLUM[i]) = A.col(i)
    MatrixXd A_inv = MatrixXd::Zero(N, N);
    for (int i = 0; i < N; ++i) {
        A_inv.col(ICOLUM[i]) = A.col(i);
    }
    A = A_inv;

    // Also permute B's columns? No — B was solved against permuted A; but since we only solve one RHS (M=1),
    // and row swaps were applied to B, the solution is correct as-is. No column permute needed for B.
    // So B is already the solution to original system.

    // Final singularity check: all rows/columns used?
    for (int i = 0; i < N; ++i) {
        if (!ROW_USED[i] || !COL_USED[i]) {
            ISING = 2;
            return;
        }
    }
    ISING = 1;
}

int main() {
    // Test case: 4x4 SPD matrix K
    // K = [200,-100,0,0; -100,200,-100,0; 0,-100,200,-100; 0,0,-100,100]
    MatrixXd K(4, 4);
    K << 200, -100,    0,    0,
         -100,  200, -100,    0,
            0, -100,  200, -100,
            0,    0, -100,  100;

    // b = K * [1,2,3,4]^T
    VectorXd x_true(4);
    x_true << 1, 2, 3, 4;
    VectorXd b = K * x_true;

    // Make copies for INVERD (which modifies A and B)
    MatrixXd A = K;
    MatrixXd B(4, 1); // Ensure B is MatrixXd (4x1)
    B.col(0) = b;

    double DETERM;
    int ISING;

    // Call INVERD: invert A, solve A*x = B
    INVERD(A, 4, B, 1, DETERM, ISING);

    // A now holds the inverse of K
    MatrixXd K_inv = A;

    // Solution is in B (first column)
    VectorXd solution = B.col(0);

    // Build JSON output manually
    json::object result;
    result["test"] = json::to_string(std::string("INVERD"));
    
    std::vector<double> sol_vec = {solution(0), solution(1), solution(2), solution(3)};
    result["solution"] = json::to_string(sol_vec);
    
    result["determinant"] = json::to_string(DETERM);

    // Convert inverse to vector of vectors
    std::vector<std::vector<double>> inv_vec;
    for (int i = 0; i < 4; ++i) {
        std::vector<double> row;
        for (int j = 0; j < 4; ++j) {
            row.push_back(K_inv(i, j));
        }
        inv_vec.push_back(row);
    }
    result["inverse"] = json::to_string(inv_vec);

    std::cout << result.dump(2) << std::endl;

    return 0;
}