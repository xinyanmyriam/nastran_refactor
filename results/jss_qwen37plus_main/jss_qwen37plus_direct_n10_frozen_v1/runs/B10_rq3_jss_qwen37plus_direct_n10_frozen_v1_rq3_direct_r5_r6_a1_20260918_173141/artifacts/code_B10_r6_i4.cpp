#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <array>
#include <string>
#include <sstream>

using MatrixXd = Eigen::MatrixXd;
using VectorXd = Eigen::VectorXd;

// Minimal JSON serializer (STL-only)
std::string json_escape(const std::string& s) {
    std::ostringstream oss;
    oss << '"';
    for (char c : s) {
        switch (c) {
            case '"':  oss << "\\\""; break;
            case '\\': oss << "\\\\"; break;
            case '\b': oss << "\\b";  break;
            case '\f': oss << "\\f";  break;
            case '\n': oss << "\\n";  break;
            case '\r': oss << "\\r";  break;
            case '\t': oss << "\\t";  break;
            default:
                if (static_cast<unsigned char>(c) < 0x20) {
                    oss << "\\u" << std::hex << std::setw(4) << std::setfill('0') << static_cast<int>(c);
                } else {
                    oss << c;
                }
        }
    }
    oss << '"';
    return oss.str();
}

std::string json_string(const std::string& s) {
    return json_escape(s);
}

std::string json_double(double x) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(15) << x;
    std::string s = oss.str();
    // Remove trailing zeros and potential trailing decimal point
    size_t dot_pos = s.find('.');
    if (dot_pos != std::string::npos) {
        // Trim trailing zeros
        s.erase(s.find_last_not_of('0') + 1, std::string::npos);
        if (s.back() == '.') {
            s.pop_back();
        }
    }
    return s;
}

std::string json_array(const std::vector<double>& v, int indent = 0) {
    if (v.empty()) return "[]";
    std::string indent_str(indent * 2, ' ');
    std::ostringstream oss;
    oss << "[\n";
    for (size_t i = 0; i < v.size(); ++i) {
        oss << indent_str << "  ";
        oss << json_double(v[i]);
        if (i != v.size() - 1) oss << ",";
        oss << "\n";
    }
    oss << indent_str << "]";
    return oss.str();
}

std::string json_array(const std::vector<std::vector<double>>& v, int indent = 0) {
    if (v.empty()) return "[]";
    std::string indent_str(indent * 2, ' ');
    std::ostringstream oss;
    oss << "[\n";
    for (size_t i = 0; i < v.size(); ++i) {
        oss << indent_str << "  [\n";
        for (size_t j = 0; j < v[i].size(); ++j) {
            oss << indent_str << "    " << json_double(v[i][j]);
            if (j != v[i].size() - 1) oss << ",";
            oss << "\n";
        }
        oss << indent_str << "  ]";
        if (i != v.size() - 1) oss << ",";
        oss << "\n";
    }
    oss << indent_str << "]";
    return oss.str();
}

std::string json_object(const std::vector<std::pair<std::string, std::string>>& fields, int indent = 0) {
    std::string indent_str(indent * 2, ' ');
    std::ostringstream oss;
    oss << "{\n";
    for (size_t i = 0; i < fields.size(); ++i) {
        oss << indent_str << "  " << json_string(fields[i].first) << ": " << fields[i].second;
        if (i != fields.size() - 1) oss << ",";
        oss << "\n";
    }
    oss << indent_str << "}";
    return oss.str();
}

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
    const double EPSI = (sizeof(long double) == 10 || sizeof(long double) == 16) ? 1.0e-18 : 1.0e-36;
    
    DETERM = 1.0;
    ISING = 1;

    // INDEX(N,3): index[i] = {irow, icolum, used}
    std::vector<std::array<int, 3>> INDEX(N);
    for (int j = 0; j < N; ++j) {
        INDEX[j][2] = 0; // used flag
    }

    // Main loop over pivot steps
    for (int i = 0; i < N; ++i) {
        // Search for pivot: max |A(j,k)| where j,k not yet used
        double AMAX = 0.0;
        int IROW = -1, ICOLUM = -1;

        for (int j = 0; j < N; ++j) {
            if (INDEX[j][2] == 1) continue;
            for (int k = 0; k < N; ++k) {
                if (INDEX[k][2] == 1) continue;
                double abs_val = std::abs(A(j, k));
                if (abs_val > AMAX) {
                    AMAX = abs_val;
                    IROW = j;
                    ICOLUM = k;
                }
            }
        }

        if (IROW == -1 || ICOLUM == -1 || AMAX < EPSI) {
            ISING = 2;
            return;
        }

        // Mark column as used
        INDEX[ICOLUM][2] = 1;
        INDEX[i][0] = IROW; // store row permutation
        INDEX[i][1] = ICOLUM; // store column permutation

        // Interchange rows IROW and ICOLUM to bring pivot to (ICOLUM, ICOLUM)
        if (IROW != ICOLUM) {
            DETERM = -DETERM;
            A.row(IROW).swap(A.row(ICOLUM));
            if (M > 0) {
                B.row(IROW).swap(B.row(ICOLUM));
            }
        }

        // Pivot element
        double PIVOT = A(ICOLUM, ICOLUM);
        DETERM *= PIVOT;

        if (std::abs(PIVOT) < EPSI) {
            ISING = 2;
            return;
        }

        // Normalize pivot row
        A.row(ICOLUM) /= PIVOT;
        if (M > 0) {
            B.row(ICOLUM) /= PIVOT;
        }

        // Eliminate other rows
        for (int l1 = 0; l1 < N; ++l1) {
            if (l1 == ICOLUM) continue;
            double T = A(l1, ICOLUM);
            A(l1, ICOLUM) = 0.0;
            if (std::abs(T) < EPSI) continue;

            A.row(l1) -= T * A.row(ICOLUM);
            if (M > 0) {
                B.row(l1) -= T * B.row(ICOLUM);
            }
        }
    }

    // Now apply column interchanges (reverse order) to get true inverse
    for (int i = N-1; i >= 0; --i) {
        int L_idx = N - 1 - i;
        if (INDEX[L_idx][0] == INDEX[L_idx][1]) continue;

        int JROW = INDEX[L_idx][0];
        int JCOLUM = INDEX[L_idx][1];

        // Swap columns JROW and JCOLUM
        A.col(JROW).swap(A.col(JCOLUM));
    }

    // Check if all columns were used (i.e., matrix is non-singular)
    for (int k = 0; k < N; ++k) {
        if (INDEX[k][2] != 1) {
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
    MatrixXd B(4, 1);
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
    std::vector<std::pair<std::string, std::string>> fields;

    fields.emplace_back("test", json_string("INVERD"));
    
    std::vector<double> sol_vec = {solution(0), solution(1), solution(2), solution(3)};
    fields.emplace_back("solution", json_array(sol_vec, 0));
    
    fields.emplace_back("determinant", json_double(DETERM));

    // Convert inverse to vector of vectors
    std::vector<std::vector<double>> inv_vec;
    for (int i = 0; i < 4; ++i) {
        std::vector<double> row;
        for (int j = 0; j < 4; ++j) {
            row.push_back(K_inv(i, j));
        }
        inv_vec.push_back(row);
    }
    fields.emplace_back("inverse", json_array(inv_vec, 0));

    std::cout << json_object(fields, 0) << std::endl;

    return 0;
}