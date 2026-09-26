#include <Eigen/Dense>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace Eigen;
using namespace std;

static const double PI = 3.14159265358979323846;

MatrixXd inverd(const MatrixXd& M, int& ising) {
    FullPivLU<MatrixXd> lu(M);
    if (!lu.isInvertible()) { ising = 2; return MatrixXd::Zero(M.rows(), M.cols()); }
    ising = 1;
    return lu.inverse();
}

struct KTRBSCResult {
    MatrixXd KU[3][3];
    MatrixXd Hinv;
    MatrixXd S;
    double XSUBB, XSUBC, YSUBC;
    bool ok;