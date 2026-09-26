Matrix3d blockFromFlat(const MatrixXd& A9, int flat1) {
    Matrix3d B;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++)
        B(i,j) = A9((flat1-1) + i + 3*j);
    return B;
}