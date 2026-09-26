for (int J = 1; J <= 4; J++) {
    MatrixXd TJTE = E;
    MatrixXd KSUMJ = KSUM.block<3,3>((J-1)*3, (J-1)*3);
    MatrixXd TEMP18 = KSUMJ * TJTE.transpose();
    MatrixXd KOUT = TITE * TEMP18;
    ...
}