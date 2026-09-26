for(int a=0;a<3;a++){
    for(int b=0;b<3;b++){
        int ga=2*M[j][a];
        int gb=2*M[j][b];
        for(int p=0;p<2;p++){
            for(int q=0;q<2;q++){
                K[ga+p][gb+q] += Ke[2*a+p][2*b+q];
            }
        }
    }
}