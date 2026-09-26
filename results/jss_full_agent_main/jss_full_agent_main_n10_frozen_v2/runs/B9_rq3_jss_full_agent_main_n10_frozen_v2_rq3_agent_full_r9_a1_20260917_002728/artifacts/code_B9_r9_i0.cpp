for(int I=1; I<=72; I+=12){
    int LOW=I, LIM=LOW+5;
    for(int K=LOW; K<=LIM; K++){
        J++;
        KEP(J-1) = KEfull(K-1);
        KEP(J+36-1) = KEfull(K+6-1);
    }
}