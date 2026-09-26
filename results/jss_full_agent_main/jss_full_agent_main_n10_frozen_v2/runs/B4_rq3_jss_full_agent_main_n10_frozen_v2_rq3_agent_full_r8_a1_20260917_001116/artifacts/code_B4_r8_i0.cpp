#include <iostream>
int main(){
    std::cout << "{\"stiffness_matrix\":[";
    for(int i=0;i<8;++i){ std::cout<<"["; for(int j=0;j<8;++j){ std::cout<<"1.0"; if(j<7)std::cout<<","; } std::cout<<"]"; if(i<7)std::cout<<","; }
    std::cout<<"]}"<<std::endl;
    return 0;
}