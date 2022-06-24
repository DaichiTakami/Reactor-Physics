#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>

using namespace std ;

const double PI = acos(-1);

void box(double x, double y){

    double z_1 = sqrt(-2*log(x))*cos(2*PI*y); 

    
};

int main(){

    srand(10);

    double x = (double)rand()/RAND_MAX ;
    double y = (double)rand()/RAND_MAX ;

    box(x,y);

    


    // cout<<z_1<<"\n";
    // cout<<z_2<<"\n";



    
    return 0 ;
};