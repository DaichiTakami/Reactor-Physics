#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<time.h>

using namespace std;


int main(){

    srand(10) ;

    double mu = 5.38 ;
    double sigma = 0.55 ;

    double x = rand()/RAND_MAX ;
    double y = rand()/RAND_MAX ;

    double z_1 = sqrt(-2*log(x))*cos(2*M_PI*y) ;        
    double z_2 = sqrt(-2*log(x))*sin(2*M_PI*y) ;

    z_1 = z_1*sigma + mu ;
    z_2 = z_2*sigma + mu ;
    
            


    return 0 ;
};