#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    double a = 1.;
    double b = -4. ;
    double c = -1. ;
    double d = 4. ;

    double x_init = 10.;
    int iter_max = 20 ;
    double threshold = 1e-6 ;

    double x_1 = x_init ;
    double fx_1 ;
    double dfx_1 ;

    double x_2 ;

    for (int i = 0; i < iter_max; i++)
    {
        fx_1 = a*(x_1*x_1*x_1)+b*(x_1*x_1)+c*x_1+d ;
        dfx_1 = 3*a*(x_1*x_1)+2*b*x_1+c ;

        x_2 = x_1 - fx_1/dfx_1;

        if (abs(x_2-x_1)<threshold)break;

        x_1 = x_2 ;
    } ;
    
    cout<<x_2<<"\n" ;



    return 0 ;
};