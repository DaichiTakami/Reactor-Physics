#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    double a = 2. ;
    double b = -1. ;
    double c = -1. ;
    double d = 2. ;

    double x_1 = -0.5773502 ;
    double x_2 = 0.5773502 ;

    double w_1 =  1. ;
    double w_2 =  1. ;

    double fx_1 = a*(x_1*x_1*x_1)+b*(x_1*x_1)+c*x_1+d ;
    double fx_2 = a*(x_2*x_2*x_2)+b*(x_2*x_2)+c*x_2+d ;

    double intfx = w_1*fx_1 + w_2*fx_2 ;

    cout<<intfx<<"\n" ;

    return 0 ;
};