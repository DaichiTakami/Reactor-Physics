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

    double x_i = -1. ;
    double x_f = 1. ;

    int div = 20 ;

    double iter = (x_f - x_i)/div ;

    double x_1 = x_i;
    double x_2 = x_1 + iter ;
    
    double intfx = 0. ;

    for (int i = 0; i < div; i++)
    {
        double fx_1 = a*(x_1*x_1*x_1)+b*(x_1*x_1)+c*x_1+d ;
        double fx_2 = a*(x_2*x_2*x_2)+b*(x_2*x_2)+c*x_2+d ;

        intfx += (fx_1+fx_2)*iter/2 ;

        x_1 += iter ;
        x_2 += iter ;
    }

    cout<<intfx<<"\n";
    

    return 0 ;
};