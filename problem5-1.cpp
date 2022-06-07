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
    double x_3 = (x_1+x_2)/2 ;
    
    double intfx_1 = 0. ;
    double intfx_2 = 0. ;
    double intfx_3 = 0. ;

    /*左端を用いた計算*/
    for (int i = 0; i < div; i++)
    {
        double fx_1 = a*(x_1*x_1*x_1)+b*(x_1*x_1)+c*x_1+d ;
        intfx_1 += fx_1*iter;

        x_1 += iter ;

    }

    cout<<intfx_1<<"\n" ;

    /*右端を用いた計算*/
    for (int i = 0; i < div; i++)
    {
        double fx_2 = a*(x_2*x_2*x_2)+b*(x_2*x_2)+c*x_2+d;
        intfx_2 += fx_2*iter;

        x_2 += iter ;
    }

    cout<<intfx_2<<"\n";

    /*中央を用いた計算*/
    for (int i = 0; i < div; i++)
    {
        double fx_3 = a*(x_3*x_3*x_3)+b*(x_3*x_3)+c*x_3+d;
        intfx_3 += fx_3*iter;

        x_3 += iter ;
    }
    
    cout<<intfx_3<<"\n";

    
    

    

    return 0 ;
};