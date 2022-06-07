#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    double a = 1. ;
    double b = -1. ;
    double c = -1.9 ;

    double x_init = -5 ;
    double x_fin = 5. ;

    double dx = 0.01 ;

    int N = (x_fin - x_init)/dx ;

    double fx_1, fx_2, x_2 ;

    double x_1 = x_init ;

    fx_1 = a*(x_1*x_1)+b*x_1+c ;

    for (int i = 1 ; i <= N ; i++){

        x_2 = x_1 + dx ;

        fx_2 = a*(x_2*x_2)+b*x_2+c ;

        if(fx_1*fx_2 < 0){
            cout<<"["<<x_1<<" , "<<x_2<<"]"<<"\n" ;
        }

        x_1 = x_2 ;
        fx_1 = fx_2 ;
    }

    return 0 ;
};