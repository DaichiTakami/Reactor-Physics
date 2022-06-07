#include<cmath>
#include<iostream>

using namespace std ;

int main(){

    double a = 1. ;
    double b = -1. ;
    double c = -2.0 ;

    double x = -5 ;
    for (int i = 1 ; i <= 1000 ; i++){
        double fx = a*x*x+b*x+c ;
        if(abs(fx) < 0.015){
            cout<<x<<"\n" ;
        };
    x+=0.01;
    };

    return 0 ;
};