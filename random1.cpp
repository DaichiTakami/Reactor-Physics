#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<time.h>

using namespace std;



int main(){

    srand(time(NULL)) ;

    cout<<RAND_MAX<<"\n" ;
    // cout<<rand()<<"\n" ;
    // cout<<random()<<"\n" ;

    int iter_max = 1000000 ;

    double a, b, pi ;

    for (int i = 0; i < iter_max; i++)
    {
        a = (double)rand()/RAND_MAX ;
        b = (double)rand()/RAND_MAX ;
        double x = a*a + b*b ;
        if (x < 1.)
        {
            pi ++ ;
        }
        // cout<<a<<"\n";        
    }
    
    pi = 4*pi/iter_max ;

    cout<<" pi = "<<pi<<"\n" ;

    return 0 ;
};