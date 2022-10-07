#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>

using namespace std ;

const double PI = acos(-1);

const double sigma = 0.55 ;
const double ave = 5.38 ;

int main(){

    srand(10);

    int num = 10 ;

    vector<double> a(num) ;
    double b = 0. ;
    double v = 0. ;
    double sk = 0. ;
    double ku = 0. ;

    int c = 10000 ;

    double range = 6*sigma ;

    double dx = range/c ;
    
    for (int i = 0; i < num; i++)
    {
        double x = (double)rand()/RAND_MAX ;
        double y = (double)rand()/RAND_MAX ;

        double w = sqrt(-2*log(x))*cos(2*PI*y);
        double z = sigma*w + ave ;

        a[i] = z ;
        b += z ;

        // cout<<z<<"\n";

        


        
    }
    
    b = b/num ;

    for (int i = 0; i < num; i++){

        double x = a[i]-b ;
        
        v += x*x/(num - 1) ;
        sk += (x/sigma)*(x/sigma)*(x/sigma) *num/((num-1)*(num-2)) ;
        ku += (x/sigma)*(x/sigma)*(x/sigma)*(x/sigma) * num * (num+1)/((num-1)*(num-2)*(num-3));

    }

    ku = ku - 3*(num-1)*(num-1)/((num-2)*(num-3)) ;




    
    

        


    // cout<<z_1<<"\n";
    // cout<<z_2<<"\n";



    
    return 0 ;
};