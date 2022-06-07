#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    int N = 2 ;

    vector<double> x_1(N) ;
    x_1[0] = 1. ;
    x_1[2] = 1. ;

    vector<double> x_2(N) ;

    vector<vector<double> > a(N) ;
    for (int i = 0; i < N; i++)
    {
        a[i].resize(N) ;
    }
    a[0][0]= 1. ;
    a[0][1] = -4. ;
    a[1][0] = -3. ;
    a[1][1] = 5. ;

    int iter_max = 100 ;
    double threshold = 1e-6 ;

    double rambda , norm, e ;




    for (int i = 0; i < iter_max; i++)
    {

        norm = 0. ;

        e = 0. ;

        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                x_2[j] += a[j][k]*x_1[k] ;
            }
            
        }
        
        rambda = x_1[1]*x_2[1] + x_1[0]*x_2[0] ;

        for (int j = 0; j < N; j++)
        {
            norm += x_2[j]*x_2[j] ;
        }

        norm = sqrt(norm) ;

        for (int j = 0; j < N; j++)
        {
            x_2[j] = x_2[j]/norm ;
            e += abs(x_2[j]-x_1[j]) ;
            x_1[j] = x_2[j] ;
            x_2[j] = 0. ;
        }

        if (e<threshold)break;
       
        

       
        
    }

    x_1[1] = x_1[1]/x_1[0];
    x_1[0] = x_1[0]/x_1[0];

    cout<<"x_1 = ( "<<x_1[0]<<" , "<<x_1[1]<<" ) , rambda_1 = "<<rambda<<"\n" ;
    
    


    

    return 0 ;
};