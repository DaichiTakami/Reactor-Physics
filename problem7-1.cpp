#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    int N = 2 ;

    vector<vector<double> > a(N) ;
    for (int i = 0; i < N; i++)
    {
        a[i].resize(N);
    }
    
    a[0][0] = 1. ;
    a[0][1] = -4. ;
    a[1][0] = -3. ;
    a[1][1] = 5. ;

    double lambda = -1. ;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
            {
                a[i][j] -= lambda ;
            }
            
        }
        
    }

    a[0][0] = -a[0][0]/a[0][1] ;
    a[0][1] = a[0][1]/a[0][1];

    cout<<" x_2 = ( "<<a[0][1]<<" , "<<a[0][0]<<" ) \n";
    

    return 0 ;
};