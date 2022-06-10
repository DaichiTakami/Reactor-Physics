#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<time.h>

using namespace std;

double average ;

int main(){

    srand(time(NULL)) ;
    
    int dice_count = 5 ;
    double range = 1./dice_count ;
    int repeat = 10000 ;

    int N = 5*dice_count +1 ;

    vector<int>a(N) ;
    int x  ;
    for (int i = 0; i < repeat; i++)
    {
        x = 0 ;
        for (int j = 0; j < dice_count; j++)
        {
            x += rand()%6 + 1 ;
        }
        double y = (double)x/dice_count ;
        // cout<<y<<"\n";

        for (int j = 0; j < N; j++)
        {
            if (abs(1 + range*j - y) < 1e-6)
            {
                a[j] += 1 ;
                break;
            }
            
        }
                
    }
    for (int j = 0; j < N; j++)
    {
        cout<<" "<<1+j*range<<" : "<<a[j]<<" 回"<<"\n"  ;
    }
        

    return 0 ;
};