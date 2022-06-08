#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<time.h>

using namespace std;



int main(){

    srand(time(NULL)) ;

    int a = 0 ;
    int b = 0 ;
    int c = 0 ;
    int d = 0 ;
    int e = 0 ;
    int f = 0 ;
    int x;
    int max = 10000 ;

    for (int i = 0; i < max ; i++)
    {
        x = rand()%6 + 1 ;

        cout<<x<<"\n";

        if (x==1)
        {
            a++ ;
        }
        if (x==2)
        {
            b++ ;  
        }
        if (x==3)
        {
            c ++ ;
        }
        if (x==4)
        {
            d ++ ;
        }
        if(x==5){
            e ++ ;
        }
        if(x==6)
        {
            f ++ ;
        }
        
        
        
        
        
        
        
        
        
        
    }
    
    cout<<" 1 : "<<a<<" 回 "<<"\n" ;
    cout<<" 2 : "<<b<<" 回 "<<"\n" ;
    cout<<" 3 : "<<c<<" 回 "<<"\n" ;
    cout<<" 4 : "<<d<<" 回 "<<"\n" ;
    cout<<" 5 : "<<e<<" 回 "<<"\n" ;
    cout<<" 6 : "<<f<<" 回 "<<"\n" ;


    return 0 ;
};