#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    int N = 3 ;

    vector<vector<double> > a(N) , d(N) , e(N) , invA(N) ;
    for (int i = 0; i < N; i++)
    {
        a[i].resize(N);
        d[i].resize(N);
        e[i].resize(N);
        invA[i].resize(N);
    }
    // a(N)は計算用、d(N)は計算時のデータの一時保存用、e(N)はAをそのまま残す用

    vector<double> b(N) ;

    double c ;

    a[0][0] = 1.;
    a[0][1] = 2.;
    a[0][2] = 1.;
    a[1][0] = 4.;
    a[1][1] = 5.;
    a[1][2] = 6.;
    a[2][0] = 1.;
    a[2][1] = 8.;
    a[2][2] = 9.;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            d[i][j] = a[i][j] ;
            e[i][j] = a[i][j] ;
            if (i == j)
            {
                invA[i][j] = 1.;
            }else{
                invA[i][j] = 0.;
            }
        }
        
    }

    b[0] = 5.;
    b[1] = 4.;
    b[2] = 3.;

    double pivot ;

    double threshold = 1e-6 ;

    
    for (int i = 0; i < N; i++)
    {
        //pivotの選択
        pivot = a[i][i] ;
        if (abs(pivot)<threshold && i<2)  
        {
            for (int k = 0; k < N; k++)
            {
                c = a[i][k] ;
                a[i][k] = a[i+1][k] ;
                a[i+1][k] = c ;
                d[i][k] = d[i+1][k] ;
                d[i+1][k] = c ;
                c = b[k] ;
                b[k] = b[k+1] ;
                b[k+1] = c ;

            }
            pivot = a[i][i] ;
                
        }

        if (abs(pivot)<threshold && i<1)  
        {
            for (int k = 0; k < N; k++)
            {
                c = a[i][k] ;
                a[i][k] = a[i+2][k] ;
                a[i+2][k] = c ;
                d[i][k] = d[i+2][k] ;
                d[i+2][k] = c ;
                c = b[k] ;
                b[k] = b[k+2] ;
                b[k+2] = c ;
            }
            pivot = a[i][i] ;
        }

        if (abs(pivot)<threshold)
        {
            cout<<"一意に定まらない"<<"\n" ;
        }



        b[i] = b[i]/pivot ;

        for (int j = i + 1; j < N; j++)
        {
            b[j] -= b[i]*d[j][i] ;
        }

        for (int j = 0; j < N; j++)
        {
            invA[i][j] = invA[i][j]/pivot ;
            for (int k = i+1; k < N; k++)
            {
                invA[k][j] -= invA[i][j]*a[k][i];
            }
            
        }
        
        
        

        for (int j = i; j < N; j++)
        {
            a[i][j] = a[i][j]/pivot ;
        }

            for (int k = i+1 ; k < N; k++)
            {
                double tmp = a[k][i] ;
                for (int j = i; j < N; j++)
                {
                    a[k][j] -= tmp*a[i][j] ;
                }
                
                //a[k][j] -= a[i][j]*d[k][i];
                //a[k][j] -= a[i][j]*a[k][i];
            }
        



        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                d[j][k] = a[j][k] ;
            }
            
        }

    }

    for (int i = N-1; i > 0; i--)
    {
        for (int j = i-1; j >-1; j--)
        {
            pivot = a[j][i];
            for (int k = 0; k < N; k++)
            {
                invA[j][k] -= pivot*invA[i][k];
            }
             
        }
        
    }
    

    for (int i = N-1; i > 0; i--)
    {
       for (int j = 0; j < i; j++)
       {
           b[j] -= a[j][i]*b[i] ;
           a[j][i] -= a[j][i]*a[j][j] ;
       }
        
    }

    cout<<" b = ( ";
    
    
    for (int i = 0; i < N; i++)
    {
        cout<<b[i];
        if (i == N-1) break;
        cout<<" , ";
        
    }

    cout<<" ) "<<"\n" ;

    cout<<" invA = "<<"\n" ;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            cout<<invA[i][j];

            if (j == N-1)
            {
                cout<<"\n" ;
            }else{
                cout<<" , ";
            }
            
        }
        
    }

    cout<<"A invA = "<<"\n";

    for(int i=0;i<N;i++){
        for (int j = 0; j < N; j++){
            a[i][j] = 0 ;
            for (int k = 0; k < N; k++){
                a[i][j] += e[i][k]*invA[k][j];
                
                
                
            }

            if (a[i][j]<threshold)
                {
                    a[i][j] = 0 ;
            }

            cout<<a[i][j] ;

            if (j == N-1)
            {
                cout<<"\n" ;
            }else{
                cout<<" , ";
            }
        
        }
    }
    
    





    return 0 ;
};