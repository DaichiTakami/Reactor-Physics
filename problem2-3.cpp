#include<cmath>
#include<iostream>
#include<vector>

using namespace std ;

int main(){

    int N = 2 ;

    vector<vector<double> >a(N) , b(N) , c(N);
    for (int i = 0; i < N; i++)
    {
        a[i]. resize(N);
        b[i]. resize(N);
        c[i]. resize(N);
    };
    

    a[0][0] = 3. ;
    a[0][1] = -1. ;
    a[1][0] = 5. ;
    a[1][1] = 6. ;

    b[0][0] = 6. ;
    b[1][0] = -1. ;
    b[0][1] = 4. ;
    b[1][1] = 3. ;

    cout<<" a + b = "<<"\n" ;

    for(int i=0;i<N;i++){
        for (int j = 0; j < N; j++)
        {
            c[i][j] = a[i][j] + b[i][j] ;
            cout<<" "<<c[i][j]<<" " ;
        }
        cout<<"\n" ;
    }
    cout<<"\n" ;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            c[i][j] = 0. ;
        }
        
    }

    cout<<" a - b = "<<"\n" ;

    for(int i=0;i<N;i++){
        for (int j = 0; j < N; j++)
        {
            c[i][j] = a[i][j] - b[i][j] ;
            cout<<" "<<c[i][j]<<" " ;
        }
        cout<<"\n" ;
    }

    cout<<"\n" ;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            c[i][j] = 0. ;
        }
        
    }

    cout<<" a b = "<<"\n" ;

    for(int i=0;i<N;i++){
        for (int j = 0; j < N; j++){
            for (int k = 0; k < N; k++){
                c[i][j] += a[i][k]*b[k][j];
            }
            cout<<" "<<c[i][j]<<" " ;
        
        }
        cout<<"\n" ;
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            c[i][j] = 0. ;
        }
        
    }


    return 0 ;
};