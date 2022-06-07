#include<cmath>
#include<iostream>
#include<vector>

using namespace std ;

int main(){

    vector<vector<double> >a(2) , b(2) , c(2) ;
    for (int i = 0; i < 2; i++)
    {
        a[i]. resize(2);
        b[i]. resize(2);
        c[i]. resize(2);

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

    for(int i=0;i<2;i++){
        for (int j = 0; j < 2; j++)
        {
            c[i][j] = a[i][j] + b[i][j] ;
            cout<<" "<<c[i][j]<<" " ;
        }
        cout<<"\n" ;
    }
    cout<<"\n" ;

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            c[i][j] = 0. ;
        }
        
    }
    

    cout<<" a - b = "<<"\n" ;

    for(int i=0;i<2;i++){
        for (int j = 0; j < 2; j++)
        {
            /* code */
            c[i][j] = a[i][j] - b[i][j] ;
            cout<<" "<<c[i][j]<<" " ;
        }
        cout<<"\n" ;
    }
    cout<<"\n" ;

    /*for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            c[i][j] = 0. ;
        }
        
    }*/

    cout<<" a b = "<<"\n" ;

    for(int i=0;i<2;i++){
        for (int j = 0; j < 2; j++){
            
            double d = 0. ;
            for (int k = 0; k < 2; k++){
                d += a[i][k]*b[k][j];
            }
            c[i][j] = d ;
            cout<<" "<<c[i][j]<<" " ;
        
        }
        cout<<"\n" ;
    }

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            c[i][j] = 0. ;
        }
        
    }


    return 0 ;
};