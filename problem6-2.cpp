#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    int N = 3 ;

    vector<double> x_1(N) ;
    for (int i = 0; i < N; i++)
    {
        x_1[i] = 0. ;
    }

    vector<double> x_2(N) ;
    
    vector<vector<double> > a(N) ;
    for (int i = 0; i < N; i++)
    {
        a[i].resize(N) ;
    }
    a[0][0] = 10. ;
    a[0][1] = 2. ;
    a[0][2] = 0. ;
    a[1][0] = 4. ;
    a[1][1] = 10. ;
    a[1][2] = 6. ;
    a[2][0] = 0. ;
    a[2][1] = 8. ;
    a[2][2] = 10. ;
    
    vector<double> b(N) ;
    b[0] = 5. ;
    b[1] = 4. ;
    b[2] = 3. ;

    int iter_max = 50 ;
    double threshold = 1e-6 ;

    double e ;

    int count = 0 ;

    for (int i = 0; i < iter_max; i++)
    {
        e = 0. ;

        for (int j = 0; j < N; j++)
        {
            double x_1 = x_2[j] ;

            x_2[j] = b[j] ;

            for (int k = 0; k < j; k++)
            {
                x_2[j] -= a[j][k]*x_2[k];
            }
            
            for (int k = j+1; k < N; k++)
            {
                x_2[j] -= a[j][k]*x_2[k];
            }

            x_2[j] = x_2[j]/a[j][j] ;

            e += abs(x_2[j] - x_1) ;


            
        }

        cout<<x_2[0]<<" , "<<x_2[1]<<" , "<<x_2[2]<<"\n" ;


        if (e < threshold)break;
        
        for (int j = 0; j < N; j++)
        {
            x_1[j] = x_2[j] ;
        }


        
        count += 1 ;

        
    }

    cout<<count<<"回"<<"\n" ;




    



    

    return 0 ;
};