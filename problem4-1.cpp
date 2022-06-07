#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;



int main(){

    double a = 1. ;
    double b = -4. ;
    double c = -1. ;
    double d = 4. ;

    double x_i = -3. ;
    double dx_i = 3. ;

    int iter_max = 20 ;

    double threshold = 1e-6 ;

    double x_1 = x_i ;
    double fx_1 = a*(x_1*x_1*x_1)+b*(x_1*x_1)+c*x_1+d ;

    double x_2 = x_1 + dx_i ;
    double fx_2 = a*(x_2*x_2*x_2)+b*(x_2*x_2)+c*x_2+d ;

    if (fx_1*fx_2>0)
    {
        cout<<"範囲を変更してください"<<"\n";
        exit(0);
    }
    
    if (abs(fx_1)<threshold)
    {
        cout<<" x = "<<x_1<<"\n";
        exit(0);
    }
    if (abs(fx_2)<threshold)
    {
        cout<<" x = "<<x_2<<"\n";
        exit(0);
    }

    double x_3 ;
    double fx_3 ;

    for (int i = 0; i < iter_max; i++ )
    {
        x_3 = (x_1+x_2)/2 ;
        fx_3 = a*(x_3*x_3*x_3)+b*(x_3*x_3)+c*x_3+d ;

        if (fx_1*fx_3>0)
        {
            x_1 = x_3 ;
            fx_1 = fx_3;
        }else if (fx_2*fx_3>0)
        {
            x_2 = x_3 ;
            fx_2 = fx_3;
        }
        
        if (abs(fx_3)<threshold)break;


        
        
        

    }

    cout<<" x = "<<x_3<<"\n";
    

    

    return 0 ;
};