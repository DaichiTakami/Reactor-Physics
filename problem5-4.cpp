#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>

using namespace std;

void cal(int N ,double al , vector<double> a,vector<double> b){

    double intfx_N = 0. ;

    double intfx = atan(1/al)*2/al ;

    for (int i = 0; i < N; i++)
    {
        intfx_N += b[i]/(a[i]*a[i]+al*al) ;
    }
    
    double RE = (intfx_N - intfx)/intfx ;

    cout<<"intfx_"<<N<<" = "<<intfx_N<<"\n" ;
    cout<<"Relative Error = "<<RE<<"\n" ;
    cout<<"\n" ;


}

int main(){

    int N ;

    //alpha = 0.1

    double a = 0.1 ;
    cout<<" alpha = "<<a<<"\n";

    // 4-point 

    N = 4 ;

    vector<double> x_4 = {-0.861136,-0.339981,0.339981,0.861136};
    vector<double> w_4 = {0.347855,0.652145,0.652145,0.347855} ;

    cal(N,a,x_4,w_4);

    // 8-point 

    N = 8 ;

    vector<double> x_8 = {-0.96029, -0.796666, -0.525532, -0.183435, 0.183435, 0.525532, 0.796666, 0.96029};
    vector<double> w_8 = {0.101229, 0.222381, 0.313707, 0.362684, 0.362684, 0.313707, 0.222381, 0.101229} ;

    cal(N,a,x_8,w_8);

    // 12-point

    N = 12 ;

    vector<double> x_12 = {-0.981561, -0.904117, -0.769903, -0.587318, -0.367831, -0.125233, 0.125233, 0.367831, 0.587318, 0.769903, 0.904117, 0.981561};
    vector<double> w_12 = {0.0471753, 0.106939, 0.160078, 0.203167, 0.233493, 0.249147, 0.249147, 0.233493, 0.203167, 0.160078, 0.106939, 0.0471753} ;

    cal(N,a,x_12,w_12);

    // 16-point 
    
    N = 16 ;

    vector<double> x_16 = {-0.989401, -0.944575, -0.865631, -0.755404, -0.617876, -0.458017, -0.281604, -0.0950125, 0.0950125, 0.281604, 0.458017, 0.617876, 0.755404, 0.865631, 0.944575, 0.989401};
    vector<double> w_16 = {0.0271525, 0.0622535, 0.0951585, 0.124629, 0.149596, 0.169157, 0.182603, 0.189451, 0.189451, 0.182603, 0.169157, 0.149596, 0.124629, 0.0951585, 0.0622535, 0.0271525} ;

    cal(N,a,x_16,w_16);

    //alpha = 0.01

    a = 0.01 ;
    cout<<" alpha = "<<a<<"\n";


    // 4-point

    N = 4 ;

    cal(N,a,x_4,w_4) ;

    // 8-point

    N = 8 ;

    cal(N,a,x_8,w_8);

    // 12-point

    N = 12 ;

    cal(N,a,x_12,w_12);

    // 16-point

    N = 16 ;

    cal(N,a,x_16,w_16);

    return 0 ;
};