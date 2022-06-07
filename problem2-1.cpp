#include<cmath>
#include<iostream>
#include<vector>

using namespace std ;

int main(){

    vector<float>a(2);
    vector<float>b(2);
    vector<float>c(2);
    double d ;

    a[0] = 3. ;
    a[1] = 5. ;

    b[0] = 6. ;
    b[1] = -1. ;

    for(int i=0;i<2;i++){
        c[i] = a[i] + b[i] ;
        cout<<" a["<<i<<"] + b["<<i<<"] = "<<c[i]<<"\n" ;
    
    }

    for(int i=0;i<2;i++){
        c[i] = a[i] - b[i] ;
        cout<<" a["<<i<<"] - b["<<i<<"] = "<<c[i]<<"\n" ;
    }

    d = 0. ;

    for(int i=0;i<2;i++){
        d +=a[i]*b[i] ;
    }

    cout<<" ab = "<<d<<"\n" ;
    
    return 0 ;
};