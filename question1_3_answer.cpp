#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>
#include<string>
#include<sstream>


using namespace std ;

int main(){
    ifstream ifs ;
    ifs.open("sigma_a.csv") ;

    int group = 0 ;
    string s ;
    istringstream ss ;

    for (;;)
    {
        if (getline(ifs,s))
        {
            group ++ ;
        }
        else{
            break;
        }
    }
    
    ifs.close() ;
    // cout<<group<<"\n" ;

    vector<double> sigma_a(group) ;
    vector<double> ex_source(group) ;
    vector< vector<double> > sigma_s(group) ;
    for (int ii = 0; ii < group; ii++)
    {
        sigma_s[ii].resize(group) ;
    }
    // sigma_s[0][0] = 0. ;
    // sigma_s[0][1] = 0.02544 ;
    // sigma_s[0][2] = 0.00056 ;
    // sigma_s[1][0] = 0. ;
    // sigma_s[1][1] = 0. ;
    // sigma_s[1][2] = 0.00655 ;
    // sigma_s[2][0] = 0. ;
    // sigma_s[2][1] = 0. ;
    // sigma_s[2][2] = 0. ;

    vector<double> sigma_r(group) ;
    vector<double> phi(group) ;
    vector<double> source(group) ;
    for (int ii = 0; ii < group; ii++)
    {
        sigma_r[ii] = 0. ;
        phi[ii] = 0. ;
        source[ii] = 0. ;
    }
    
    ifs.open("sigma_a.csv") ;
    for (int ii = 0; ii < group; ii++)
    {
        getline(ifs,s) ;
        sigma_a[ii] = stod(s) ;
        cout<<sigma_a[ii]<<"\n" ;
    }
    ifs.close() ;

    ifs.open("ex_source.csv") ;
    for (int ii = 0; ii < group; ii++)
    {
        getline(ifs,s) ;
        ex_source[ii] = stod(s) ;
        cout<<ex_source[ii]<<"\n" ;
    }
    ifs.close() ;
    
    ifs.open("sigma_s.csv") ;
    for (int ii = 0; ii < group; ii++)
    {
        for (int jj = 0; jj < group; jj++)
        {
            ifs >> sigma_s[ii][jj] ;
            cout<<sigma_s[ii][jj]<<"\n" ;
        }
        
    }
    ifs.close() ;
            
            
        
        
    // }

    for (int ii = 0; ii < group; ii++)
    {
        sigma_r[ii] = sigma_a[ii] ;
        for (int jj = ii; jj < group; jj++)
        {
            sigma_r[ii] += sigma_s[ii][jj] ;
            
        }
        cout<<sigma_r[ii]<<"\n";
    }
    
    for (int ii = 0; ii < group; ii++)
    {
        double source_t = source[ii] + ex_source[ii] ;
        cout<<source_t<<"\n" ;
        phi[ii] = source_t/sigma_r[ii] ;
        for (int jj = ii+1; jj < group; jj++)
        {
            source[jj] = source[jj] +sigma_s[ii][jj]*phi[ii] ;
        }
        
    }
    
    for (int ii = 0; ii < group; ii++)
    {
        cout<<phi[ii]<<"\n";
    }
    

    return 0 ;
}

