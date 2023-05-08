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
    // ofstream ofs("question3-3_200.csv") ;

    char buf[256] ;

    ifs.open("a.csv") ;
    if (ifs.fail())
    {
        cout<<"Cannot open a.csv"<<"\n" ;
    }
    
    ifs.getline(buf,sizeof(buf)) ;
    double a = atof(buf) ;
    ifs.close() ;

    ifs.open("D.csv") ;
    ifs.getline(buf,sizeof(buf)) ;
    double D = atof(buf) ;
    ifs.close() ;

    // double sigma_a = 0.1 ;
    ifs.open("sigma_a.csv") ;
    ifs.getline(buf,sizeof(buf)) ;
    double sigma_a = atof(buf) ;
    ifs.close() ;
    // double S = 1. ;
    // double d_x = 0.25 ;

    ifs.open("meshNum.csv") ;
    ifs.getline(buf,sizeof(buf)) ;
    int meshNum = atoi(buf) ;
    ifs.close() ;

    double d_x = a/(double)meshNum ;
    vector<double> S(meshNum) ;
    // for (int ii = 0; ii < meshNum; ii++)
    // {
    //     S[ii] = 1. ;
    // }

    ifs.open("ex_source3-7.csv") ;
    for (int jj = 0; jj < meshNum; jj++)
        {
            ifs.getline(buf,sizeof(buf)) ;
            S[jj] = atof(buf) ;
            // cout<<S[jj]<<"\n" ;
        }

    ifs.close() ;
    

    vector<double> phi_1(meshNum) ;
    for (int ii = 0; ii < meshNum; ii++)
    {
        phi_1[ii] = 0. ;
    }
    
    vector<double> phi_2(meshNum) ;

    double eps ;

    int repeat = 0 ;
    int repeat_max = 10000 ;

    for (int ii = 0; ii < repeat_max; ii++)
    {
        repeat ++ ;
        for (int jj = 0; jj < meshNum; jj++)
        {
            if (jj == 0)
            {
                phi_2[jj] = (S[jj]*d_x+D*phi_1[jj+1]/d_x)/(3*D/d_x+sigma_a*d_x) ;
                eps = fabs((phi_2[jj]-phi_1[jj])/phi_1[jj]) ;
            }else if(jj == meshNum-1){
                phi_2[jj] = (S[jj]*d_x+D*phi_2[jj-1]/d_x)/(3*D/d_x+sigma_a*d_x) ;
                if (fabs((phi_2[jj]-phi_1[jj])/phi_1[jj])>eps)
                {
                    eps = fabs((phi_2[jj]-phi_1[jj])/phi_1[jj]) ;
                }
                
            }else{
                phi_2[jj] = (S[jj]*d_x+D*phi_1[jj+1]/d_x+D*phi_2[jj-1]/d_x)/(2*D/d_x+sigma_a*d_x) ;
                if (fabs((phi_2[jj]-phi_1[jj])/phi_1[jj])>eps)
                {
                    eps = fabs((phi_2[jj]-phi_1[jj])/phi_1[jj]) ;
                }
            }
            
        }
        
        if (eps<1e-5)
        {
            break;
        }

        for (int jj = 0; jj < meshNum; jj++)
        {
            phi_1[jj] = phi_2[jj] ;
        }
        
        // cout<<eps<<"\n" ;
    }

    cout<<"#repeat : "<<repeat<<"\n" ;

    
    for (int ii = 0; ii < meshNum; ii++)
    {
        cout<<(ii+0.5)*d_x<<" "<<phi_2[ii]<<"\n";
    }
    

   

}