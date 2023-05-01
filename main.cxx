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

    char buf[256] ;
    
    ifs.open("meshNum.csv") ;
    ifs.getline(buf,sizeof(buf)) ;
    int meshNum = atoi(buf) ;
    ifs.close() ;

    // double d_x = a/(double)meshNum ;
    vector<double> D(meshNum);
    vector<double> sigma_a(meshNum) ;
    vector<double> S(meshNum) ;

    ifs.open("a.csv") ;
    if (ifs.fail())
    {
        cout<<"Cannot open a.csv"<<"\n" ;
    }
    ifs.getline(buf,sizeof(buf)) ;
    double a = atof(buf) ;
    ifs.close() ;

    double d_x = a/(double)meshNum ;

    ifs.open("D.csv") ;
    for (int jj = 0; jj < meshNum; jj++)
    {
        ifs.getline(buf,sizeof(buf)) ;
        D[jj] = atof(buf) ;
    }
    ifs.close() ;

    ifs.open("sigma_a.csv") ;
    for (int jj = 0; jj < meshNum; jj++)
    {
        ifs.getline(buf,sizeof(buf)) ;
        sigma_a[jj] = atof(buf) ;
    }
    ifs.close() ;

    ifs.open("ex_source4-1.csv") ;
    for (int jj = 0; jj < meshNum; jj++)
    {
        ifs.getline(buf,sizeof(buf)) ;
        S[jj] = atof(buf) ;
    }
    ifs.close() ;

    vector<double> phi_1(meshNum) ;
    for (int ii = 0; ii < meshNum; ii++)
    {
        phi_1[ii] = 0. ;
    }

    int repeat = 0 ;
    int repeat_max = 10000 ;

    for (int ii = 0; ii < repeat_max; ii++)
    {
        double temp = 0. ;
        double eps = 0. ;
        double eps_max = 0. ;

        repeat ++ ;

        for (int jj = 0; jj < meshNum; jj++)
        {
            if (jj == 0)
            {
                temp = (S[jj]*d_x+D[jj]*phi_1[jj+1]/d_x)/(3*D[jj]/d_x+sigma_a[jj]*d_x) ;
            }else if(jj == meshNum-1){
                temp = (S[jj]*d_x+D[jj]*phi_1[jj-1]/d_x)/(3*D[jj]/d_x+sigma_a[jj]*d_x) ;
            }else{
                double cn_1 = (2*D[jj]*D[jj+1])/(d_x*D[jj+1]+d_x*D[jj]) ;
                double cn_2 = (2*D[jj]*D[jj-1])/(d_x*D[jj-1]+d_x*D[jj]) ;
                temp = (cn_1*phi_1[jj+1]+cn_2*phi_1[jj-1]+S[jj]*d_x)/(cn_1+cn_2+sigma_a[jj]*d_x) ;
            }

            if (ii!=0)
            {
                eps = fabs((temp-phi_1[jj])/phi_1[jj]) ;
                if (eps>eps_max)
                {
                    eps_max = eps ;
                }
            }
            
            phi_1[jj] = temp ;
        }
        
        if (ii!=0&&eps_max<1e-5)
        {
            break;
        }

    }

    cout<<"#repeat : "<<repeat<<"\n" ;

    
    for (int ii = 0; ii < meshNum; ii++)
    {
        cout<<(ii+0.5)*d_x<<" "<<phi_1[ii]<<"\n";
    }

    return 0 ;
};