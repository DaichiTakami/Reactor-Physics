#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>

using namespace std ;

int random(int fileNum){
    return rand()%fileNum ;
}

int main(){

    ifstream ifs ;

    const char *inputfile = "covariance.dat" ;
    ofstream ofs(inputfile) ;

    const int fileNum = 10000 ;
    const int sampleNum = 100 ;
    const int cal = fileNum/sampleNum ;

    char filename[256] ;
    char buf[256] ;

    const int targetNum = 21 ;
    const int mockupNum = 2 ;

    int target[] = {0,0.1,1.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0} ;

    double kt[fileNum] ;

    double ktSum ;
    vector<double> kmSum(mockupNum) ;

    double ktSum2 ;
    double kmSum2[mockupNum] ;

    double ktSum4 ;
    vector<double> kmSum4(mockupNum) ;

    double ktkmSum[mockupNum] ;
    double ktkmSum2[mockupNum] ;

    double mu_t_cas[targetNum][cal] ;
    double mu_t_est_cas[targetNum][cal] ;
    double std_t_cas[targetNum][cal] ;
    double std_t_est_cas[targetNum][cal] ;

    for (int ii = 0; ii < targetNum; ii++)
    {
        ktSum[ii] = 0. ;
        ktSum2[ii] = 0. ;
        ktSum4[ii] = 0. ;
    }
    
    for (int ii = 0; ii < mockupNum; ii++)
    {
        kmSum[ii] = 0. ;
        kmSum4[ii] = 0. ;
    }

    for (int ii = 0; ii < mockupNum; ii++)
    {
        kmSum2[ii] = 0. ;
        
    }
    

    for (int ii = 0; ii < targetNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            ktkmSum[ii][jj] = 0. ;
            ktkmSum2[ii][jj] = 0. ;
        }
        
    }
    
    
    ifs.open("SAMPLE/sample_ref");

    int linenum = 0 ;

    while (ifs.getline(buf,sizeof(buf)))
    {
        linenum ++ ;
    }
    
    // ofs<<"line number = "<<linenum<<"\n" ;

    ifs.clear() ;
    ifs.seekg(0, std::ios::beg) ;

    double *ref ;
    ref = new double[linenum] ;

    for (int jj = 0; jj < linenum; jj++)
    {
        ifs.getline(buf,sizeof(buf)) ;
        ref[jj] = atof(buf) ;
    }

    ifs.close() ;

    // ofs<<linenum ;

    for (int ii = 0; ii < fileNum; ii++)
    {

        sprintf(filename, "SAMPLE/sample_%d", ii);

        // ofs<<filename<<"\n" ;

        ifs.open(filename) ;
        if (ifs.fail())
        {
            ofs<<"can not open file."<<"\n" ;
            exit(1);
        }

        double *arr ;
        arr = new double[linenum] ;

        for (int jj = 0; jj < linenum; jj++)
        {
            ifs.getline(buf,sizeof(buf)) ;
            arr[jj] = atof(buf) ;
        }

        ifs.close() ;

        for (int jj = 0; jj < linenum; jj++)
        {
            arr[jj] /= ref[jj] ;
            arr[jj] -= 1. ;
        }
        
        for (int jj = 0; jj < targetNum; jj++)
        {
            ktSum[jj] += arr[jj] ;
            ktSum2[jj] += arr[jj]*arr[jj] ;
            ktSum4[jj] += pow(arr[jj],4) ;
        }
        
        for (int jj = 0; jj < mockupNum; jj++)
        {
            kmSum[jj] += arr[jj+targetNum] ;
            kmSum2[jj] += arr[jj+targetNum]*arr[jj+targetNum] ;
            kmSum4[jj] += pow(arr[jj+targetNum],4) ;
        }
        
        for (int jj = 0; jj < targetNum; jj++)
        {
            for (int kk = 0; kk < mockupNum; kk++)
            {
                ktkmSum[jj][kk] += arr[jj]*arr[kk+targetNum] ;
                ktkmSum2[jj][kk] += pow(arr[jj]*arr[kk+targetNum],2) ;
            }   
        }
    }

    double mu_tm[targetNum][mockupNum] ;
    double mu_tm2[targetNum][mockupNum] ;

    for (int ii = 0; ii < targetNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            mu_tm[ii][jj] = ktkmSum[ii][jj]/fileNum ;
            mu_tm2[ii][jj] = ktkmSum2[ii][jj]/fileNum ;
        }
    }
    
    double mu_t[targetNum] ;
    double mu_m[mockupNum] ;

    double mu_t2[targetNum] ;
    double mu_m2[targetNum] ;

    for (int ii = 0; ii < targetNum; ii++)
    {
        mu_t[ii] = ktSum[ii]/fileNum ;
        mu_t2[ii] = ktSum2[ii]/fileNum ;
    }

    for (int ii = 0; ii < mockupNum; ii++)
    {
        mu_m[ii] = kmSum[ii]/fileNum ;
        mu_m2[ii] = kmSum2[ii]/fileNum ;
    }

    double cov_tm[targetNum][mockupNum];
    double var_t[targetNum] ;
    double var_m[mockupNum] ;
    
    for (int ii = 0; ii < targetNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            cov_tm[ii][jj] = (mu_tm[ii][jj] - mu_t[ii]*mu_m[jj])*fileNum/(fileNum-1) ;
        }
        
    }
    
    for (int ii = 0; ii < targetNum; ii++)
    {
        var_t[ii] = (mu_t2[ii]-mu_t[ii]*mu_t[ii])*fileNum/(fileNum-1) ;
    }
    
    for (int ii = 0; ii < mockupNum; ii++)
    {
        var_m[ii] = (mu_m2[ii]-mu_m[ii]*mu_m[ii])*fileNum/(fileNum-1) ;
    }
    
    double std_t[targetNum] ;
    double std_m[mockupNum] ;

    for (int ii = 0; ii < targetNum; ii++)
    {
        std_t[ii] = sqrt(var_t[ii]) ;
    }
    
    for (int ii = 0; ii < mockupNum; ii++)
    {
        std_m[ii] = sqrt(var_m[ii]) ;
    }

    double corr_tm[targetNum][mockupNum] ;
    
    for (int ii = 0; ii < targetNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            corr_tm[ii][jj] = cov_tm[ii][jj]/(std_t[ii]*std_m[jj]) ;
        }
        
    }
    

    for (int ii = 0; ii < targetNum; ii++)
    {
        cout<<target[ii]<<" "<<corr_tm[ii][0]<<" "<<corr_tm[ii][1]<<"\n" ;
    }
    
    return 0 ;
};