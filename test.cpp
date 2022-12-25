#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>

#define mean 1
#define deviation 0

#define Pu9 1
#define Pu2 0

using namespace std ;

int main(){

    ifstream ifs ;

    const char *inputfile = "Pu9_mean.dat" ;
    ofstream ofs(inputfile) ;

    const int fileNum = 10000 ;
    const int sampleNum = 100 ;
    const int cal = fileNum/sampleNum ;

    char filename[256] ;
    char buf[256] ;

    const int targetNum = 21 ;
    const int mockupNum = 2 ;

    double target[] = {0,0.1,1.0,2.5,5.0,7.5,10.0,12.5,15.0,17.5,20.0,22.5,25.0,27.5,30.0,32.5,35.0,37.5,40.0,42.5,45.0} ;

    double kt[targetNum][sampleNum] ;
    double km_i[mockupNum][sampleNum] ;

    vector<double> ktSum(targetNum) ;
    vector<double> kmSum(mockupNum) ;

    vector<double> ktSum2(targetNum) ;
    double kmSum2[mockupNum][mockupNum] ;

    vector<double> ktSum4(targetNum) ;
    vector<double> kmSum4(mockupNum) ;

    double ktkmSum[targetNum][mockupNum] ;
    double ktkmSum2[targetNum][mockupNum] ;

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
        for (int jj = 0; jj < mockupNum; jj++)
        {
            kmSum2[ii][jj] = 0. ;
        }
        
    }
    

    for (int ii = 0; ii < targetNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            ktkmSum[ii][jj] = 0. ;
            ktkmSum2[ii][jj] = 0. ;
        }
        
    }
    
    
    ifs.open("/Users/takamidaichi/Documents/GitHub/Reactor-Physics/sample_ref");

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

    for (int i = 0; i < cal; i++)
    {

        for (int ii = 0; ii < sampleNum; ii++)
        {

            sprintf(filename, "/Users/takamidaichi/Documents/GitHub/Reactor-Physics/SAMPLE/sample_%d", 100*i+ii);

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
                kt[jj][ii] = arr[jj] ;
            }
            
            for (int jj = 0; jj < mockupNum; jj++)
            {
                km_i[jj][ii] = arr[jj+targetNum] ;
            }
            
            // ofs<<km_i[1][ii]<<"\n" ;
            
            for (int jj = 0; jj < targetNum; jj++)
            {
                ktSum[jj] += arr[jj] ;
                ktSum2[jj] += arr[jj]*arr[jj] ;
                ktSum4[jj] += pow(arr[jj],4) ;
            }
            
            for (int jj = 0; jj < mockupNum; jj++)
            {
                kmSum[jj] += arr[jj+targetNum] ;
                kmSum4[jj] += pow(arr[jj+targetNum],4) ;
            }

            for (int jj = 0; jj < mockupNum; jj++)
            {
                for (int kk = 0; kk < mockupNum; kk++)
                {
                    kmSum2[jj][kk] += arr[jj+targetNum]*arr[kk+targetNum] ;
                }
                
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

        double mu_m[mockupNum] ;
        double mu_m2[mockupNum][mockupNum] ;

        double var_m2[mockupNum] ;

        double cov_mm[mockupNum][mockupNum] ;

        for (int ii = 0; ii < mockupNum; ii++)
        {
            mu_m[ii] = kmSum[ii]/sampleNum ;

            // ofs<<mu_m[ii]<<"\n";

        }

        for (int ii = 0; ii < mockupNum; ii++)
        {
            for (int jj = 0; jj < mockupNum; jj++)
            {
                mu_m2[ii][jj] = kmSum2[ii][jj]/sampleNum ;

                // ofs<<mu_m2[ii][jj]<<"\n" ;

                cov_mm[ii][jj] = (mu_m2[ii][jj] - mu_m[ii]*mu_m[jj])*sampleNum/(sampleNum-1) ;
            }
        }

        // for (int ii = 0; ii < mockupNum; ii++)
        // {
        //     for (int jj = 0; jj < mockupNum; jj++)
        //     {
        //         ofs<<cov_mm[ii][jj]<<"\n" ;
        //     }
            
        // }
        

        for (int ii = 0; ii < mockupNum; ii++)
        {
            var_m2[ii] = (kmSum4[ii]/sampleNum - pow(mu_m2[ii][ii],2))*sampleNum/(sampleNum-1) ;
        }

        double invM[mockupNum][mockupNum] ;
        double buffer; //一時的なデータを蓄える
        
        //単位行列を作る
        for (int ii = 0; ii < mockupNum; ii++){
            for (int jj = 0; jj < mockupNum; jj++){
                invM[ii][jj] = (ii == jj)?1.0:0.0;
            }
        }
        //掃き出し法
        for (int ii = 0; ii < mockupNum; ii++){
            buffer = 1/cov_mm[ii][ii];
            for(int jj = 0; jj < mockupNum; jj++){
                cov_mm[ii][jj] *= buffer;
                invM[ii][jj] *= buffer;
            }
            for(int jj = 0; jj < mockupNum; jj++){
                if(ii != jj){
                    buffer = cov_mm[jj][ii];
                    for(int kk = 0; kk < mockupNum; kk++){
                    cov_mm[jj][kk] -= cov_mm[ii][kk]*buffer;
                    invM[jj][kk] -= invM[ii][kk]*buffer ;
                    }
                }
            }
        }

        double km[sampleNum] ;
        
        for (int ii = 0; ii < targetNum; ii++)
        {

            double sum1_m = 0. ;
            double sum2_m = 0. ;
            double sum4_m = 0. ;

            double sum1_tm = 0. ;
            double sum2_tm = 0. ;

            double mu_t = ktSum[ii]/sampleNum ;
            double mu_t2 = ktSum2[ii]/sampleNum ;

            double var_t = (mu_t2-pow(mu_t,2))*sampleNum/(sampleNum-1) ;

            double var_t2 = (ktSum4[ii]/sampleNum-pow(mu_t2,2))*sampleNum/(sampleNum-1);

            double cov_tm[mockupNum] ;
            double cov_tm2[mockupNum] ;
            for (int jj = 0; jj < mockupNum; jj++)
            {
                cov_tm[jj] = (ktkmSum[ii][jj]/sampleNum - mu_t*mu_m[jj])*sampleNum/(sampleNum-1) ;
                cov_tm2[jj] = (ktkmSum2[ii][jj]/sampleNum - mu_t2*mu_m2[jj][jj])*sampleNum/(sampleNum-1) ;
            }

            double a[mockupNum] ;
            for (int jj = 0; jj < mockupNum; jj++)
            {
                a[jj] = 0. ;
                for (int kk = 0; kk < mockupNum; kk++)
                {
                    a[jj] += invM[jj][kk]*cov_tm[kk] ;
                }
                
            }

            #if Pu9
                a[0] = 1. ;
                a[1] = 0. ;
            #endif

            #if Pu2
                a[0] = 0. ;
                a[1] = 1. ;
            #endif

            // for (int jj = 0; jj < mockupNum; jj++)
            // {
            //     ofs<<a[jj]<<"\n";
            // }
            

            for (int jj = 0; jj < sampleNum; jj++)
            {

                
                
                
                km[jj] = 0. ;

                for (int kk = 0; kk < mockupNum; kk++)
                {
                    km[jj] += a[kk]*km_i[kk][jj] ;
                }
                
                sum1_m += km[jj] ;
                sum2_m += pow(km[jj],2) ;
                sum4_m += pow(km[jj],4) ;

                sum1_tm += kt[ii][jj]*km[jj] ;
                sum2_tm += pow(kt[ii][jj]*km[jj],2) ;

            }
            
            double mu_m_vi = sum1_m/sampleNum ;
            double mu_m2_vi = sum2_m/sampleNum ;

            double var_m_vi = (mu_m2_vi-pow(mu_m_vi,2))*sampleNum/(sampleNum-1) ;
            double var_m2_vi = (sum4_m/sampleNum - pow(mu_m2_vi,2))*sampleNum/(sampleNum-1) ;

            double cov_tm_vi = (sum1_tm/sampleNum - mu_t*mu_m_vi)*sampleNum/(sampleNum-1) ;
            double corr_tm_vi = cov_tm_vi/(sqrt(var_t)*sqrt(var_m_vi)) ;

            double cov_t2m2_vi = (sum2_tm/sampleNum-mu_t2*mu_m2_vi)*sampleNum/(sampleNum-1) ;
            double corr_t2m2_vi = cov_t2m2_vi/(sqrt(var_t2)*sqrt(var_m2_vi)) ;

            double alpha = cov_tm_vi/var_m_vi; 
            double beta = cov_t2m2_vi/var_m2_vi;

            ofs<<corr_tm_vi<<"\n";

            double sample_H[sampleNum] ;
            double sample_Hbar[sampleNum] ;

            for (int jj = 0; jj < sampleNum; jj++)
            {
                sample_H[jj] = kt[ii][jj] - alpha*km[jj] ;
                sample_Hbar[jj] = pow(kt[ii][jj],2) - beta*pow(km[jj],2) ;
            }
            
            double sum1_h = 0.;
            double sum2_h = 0.;
            double sum1_hbar = 0.;
            double sum2_hbar = 0.;

            for (int jj = 0; jj < sampleNum; jj++)
            {
                sum1_h += sample_H[jj];
                sum2_h += pow(sample_H[jj],2) ;

                sum1_hbar += sample_Hbar[jj] ;
                sum2_hbar += pow(sample_Hbar[jj],2) ;
            }
            
            double mu_h = sum1_h/sampleNum;
            double mu_hbar = sum1_hbar/sampleNum;

            // ofs<<mu_h<<" "<<mu_hbar<<"\n" ;

            double mu_t_est = mu_h ;
            double mu_m2_rigorous  = a[0]*a[0]*2.244672e-04+a[1]*a[1]*7.632769e-03+2*a[0]*a[1]*9.526347e-06;

            double mu_t2_est = mu_hbar+beta*mu_m2_rigorous;
            double var_t_est = mu_t2_est-pow(mu_t_est,2);
            if(var_t_est<0)var_t_est = 0.;

            mu_t_cas[ii][i] = mu_t ;
            mu_t_est_cas[ii][i] = mu_t_est ;
            std_t_cas[ii][i] = sqrt(var_t) ;
            std_t_est_cas[ii][i] = sqrt(var_t_est) ;
        }
    
    }

    for (int ii = 0; ii < targetNum; ii++)
    {
        double sum1_mut = 0.;
        double sum2_mut = 0.;
        double sum1_mutest = 0.;
        double sum2_mutest = 0.;
        double sum1_stdt = 0.;
        double sum2_stdt = 0.;
        double sum1_stdtest = 0.;
        double sum2_stdtest = 0.;
        for(int j = 0; j < cal; j++){
            sum1_mut += mu_t_cas[ii][j];
            sum2_mut += pow(mu_t_cas[ii][j],2);
            sum1_mutest += mu_t_est_cas[ii][j];
            sum2_mutest += pow(mu_t_est_cas[ii][j],2);
            sum1_stdt += std_t_cas[ii][j];
            sum2_stdt += pow(std_t_cas[ii][j],2);
            sum1_stdtest += std_t_est_cas[ii][j];
            sum2_stdtest += pow(std_t_est_cas[ii][j],2);
        }

        sum1_mut /= cal;
        sum1_mutest /= cal;
        sum1_stdt /= cal;
        sum1_stdtest /= cal;

        double var_mut = (sum2_mut/cal-pow(sum1_mut,2))*cal/(cal-1);
        double var_mutest = (sum2_mutest/cal-pow(sum1_mutest,2))*cal/(cal-1);
        double var_stdt = (sum2_stdt/cal-pow(sum1_stdt,2))*cal/(cal-1);
        double var_stdtest = (sum2_stdtest/cal-pow(sum1_stdtest,2))*cal/(cal-1);
        if(var_mutest < 0.)var_mutest = 0.;
        if(var_stdtest < 0.)var_stdtest = 0.;    

        double ur_mu = sqrt(var_mutest)/sqrt(var_mut);    
        double ur_std = sqrt(var_stdtest)/sqrt(var_stdt); 

        #if mean
        ofs<<target[ii]<<" "<<sum1_mut<<"\n";
        #endif

        #if deviation
        ofs<<target[ii]<<" "<<ur_std<<"\n";
        #endif
    }
    

    return 0 ;
};