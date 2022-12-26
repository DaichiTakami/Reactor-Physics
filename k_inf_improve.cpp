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

int random(int fileNum){
    return rand()%fileNum ;
}

vector<string> split(string& input, char delimiter)
{
    istringstream stream(input);
    string field;
    vector<string> result;
    while (getline(stream, field, delimiter)) {
        result.push_back(field);
    }
    return result;
}

int main(){

    srand(10) ;

    ifstream ifs ;

    char output[256] ;
    cout<<"Enter a name for the output file."<<"\n" ;
    cin>>output;

    const char *outputfile = output ;
    ofstream ofs(outputfile) ;

    int fileNum = 0 ;

    int sampleNum ;
    cout<<"Enter the number of samples."<<"\n" ;
    cin>>sampleNum ;

    int num[sampleNum] ;

    int cal = 100 ;

    int loopMax = 1000000000 ;

    char filename[256] ;
    char buf[256] ;

    int targetNum = 0 ;

    ifs.open("SAMPLE/cal_point");

    while (ifs.getline(buf,sizeof(buf)))
    {
        targetNum ++ ;
    }
    
    // ofs<<"line number = "<<linenum<<"\n" ;

    ifs.clear() ;
    ifs.seekg(0, std::ios::beg) ;

    double *cal_point ;
    cal_point = new double[targetNum] ;

    for (int jj = 0; jj < targetNum; jj++)
    {
        ifs.getline(buf,sizeof(buf)) ;
        cal_point[jj] = atof(buf) ;
    }

    ifs.close() ;

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

        // cout<<ref[jj]<<"\n" ;
    }

    ifs.close() ;


    int mockupNum = linenum - targetNum ;

    
    double ktSum ;
    vector<double> km_iSum(mockupNum) ;

    double ktSum2 ;
    double km_iSum2[mockupNum][mockupNum] ;

    double ktSum4 ;
    vector<double> km_iSum4(mockupNum) ;

    double ktkm_iSum[mockupNum] ;
    double ktkm_iSum2[mockupNum] ;

    double mu_t_cas[cal] ;
    double mu_t_est_cas[cal] ;
    double std_t_cas[cal] ;
    double std_t_est_cas[cal] ;

    for (int ii = 0; ; ii++)
    {
        
        sprintf(filename, "SAMPLE/sample_%d", ii);

        ifs.open(filename) ;
        if (ifs.fail())
        {
            break ;
        }

        ifs.close() ;

        fileNum ++ ;
    }

    cout<<"Number of files : "<<fileNum<<"\n" ;

    double kt[targetNum][fileNum] ;
    double km_i[mockupNum][fileNum] ;

    for (int ii = 0; ii < fileNum; ii++)
    {
        
        sprintf(filename, "SAMPLE/sample_%d", ii);

        ifs.open(filename) ;
        if (ifs.fail())
        {
            break ;
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

        // cout<<km_i[0][ii]<<" "<<km_i[1][ii]<<"\n" ;
    }
    
    for (int i = 0; i < targetNum; i++)
    {
        for (int j = 0; j < cal; j++)
        {
            ktSum = 0. ;
            ktSum2 = 0. ;
            ktSum4 = 0. ;

            for (int kk = 0; kk < mockupNum; kk++)
            {
                km_iSum[kk] = 0. ;
                km_iSum4[kk] = 0. ;
                for (int ii = 0; ii < mockupNum; ii++)
                {
                    km_iSum2[kk][ii] = 0. ;
                }
                ktkm_iSum[kk] = 0. ;
                ktkm_iSum2 [kk] = 0. ;
            }

            for (int ii = 0; ii < sampleNum; ii++)
            {
                num[ii] = random(fileNum) ;
                
                
                ktSum += kt[i][num[ii]] ;
                ktSum2 += pow(kt[i][num[ii]],2) ;
                ktSum4 += pow(kt[i][num[ii]],4) ;
                
                
                for (int jj = 0; jj < mockupNum; jj++)
                {
                    km_iSum[jj] += km_i[jj][num[ii]] ;
                    km_iSum4[jj] += pow(km_i[jj][num[ii]],4) ;

                    for (int kk = 0; kk < mockupNum; kk++)
                    {
                        km_iSum2[jj][kk] += km_i[jj][num[ii]]*km_i[kk][num[ii]] ;
                    }
                }

                
                for (int kk = 0; kk < mockupNum; kk++)
                {
                    ktkm_iSum[kk] += kt[i][num[ii]]*km_i[i][num[ii]] ;
                    ktkm_iSum2[kk] += pow(kt[i][num[ii]]*km_i[i][num[ii]],2) ;
                }   
                
            }

            double mu_m[mockupNum] ;
            double mu_m2[mockupNum][mockupNum] ;

            double var_m2[mockupNum] ;

            double cov_mm[mockupNum][mockupNum] ;

            for (int ii = 0; ii < mockupNum; ii++)
            {
                mu_m[ii] = km_iSum[ii]/sampleNum ;
            }

            for (int ii = 0; ii < mockupNum; ii++)
            {
                for (int jj = 0; jj < mockupNum; jj++)
                {
                    mu_m2[ii][jj] = km_iSum2[ii][jj]/sampleNum ;

                    cov_mm[ii][jj] = (mu_m2[ii][jj] - mu_m[ii]*mu_m[jj])*sampleNum/(sampleNum-1) ;
                }
            }
            
            for (int ii = 0; ii < mockupNum; ii++)
            {
                var_m2[ii] = (km_iSum4[ii]/sampleNum - pow(mu_m2[ii][ii],2))*sampleNum/(sampleNum-1) ;
            }

            double invM[mockupNum][mockupNum] ;
            double buffer; //一時的なデータを蓄える
            
            //単位行列を作る
            for (int ii = 0; ii < mockupNum; ii++)
            {
                for (int jj = 0; jj < mockupNum; jj++){
                    invM[ii][jj] = (ii == jj)?1.0:0.0;
                }
            }
            //掃き出し法
            for (int ii = 0; ii < mockupNum; ii++)
            {
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

            double sum1_m = 0. ;
            double sum2_m = 0. ;
            double sum4_m = 0. ;

            double sum1_tm = 0. ;
            double sum2_tm = 0. ;

            double mu_t = ktSum/sampleNum ;
            double mu_t2 = ktSum2/sampleNum ;

            double var_t = (mu_t2-pow(mu_t,2))*sampleNum/(sampleNum-1) ;

            double var_t2 = (ktSum4/sampleNum-pow(mu_t2,2))*sampleNum/(sampleNum-1);

            double cov_tm[mockupNum] ;
            double cov_tm2[mockupNum] ;
            for (int jj = 0; jj < mockupNum; jj++)
            {
                cov_tm[jj] = (ktkm_iSum[jj]/sampleNum - mu_t*mu_m[jj])*sampleNum/(sampleNum-1) ;
                cov_tm2[jj] = (ktkm_iSum2[jj]/sampleNum - mu_t2*mu_m2[jj][jj])*sampleNum/(sampleNum-1) ;
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

            for (int jj = 0; jj < sampleNum; jj++)
            {
                km[num[jj]] = 0. ;

                for (int kk = 0; kk < mockupNum; kk++)
                {
                    km[jj] += a[kk]*km_i[kk][num[jj]] ;
                }
                
                sum1_m += km[jj] ;
                sum2_m += pow(km[jj],2) ;
                sum4_m += pow(km[jj],4) ;

                sum1_tm += kt[i][num[jj]]*km[jj] ;
                sum2_tm += pow(kt[i][num[jj]]*km[jj],2) ;

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

            // cout<<corr_tm_vi<<"\n" ;
 
            double sample_H[sampleNum] ;
            double sample_Hbar[sampleNum] ;

            for (int jj = 0; jj < sampleNum; jj++)
            {
                sample_H[jj] = kt[i][num[jj]] - alpha*km[jj] ;
                sample_Hbar[jj] = pow(kt[i][num[jj]],2) - beta*pow(km[jj],2) ;
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

            // cout<<mu_h<<" "<<mu_hbar<<"\n" ;

            double mu_t_est = mu_h ;
            double mu_m2_rigorous  = a[0]*a[0]*2.244672e-04+a[1]*a[1]*7.632769e-03+2*a[0]*a[1]*9.526347e-06;

            double mu_t2_est = mu_hbar+beta*mu_m2_rigorous;
            double var_t_est = mu_t2_est-pow(mu_t_est,2);
            if(var_t_est<0)var_t_est = 0.;

            mu_t_cas[j] = mu_t ;
            mu_t_est_cas[j] = mu_t_est ;
            std_t_cas[j] = sqrt(var_t) ;
            std_t_est_cas[j] = sqrt(var_t_est) ;
            
        
        }

        
        double sum1_mut = 0.;
        double sum2_mut = 0.;
        double sum1_mutest = 0.;
        double sum2_mutest = 0.;
        double sum1_stdt = 0.;
        double sum2_stdt = 0.;
        double sum1_stdtest = 0.;
        double sum2_stdtest = 0.;
        for(int j = 0; j < cal; j++){
            sum1_mut += mu_t_cas[j];
            sum2_mut += pow(mu_t_cas[j],2);
            sum1_mutest += mu_t_est_cas[j];
            sum2_mutest += pow(mu_t_est_cas[j],2);
            sum1_stdt += std_t_cas[j];
            sum2_stdt += pow(std_t_cas[j],2);
            sum1_stdtest += std_t_est_cas[j];
            sum2_stdtest += pow(std_t_est_cas[j],2);
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

        ofs<<cal_point[i]<<" "<<ur_mu<<"\n" ;
    }
    
        

        
    

    return 0 ;
};