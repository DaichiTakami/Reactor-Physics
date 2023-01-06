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

    char output[256] ;
    cout<<"Please enter a name for the output file."<<"\n" ;
    cin>>output ;

    const char *outputfile = output ;
    ofstream ofs(outputfile) ;

    srand(10) ;

    int fileNum = 0 ;
    int sampleNum[] = {10,15,20,25,30,35,40,50,60,80,100,150,200,250,300,350,400,500,600,800,1000,1500,2000,2500,3000,3500,4000,5000,6000,8000,10000} ;
    int cal = 100 ;

    int loopMax = 100000000 ;

    char filename[256] ;
    char buf[256] ;

    int targetNum = sizeof(sampleNum)/sizeof(int) ;

    // ofs<<targetNum;

    int target = 0 ;

    ifs.open("SAMPLE/cal_point");

    while (ifs.getline(buf,sizeof(buf)))
    {
        target ++ ;
    }
    
    // ofs<<"line number = "<<linenum<<"\n" ;

    ifs.close() ;

    ifs.open("SAMPLE/sample_ref") ;
    
    int lineNum = 0 ;

    while (ifs.getline(buf,sizeof(buf)))
    {
        lineNum ++ ;
    }
    
    ifs.clear() ;
    ifs.seekg(0,std::ios::beg) ;

    double *ref ;
    ref = new double[lineNum] ;

    for (int kk = 0; kk < lineNum; kk++)
    {
        ifs.getline(buf,sizeof(buf)) ;
        ref[kk] = atof(buf) ;
    }

    ifs.close() ;

    int mockupNum = lineNum - target ;

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

        // ofs<<filename<<"\n" ;

        ifs.open(filename) ;
        if (ifs.fail())
        {
            break ;
        }

        ifs.close() ;

        fileNum ++ ;
    }
    
    double kt[target][fileNum] ;
    double km_i[mockupNum][fileNum] ;

    for (int ii = 0; ii < fileNum; ii++)
    {
        
        sprintf(filename, "SAMPLE/sample_%d", ii);

        // ofs<<filename<<"\n" ;

        ifs.open(filename) ;
        if (ifs.fail())
        {
            break ;
        }

        double *arr ;
        arr = new double[lineNum] ;

        for (int jj = 0; jj < lineNum; jj++)
        {
            ifs.getline(buf,sizeof(buf)) ;
            arr[jj] = atof(buf) ;
        }

        ifs.close() ;

        for (int jj = 0; jj < lineNum; jj++)
        {
            arr[jj] /= ref[jj] ;
            arr[jj] -= 1. ;
        }

        for (int jj = 0; jj < target; jj++)
        {
            kt[jj][ii] = arr[jj] ;
        }
        
        for (int jj = 0; jj < mockupNum; jj++)
        {
            km_i[jj][ii] = arr[jj+target] ;
        }


    }

    int calPoint ;
    cout<<"Enter calculation points."<<"\n" ;
    cin>>calPoint;
    calPoint -- ;

    int calData ;
    cout<<"Enter calculation information."<<"\n" ;
    cout<<"If 0 is entered, the data is calculated in the conventional method."<<"\n" ;
    cin>>calData;


    for (int ii = 0; ii < targetNum; ii++)
    {
        int randFig[sampleNum[ii]] ;

        for (int jj = 0; jj < cal; jj++)
        {
            ktSum = 0. ;
            ktSum2 = 0. ;
            ktSum4 = 0. ;
            
            for (int kk = 0; kk < mockupNum; kk++)
            {
                km_iSum[kk] = 0. ;
                km_iSum4[kk] = 0. ;
                for (int i = 0; i < mockupNum; i++)
                {
                    km_iSum2[kk][i] = 0. ;
                }
                ktkm_iSum[kk] = 0. ;
                ktkm_iSum2 [kk] = 0. ;
            }
            

            for (int kk = 0; kk < sampleNum[ii]; kk++)
            {
                randFig[kk] = random(fileNum) ;

                ktSum += kt[calPoint][randFig[kk]] ;
                ktSum2 += pow(kt[calPoint][randFig[kk]],2) ;
                ktSum4 += pow(kt[calPoint][randFig[kk]],4) ;

                for (int i = 0; i < mockupNum; i++)
                {
                    km_iSum[i] += km_i[i][randFig[kk]] ;
                    km_iSum4[i] += pow(km_i[i][randFig[kk]],2) ;

                    for (int j = 0; j < mockupNum; j++)
                    {
                        km_iSum2[i][j] += km_i[i][randFig[kk]]*km_i[j][randFig[kk]] ;
                    }

                    ktkm_iSum[i] += kt[calPoint][randFig[kk]]*km_i[i][randFig[kk]] ;
                    ktkm_iSum2[i] += pow(kt[calPoint][randFig[kk]]*km_i[i][randFig[kk]],2) ;
                    
                }


                
            }

            double mu_m[mockupNum] ;
            double mu_m2[mockupNum][mockupNum] ;

            double var_m2[mockupNum] ;

            double cov_mm[mockupNum][mockupNum] ;

            for (int kk = 0; kk < mockupNum; kk++)
            {
                mu_m[kk] = km_iSum[kk]/sampleNum[ii] ;
            }

            for (int kk = 0; kk < mockupNum; kk++)
            {
                for (int i = 0; i < mockupNum; i++)
                {
                    mu_m2[kk][i] = km_iSum2[kk][i]/sampleNum[ii] ;

                    cov_mm[kk][i] = (mu_m2[kk][i]-mu_m[kk]*mu_m[i])*sampleNum[ii]/(sampleNum[ii]-1) ;

                }            
            }
            
            for (int kk = 0; kk < mockupNum; kk++)
            {
                var_m2[kk] = (km_iSum4[kk]/sampleNum[ii] - pow(mu_m2[kk][kk],2))*sampleNum[ii]/(sampleNum[ii]-1) ;
            }
            
            double invM[mockupNum][mockupNum] ;
            double buffer; //一時的なデータを蓄える
            
            //単位行列を作る
            for (int kk = 0; kk < mockupNum; kk++){
                for (int i = 0; i < mockupNum; i++){
                    invM[kk][i] = (kk == i)?1.0:0.0;
                }
            }
            //掃き出し法
            for (int kk = 0; kk < mockupNum; kk++){
                buffer = 1/cov_mm[kk][kk];
                for(int i = 0; i < mockupNum; i++){
                    cov_mm[kk][i] *= buffer;
                    invM[kk][i] *= buffer;
                }
                for(int i = 0; i < mockupNum; i++){
                    if(kk != i){
                        buffer = cov_mm[i][kk];
                        for(int kk = 0; kk < mockupNum; kk++){
                        cov_mm[i][kk] -= cov_mm[kk][kk]*buffer;
                        invM[i][kk] -= invM[kk][kk]*buffer ;
                        }
                    }
                }
            }

            double km[sampleNum[ii]] ;
            
            

            double sum1_m = 0. ;
            double sum2_m = 0. ;
            double sum4_m = 0. ;

            double sum1_tm = 0. ;
            double sum2_tm = 0. ;

            double mu_t = ktSum/sampleNum[ii] ;
            double mu_t2 = ktSum2/sampleNum[ii] ;

            double var_t = (mu_t2-pow(mu_t,2))*sampleNum[ii]/(sampleNum[ii]-1) ;

            double var_t2 = (ktSum4/sampleNum[ii]-pow(mu_t2,2))*sampleNum[ii]/(sampleNum[ii]-1);

            double cov_tm[mockupNum] ;
            double cov_tm2[mockupNum] ;
            for (int i = 0; i < mockupNum; i++)
            {
                cov_tm[i] = (ktkm_iSum[i]/sampleNum[ii] - mu_t*mu_m[i])*sampleNum[ii]/(sampleNum[ii]-1) ;
                cov_tm2[i] = (ktkm_iSum2[i]/sampleNum[ii] - mu_t2*mu_m2[i][i])*sampleNum[ii]/(sampleNum[ii]-1) ;
            }

            double a[mockupNum] ;
            for (int i = 0; i < mockupNum; i++)
            {
                a[i] = 0. ;
                for (int j = 0; j < mockupNum; j++)
                {
                    a[i] += invM[i][j]*cov_tm[j] ;
                }
                
            }

            if (calData == 9)
            {
                a[0] = 1. ;
                a[1] = 0. ;
            }

           if (calData == 2)
           {
                a[0] = 0. ;
                a[1] = 1. ;
           }
           

            for (int i = 0; i < sampleNum[ii]; i++)
            {                       
                km[i] = 0. ;

                for (int j = 0; j < mockupNum; j++)
                {
                    km[i] += a[j]*km_i[j][randFig[i]] ;
                }
                
                sum1_m += km[i] ;
                sum2_m += pow(km[i],2) ;
                sum4_m += pow(km[i],4) ;

                sum1_tm += kt[calPoint][randFig[i]]*km[i] ;
                sum2_tm += pow(kt[calPoint][randFig[i]]*km[i],2) ;

            }
            
            double mu_m_vi = sum1_m/sampleNum[ii] ;
            double mu_m2_vi = sum2_m/sampleNum[ii] ;

            double var_m_vi = (mu_m2_vi-pow(mu_m_vi,2))*sampleNum[ii]/(sampleNum[ii]-1) ;
            double var_m2_vi = (sum4_m/sampleNum[ii] - pow(mu_m2_vi,2))*sampleNum[ii]/(sampleNum[ii]-1) ;

            double cov_tm_vi = (sum1_tm/sampleNum[ii] - mu_t*mu_m_vi)*sampleNum[ii]/(sampleNum[ii]-1) ;
            double corr_tm_vi = cov_tm_vi/(sqrt(var_t)*sqrt(var_m_vi)) ;

            double cov_t2m2_vi = (sum2_tm/sampleNum[ii]-mu_t2*mu_m2_vi)*sampleNum[ii]/(sampleNum[ii]-1) ;
            double corr_t2m2_vi = cov_t2m2_vi/(sqrt(var_t2)*sqrt(var_m2_vi)) ;

            cout<<corr_tm_vi<<" "<<corr_t2m2_vi<<"\n" ;

            double alpha = cov_tm_vi/var_m_vi; 
            double beta = cov_t2m2_vi/var_m2_vi;

            double sample_H[sampleNum[ii]] ;
            double sample_Hbar[sampleNum[ii]] ;

            for (int i = 0; i < sampleNum[ii]; i++)
            {
                sample_H[i] = kt[calPoint][randFig[i]] - alpha*km[i] ;
                sample_Hbar[i] = pow(kt[calPoint][randFig[i]],2) - beta*pow(km[i],2) ;
            }
            
            double sum1_h = 0.;
            double sum2_h = 0.;
            double sum1_hbar = 0.;
            double sum2_hbar = 0.;

            for (int i = 0; i < sampleNum[ii]; i++)
            {
                sum1_h += sample_H[i];
                sum2_h += pow(sample_H[i],2) ;

                sum1_hbar += sample_Hbar[i] ;
                sum2_hbar += pow(sample_Hbar[i],2) ;
            }
            
            double mu_h = sum1_h/sampleNum[ii];
            double mu_hbar = sum1_hbar/sampleNum[ii];

            // ofs<<mu_h<<" "<<mu_hbar<<"\n" ;

            double mu_t_est = mu_h ;
            double mu_m2_rigorous  = a[0]*a[0]*2.244672e-04+a[1]*a[1]*7.632769e-03+2*a[0]*a[1]*9.526347e-06;

            double mu_t2_est = mu_hbar+beta*mu_m2_rigorous;
            double var_t_est = mu_t2_est-pow(mu_t_est,2);
            if(var_t_est<0)var_t_est = 0.;

            mu_t_cas[jj] = mu_t ;
            mu_t_est_cas[jj] = mu_t_est ;
            std_t_cas[jj] = sqrt(var_t) ;
            std_t_est_cas[jj] = sqrt(var_t_est) ;
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

        if (calData == 0)
        {
            ofs<<sampleNum[ii]<<" "<<sum1_stdt<<" "<<sqrt(var_stdt)<<"\n" ;
        }else{
            ofs<<sampleNum[ii]<<" "<<sum1_stdtest<<" "<<sqrt(var_stdtest)<<"\n" ;
        }
        
    


            
        
        

    }
    

}