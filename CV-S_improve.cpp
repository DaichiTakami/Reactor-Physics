#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>

#define conventional 0
#define dimensionless 1
#define mean 1
#define deviation 0
#define test 0
#define case 10

using namespace std ;

double gauss (double mu, double sigma)
{
    double x = (double)rand()/(double)RAND_MAX ;
    double y = (double)rand()/(double)RAND_MAX ;

    double z = sqrt(-2*log(x))*sin(2*M_PI*y) ;

    return mu + sigma*z ;
};

int main(){

    const char *filename = "dimensionless_mean_9.dat" ;
    ofstream ofs(filename) ;
    
    srand(10) ;

    int samples = 100 ;
    int sets = 1000 ;
    int loops = 100 ;

    for (int ii = 0; ii < loops+1; ii++)
    {
        #if case==1
            // Case 1
            double k_t_bar=ii*0.01;
            double k_m_bar=0.5;
            double s_t1=0.1;    
            double s_m1=0.1;
            double s_t2=0.;
            double target=k_t_bar;
        #endif    
            
        #if case==2
            // Case 2
            double k_m_bar=ii*0.01;
            double k_t_bar=0.5;
            double s_t1=0.1;    
            double s_m1=0.1;
            double s_t2=0.;
            double target=k_m_bar;
        #endif    
            
        #if case==3
            // Case 3
            double k_t_bar=ii*0.01;
            double k_m_bar=0.5;
            double s_t1=0.1;    
            double s_m1=0.2;
            double s_t2=0.;
            double target=k_t_bar;
        #endif    
            
        #if case==4
            // Case 4
            double k_t_bar=0.5;
            double k_m_bar=ii*0.01;    
            double s_t1=0.1;
            double s_m1=0.2;
            double s_t2=0.;
            double target=k_m_bar;
        #endif    
            
        #if case==5
            // Case 5
            double k_m_bar=0.5;
            double k_t_bar=ii*0.01;    
            double s_t1=0.1;
            double s_m1=0.1;
            double s_t2=0.01;
            double target=k_t_bar;
        #endif    
            
        #if case==6
            // Case 6
            double k_t_bar=0.5;
            double k_m_bar=ii*0.01;    
            double s_t1=0.1;
            double s_m1=0.1;
            double s_t2=0.01;
            double target=k_m_bar;
        #endif    
            
        #if case==7
            // Case 7
            double k_m_bar=0.5;
            double k_t_bar=ii*0.01;    
            double s_t1=0.1;
            double s_m1=0.2;
            double s_t2=0.01;
            double target=k_t_bar;
        #endif    
            
        #if case==8
            // Case 8
            double k_t_bar=0.5;
            double k_m_bar=ii*0.01;    
            double s_t1=0.1;
            double s_m1=0.2;
            double s_t2=0.01;
            double target=k_m_bar;
        #endif    
            
        #if case==9
            // Case 9
            double k_t_bar=0.1;
            double k_m_bar=0.2;
            double s_t1=ii*0.01;    
            double s_m1=0.4;
            double s_t2=0.01;
            double target=s_t1;
        #endif    
            
        #if case==10
            // Case 10
            double k_t_bar=0.1;
            double k_m_bar=0.2;
            double s_t2=ii*0.01;    
            double s_t1=0.3;
            double s_m1=0.4;
            double target=s_t2;
        #endif    
            
        #if case==11
            // Case 11
            double k_t_bar=0.1;
            double k_m_bar=0.2;
            double s_t1=0.3;
            double s_m1=0.3;
            double s_t2=0.01*ii;
            double target=s_t2;
        #endif    
            
        #if case==12
            // Case 12
            double k_t_bar=0.1;
            double k_m_bar=0.2;
            double s_t1=0.3;
            double s_t2=0.01;
            double s_m1=0.01*ii;
            double target=s_m1;
        #endif    

        double mean_inp = 0. ;
        double std_inp = 0.01 ;

        vector<double>  mu_t_cas(samples);
        vector<double>  mu_t_est_cas(samples);
        vector<double>  std_t_cas(samples);
        vector<double>  std_t_est_cas(samples);

        for (int jj = 0; jj < sets; jj++)
        {
            
        

            double sample[samples] ;
            for (int i = 0; i < samples; i++)
            {
                sample[i] = gauss(mean_inp,std_inp) ;
            }
            
            double sum1_t = 0. ;
            double sum2_t = 0. ;
            double sum4_t = 0. ;
            double sum1_m = 0. ;
            double sum2_m = 0. ;
            double sum4_m = 0. ;
            double sum1_tm = 0. ;
            double sum2_tm = 0. ;

            for (int i = 0; i < samples; i++)
            {
                double ds_s = sample[i] ;

                #if conventional
                    double k_t_i = k_t_bar*(1.+s_t1*ds_s+s_t2*ds_s*ds_s) ;
                    double k_m_i = k_m_bar*(1.+s_m1*ds_s) ;
                #endif

                #if dimensionless
                    double k_t_i = s_t1*ds_s+s_t2*ds_s*ds_s ;
                    double k_m_i = s_m1*ds_s ;
                #endif 

                sum1_t += k_t_i ;
                sum2_t += pow(k_t_i,2) ;
                sum4_t += pow(k_t_i,4) ;
                sum1_m += k_m_i ;
                sum2_m += pow(k_m_i,2);
                sum4_m += pow(k_m_i,4) ;
                sum1_tm += k_t_i*k_m_i ;
                sum2_tm += pow(k_t_i*k_m_i,2) ;
            }
            
            double mu_t = sum1_t/samples ;
            double mu_m = sum1_m/samples ;
            double mu_t2 = sum2_t/samples ;
            double mu_m2 = sum2_m/samples ;

            double var_t = (mu_t2-pow(mu_t,2))*samples/(samples-1) ;
            double var_m = (mu_m2-pow(mu_m,2))*samples/(samples-1) ;

            double var_t2 = (sum4_t/samples-pow(mu_t2,2))*samples/(samples-1);
            double var_m2 = (sum4_m/samples-pow(mu_m2,2))*samples/(samples-1);

            double cov_tm = (sum1_tm/samples-mu_t*mu_m)*samples/(samples-1) ;
            double corr_tm = cov_tm/(sqrt(var_t)*sqrt(var_m)) ;

            double cov_t2m2 = (sum2_tm/samples-mu_t2*mu_m2)*samples/(samples-1) ;
            double corr_t2m2 = cov_t2m2/(sqrt(var_t2)*sqrt(var_m2)) ;

            double alpha = cov_tm/var_m; 
            double beta = cov_t2m2/var_m2;

            double sample_H[samples] ;
            double sample_Hbar[samples] ;

            for (int i = 0; i < samples; i++)
            {
                double ds_s = sample[i] ;

                #if conventional
                    double k_t_i = k_t_bar*(1.+s_t1*ds_s+s_t2*ds_s*ds_s) ;
                    double k_m_i = k_m_bar*(1.+s_m1*ds_s) ;
                #endif

                #if dimensionless
                    double k_t_i = s_t1*ds_s+s_t2*ds_s*ds_s ;
                    double k_m_i = s_m1*ds_s ;
                #endif 

                sample_H[i] = k_t_i-alpha*k_m_i;
                sample_Hbar[i] = pow(k_t_i,2)-beta*pow(k_m_i,2) ;
            }

            double sum1_h = 0.;
            double sum2_h = 0.;
            double sum1_hbar = 0.;
            double sum2_hbar = 0.;

            for (int i = 0; i < samples; i++)
            {
                sum1_h += sample_H[i];
                sum2_h += pow(sample_H[i],2) ;

                sum1_hbar += sample_Hbar[i] ;
                sum2_hbar += pow(sample_Hbar[i],2) ;
            }
            
            double mu_h = sum1_h/samples;
            double mu_hbar = sum1_hbar/samples;

            #if conventional
                double mu_t_est = mu_h+alpha*k_m_bar;         
                double mu_m2_rigorous = k_m_bar*k_m_bar*s_m1*s_m1*std_inp*std_inp+k_m_bar*k_m_bar;
            #endif      

            #if dimensionless      
                double mu_t_est = mu_h;         
                double mu_m2_rigorous = s_m1*s_m1*std_inp*std_inp;
            #endif      
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
        for(int j = 0; j < sets; j++){
            sum1_mut += mu_t_cas[j];
            sum2_mut += pow(mu_t_cas[j],2);
            sum1_mutest += mu_t_est_cas[j];
            sum2_mutest += pow(mu_t_est_cas[j],2);
            sum1_stdt += std_t_cas[j];
            sum2_stdt += pow(std_t_cas[j],2);
            sum1_stdtest += std_t_est_cas[j];
            sum2_stdtest += pow(std_t_est_cas[j],2);
        };

        sum1_mut /= sets;
        sum1_mutest /= sets;
        sum1_stdt /= sets;
        sum1_stdtest /= sets;

        double var_mut = (sum2_mut/sets-pow(sum1_mut,2))*sets/(sets-1);
        double var_mutest = (sum2_mutest/sets-pow(sum1_mutest,2))*sets/(sets-1);
        double var_stdt = (sum2_stdt/sets-pow(sum1_stdt,2))*sets/(sets-1);
        double var_stdtest = (sum2_stdtest/sets-pow(sum1_stdtest,2))*sets/(sets-1);
        if(var_mutest < 0.)var_mutest = 0.;
        if(var_stdtest < 0.)var_stdtest = 0.;    

        double ur_mu = sqrt(var_mutest)/sqrt(var_mut);    
        double ur_std = sqrt(var_stdtest)/sqrt(var_stdt); 

        #if mean
        ofs<<target<<" "<<ur_mu<<"\n";
        #endif

        #if deviation
        ofs<<target<<" "<<ur_std<<"\n";
        #endif

        #if test
        ofs<<target<<" "<<sum1_mut<<" "<<sum1_mutest<<" "<<sum1_stdt<<" "<<sum1_stdtest<<"\n" ;
        #endif
    }
    
    
    return 0 ;
};