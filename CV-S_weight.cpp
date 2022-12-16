#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>

#define conventional 1
#define dimensionless 0
#define mean 1
#define deviation 0
#define test 0
#define case 1

using namespace std ;

double gauss (double mu, double sigma)
{
    double x = (double)rand()/(double)RAND_MAX ;
    double y = (double)rand()/(double)RAND_MAX ;

    double z = sqrt(-2*log(x))*sin(2*M_PI*y) ;

    return mu + sigma*z ;
}

int main(){
    
    const char *filename = "cc.dat" ;
    ofstream ofs(filename) ;
    
    srand(10) ;

    // int target = 100 ;
    int set = 100 ;
    int loop = 31 ;
    int sample[] = {10,15,20,25,30,40,50,60,70,90,100,150,200,250,300,400,500,600,700,900,1000,1500,2000,2500,3000,4000,5000,6000,7000,9000,10000} ;

    for (int ii = 0; ii < loop; ii++)
    {
        #if case == 1
            // case 1
            double kt_bar = 1. ;
            double km1_bar = 1. ;
            double km2_bar = 1. ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.0 ;
            double stTilde_2 = 0.0 ;
            double sm1Hat_1 = 0.5 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            int target = sample[ii] ;
        #endif

        double mean_inp_1 = 0. ;
        double mean_inp_2 = 0. ;
        double std_inp_1 = 0.01 ;
        double std_inp_2 = 0.01 ;

        vector<double>  mu_t_cas(target);
        vector<double>  mu_t_est_cas(target);
        vector<double>  std_t_cas(target);
        vector<double>  std_t_est_cas(target);

        double sum_a1 = 0. ;
        double sum_a2 = 0. ;

        for (int jj = 0; jj < set; jj++)
        {
            

            double target_1[target] ;
            for (int i = 0; i < target; i++)
            {
                target_1[i] = gauss(mean_inp_1,std_inp_1) ;
            }
            
            double target_2[target] ;
            for (int i = 0; i < target; i++)
            {
                target_2[i] = gauss(mean_inp_2,std_inp_2) ;
            }
            
            double sum1_t = 0. ;
            double sum2_t = 0. ;
            double sum4_t = 0. ;
            double sum1_m = 0. ;
            double sum2_m = 0. ;
            double sum4_m = 0. ;
            double sum1_m_1 = 0. ;
            double sum2_m_1 = 0. ;
            double sum4_m_1 = 0. ;
            double sum1_m_2 = 0. ;
            double sum2_m_2 = 0. ;
            double sum4_m_2 = 0. ;
            double sum1_tm = 0. ;
            double sum2_tm = 0. ;
            double sum1_tm_1 = 0. ;
            double sum2_tm_1 = 0. ;
            double sum1_tm_2 = 0. ;
            double sum2_tm_2 = 0. ;
            double sum1_mm = 0 ;
            

            for (int i = 0; i < target; i++)
            {
                double ds_s_1 = target_1[i] ;
                double ds_s_2 = target_2[i] ;

                double kt_i = kt_bar*(1.+stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2) ;
                double km1_i = km1_bar*(1.+sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2) ;
                double km2_i = km2_bar*(1.+sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2) ;
                
                

                #if dimensionless
                    kt_i /= kt_bar ;
                    kt_i -= 1. ;
                    km1_i /= km1_bar ;
                    km1_i -= 1. ;
                    km2_i /= km2_bar ;
                    km2_i -= 1. ;
                #endif

                // ofs<<km_i<<" "<<kt_i<<"\n" ;

                sum1_t += kt_i ;
                sum2_t += pow(kt_i,2) ;
                sum4_t += pow(kt_i,4) ;
                sum1_m_1 += km1_i ;
                sum2_m_1 += pow(km1_i,2);
                sum4_m_1 += pow(km1_i,4) ;
                sum1_m_2 += km2_i ;
                sum2_m_2 += pow(km2_i,2);
                sum4_m_2 += pow(km2_i,4) ;
                sum1_tm_1 += kt_i*km1_i ;
                sum2_tm_1 += pow(kt_i*km1_i,2) ;
                sum1_tm_2 += kt_i*km2_i ;
                sum2_tm_2 += pow(kt_i*km2_i,2) ;
                sum1_mm += km1_i*km2_i ;
            }   
            
            double mu_t = sum1_t/target ;
            double mu_m_1 = sum1_m_1/target ;
            double mu_m_2 = sum1_m_2/target ;
            double mu_t2 = sum2_t/target ;
            double mu_m2_1 = sum2_m_1/target ;
            double mu_m2_2 = sum2_m_2/target ;

            double var_t = (mu_t2-pow(mu_t,2))*target/(target-1) ;
            double var_m_1 = (mu_m2_1-pow(mu_m_1,2))*target/(target-1) ;
            double var_m_2 = (mu_m2_2-pow(mu_m_2,2))*target/(target-1) ;

            double var_t2 = (sum4_t/target-pow(mu_t2,2))*target/(target-1);
            double var_m2_1 = (sum4_m_1/target-pow(mu_m2_1,2))*target/(target-1);
            double var_m2_2 = (sum4_m_2/target-pow(mu_m2_2,2))*target/(target-1);

            double cov_tm_1 = (sum1_tm_1/target-mu_t*mu_m_1)*target/(target-1) ;
            double cov_tm_2 = (sum1_tm_2/target-mu_t*mu_m_2)*target/(target-1) ;

            double cov_mm = (sum1_mm/target-mu_m_1*mu_m_2)*target/(target-1) ;

            double cov_t2m2_1 = (sum2_tm_1/target-mu_t2*mu_m2_1)*target/(target-1) ;
            double cov_t2m2_2 = (sum2_tm_2/target-mu_t2*mu_m2_2)*target/(target-1) ;

            double M[2][2] ;
            M[0][0] = var_m_1 ;
            M[1][1] = var_m_2 ;
            M[0][1] = cov_mm ;
            M[1][0] = cov_mm ;

            double invM[2][2] ;
            double buf; //一時的なデータを蓄える
            int n=2;  //配列の次数
            
            //単位行列を作る
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    invM[i][j]=(i==j)?1.0:0.0;
                }
            }
            //掃き出し法
            for(int i=0;i<n;i++){
                buf=1/M[i][i];
                for(int j=0;j<n;j++){
                    M[i][j]*=buf;
                    invM[i][j]*=buf;
                }
                for(int j=0;j<n;j++){
                    if(i!=j){
                        buf=M[j][i];
                        for(int k=0;k<n;k++){
                        M[j][k]-=M[i][k]*buf;
                        invM[j][k]-=invM[i][k]*buf;
                        }
                    }
                }
            }

            double B[2] ;
            B[0] = cov_tm_1 ;
            B[1] = cov_tm_2 ;

            double a[2] ;
            for (int i = 0; i < n; i++)
            {
                a[i] = 0. ;
                for (int j = 0; j < n; j++)
                {
                    a[i] += invM[i][j]*B[j] ;
                }
                
            }
            
            double a_0 = a[0] ;
            double a_1 = a[1] ;

            sum_a1 += a_0 ;
            sum_a2 += a_1 ;

            // for (int i = 0; i < target; i++)
            // {
            //     double ds_s_1 = target_1[i] ;
            //     double ds_s_2 = target_2[i] ;

            //     double kt_i = kt_bar*(1.+stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2) ;
            //     double km1_i = km1_bar*(1.+sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2) ;
            //     double km2_i = km2_bar*(1.+sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2) ;
            //     double km_i = pow(km1_i,a_0)*pow(km2_i,a_1) ;

            //     double km_bar = pow(km1_bar,a_0)*pow(km2_bar,a_1) ;

            //     #if dimensionless
            //         kt_i /= kt_bar ;
            //         kt_i -= 1. ;
            //         km1_i /= km1_bar ;
            //         km1_i -= 1. ;
            //         km2_i /= km2_bar ;
            //         km2_i -= 1. ;
            //         km_i /= km_bar ;
            //         km_i -= 1. ;
            //     #endif

            //     sum1_m += km_i ;
            //     sum2_m += pow(km_i,2);
            //     sum4_m += pow(km_i,4) ;
            //     sum1_tm += kt_i*km_i ;
            //     sum2_tm += pow(kt_i*km_i,2) ;
            // }
            
            // double mu_m = sum1_m/target ;
            // double mu_m2 = sum2_m/target ;

            // double var_m = (mu_m2-pow(mu_m,2))*target/(target-1) ;

            // double var_m2 = (sum4_m/target-pow(mu_m2,2))*target/(target-1);

            // double cov_tm = (sum1_tm/target-mu_t*mu_m)*target/(target-1) ;
            // double cov_t2m2 = (sum2_tm/target-mu_t2*mu_m2)*target/(target-1) ;

            // double alpha = cov_tm/var_m; 
            // double beta = cov_t2m2/var_m2;

            // // ofs<<target<<" "<<beta<<"\n" ;

            // double target_H[target] ;
            // double target_Hbar[target] ;

            // for (int i = 0; i < target; i++)
            // {
            //     double ds_s_1 = target_1[i] ;
            //     double ds_s_2 = target_2[i] ;

            //     double kt_i = kt_bar*(1.+stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2) ;
            //     double km1_i = km1_bar*(1.+sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2) ;
            //     double km2_i = km2_bar*(1.+sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2) ;
            //     double km_i = pow(km1_i,a_0)*pow(km2_i,a_1) ;

            //     double km_bar = pow(km1_bar,a_0)*pow(km2_bar,a_1) ;

            //     #if dimensionless
            //         kt_i /= kt_bar ;
            //         kt_i -= 1. ;
            //         km1_i /= km1_bar ;
            //         km1_i -= 1. ;
            //         km2_i /= km2_bar ;
            //         km2_i -= 1. ;
            //         km_i /= km_bar ;
            //         km_i -= 1. ;
            //     #endif

            //     target_H[i] = kt_i-alpha*km_i;
            //     target_Hbar[i] = pow(kt_i,2)-beta*pow(km_i,2) ;

            //     // ofs<<target_H[i]<<"\n" ;
            
            // }

            // double sum1_h = 0.;
            // double sum2_h = 0.;
            // double sum1_hbar = 0.;
            // double sum2_hbar = 0.;

            // for (int i = 0; i < target; i++)
            // {
            //     sum1_h += target_H[i];
            //     sum2_h += pow(target_H[i],2) ;

            //     sum1_hbar += target_Hbar[i] ;
            //     sum2_hbar += pow(target_Hbar[i],2) ;
            // }
            
            // double mu_h = sum1_h/target;
            // double mu_hbar = sum1_hbar/target;

            // #if conventional
            //     double mu_t_est = mu_h+alpha*pow(km1_bar,a_0)*pow(km2_bar,a_1);         
            //     double mu_m2_1_rigorous = km1_bar*km1_bar*pow(sm1Hat_1*std_inp_1+sm1Hat_1*std_inp_2,2)+km1_bar*km1_bar ;
            //     double mu_m2_2_rigorous = km2_bar*km2_bar*pow(sm2Hat_1*std_inp_1+sm2Hat_2*std_inp_2,2)+km2_bar*km2_bar ;
            //     double mu_m2_rigorous = pow(mu_m2_1_rigorous,a_0)*pow(mu_m2_2_rigorous,a_1) ;
            // #endif     

            // #if dimensionless
            //     double mu_t_est = mu_h ;         
            //     double mu_m2_1_rigorous = pow(sm1Hat_1*std_inp_1+sm1Hat_2*std_inp_2,2) ;
            //     double mu_m2_2_rigorous = pow(sm2Hat_1*std_inp_1+sm2Hat_2*std_inp_2,2) ;
            //     double mu_m2_rigorous = pow(mu_m2_1_rigorous,a_0)*pow(mu_m2_2_rigorous,a_1) ;
            // #endif     

            // // ofs<<mu_t_est<<" "<<mu_m2_1_rigorous<<" "<<mu_m2_2_rigorous<<"\n" ;

            // double mu_t2_est = mu_hbar+beta*mu_m2_rigorous;
            // double var_t_est = mu_t2_est-pow(mu_t_est,2);
            // if(var_t_est<0)var_t_est = 0.;

            // mu_t_cas[jj] = mu_t ;
            // mu_t_est_cas[jj] = mu_t_est ;
            // std_t_cas[jj] = sqrt(var_t) ;
            // std_t_est_cas[jj] = sqrt(var_t_est) ;

            // // ofs<<mu_m2_rigorous<<" "<<mu_t2_est<<"\n" ;

        }

        double a0_est = sum_a1/set ;
        double a1_est = sum_a2/set ;

        ofs<<target<<" "<<a1_est<<"\n";
        
        // double sum1_mut = 0.;
        // double sum2_mut = 0.;
        // double sum1_mutest = 0.;
        // double sum2_mutest = 0.;
        // double sum1_stdt = 0.;
        // double sum2_stdt = 0.;
        // double sum1_stdtest = 0.;
        // double sum2_stdtest = 0.;
        // for(int j = 0; j < set; j++){
        //     sum1_mut += mu_t_cas[j];
        //     sum2_mut += pow(mu_t_cas[j],2);
        //     sum1_mutest += mu_t_est_cas[j];
        //     sum2_mutest += pow(mu_t_est_cas[j],2);
        //     sum1_stdt += std_t_cas[j];
        //     sum2_stdt += pow(std_t_cas[j],2);
        //     sum1_stdtest += std_t_est_cas[j];
        //     sum2_stdtest += pow(std_t_est_cas[j],2);
        // }

        // sum1_mut /= set;
        // sum1_mutest /= set;
        // sum1_stdt /= set;
        // sum1_stdtest /= set;

        // double var_mut = (sum2_mut/set-pow(sum1_mut,2))*set/(set-1);
        // double var_mutest = (sum2_mutest/set-pow(sum1_mutest,2))*set/(set-1);
        // double var_stdt = (sum2_stdt/set-pow(sum1_stdt,2))*set/(set-1);
        // double var_stdtest = (sum2_stdtest/set-pow(sum1_stdtest,2))*set/(set-1);
        // if(var_mutest < 0.)var_mutest = 0.;
        // if(var_stdtest < 0.)var_stdtest = 0.;    

        // double ur_mu = sqrt(var_mutest)/sqrt(var_mut);    
        // double ur_std = sqrt(var_stdtest)/sqrt(var_stdt); 

        // #if mean
        // ofs<<target<<" "<<ur_mu<<"\n";
        // #endif

        // #if deviation
        // ofs<<target<<" "<<ur_std<<"\n";
        // #endif

        // #if test
        // ofs<<target<<" "<<sum1_mut<<" "<<sum1_mutest<<" "<<sum1_stdt<<" "<<sum1_stdtest<<"\n" ;
        // #endif


    }
    
    
    
    
    return 0 ;
};