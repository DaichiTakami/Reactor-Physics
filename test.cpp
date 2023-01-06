#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>

<<<<<<< HEAD
#define conventional 0
#define dimensionless 1
#define mean 0
#define deviation 0
#define test 1
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
    
    const char *filename = "aa.dat" ;
    ofstream ofs(filename) ;
    
    srand(10) ;

    const int sample = 100 ;
    int set = 100 ;
    int loop = 101 ;

    for (int ii = 0; ii < loop; ii++)
    {
        #if case == 1
            // case 1
            double kt_bar = 1. ;
            double km1_bar = 1. ;
            double km2_bar = 1. ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0. ;
            double stTilde_2 = 0. ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = ii+1 ;
        #endif

        #if case == 2
            // case 2
            double kt_bar = 0.5 ;
            double km1_bar = 0.01*ii ;
            double km2_bar = 0.5 ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0. ;
            double stTilde_2 = 0. ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = km1_bar ;
        #endif

        #if case == 3
            // case 3
            double kt_bar = 0.01*ii ;
            double km1_bar = 0.5 ;
            double km2_bar = 0.5 ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01 ;
            double stTilde_2 = 0. ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = kt_bar ;
        #endif

        #if case == 4
            // case 4
            double kt_bar = 0.5 ;
            double km1_bar = 0.01*ii ;
            double km2_bar = 0.5 ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01 ;
            double stTilde_2 = 0. ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = km1_bar ;
        #endif

        #if case == 5
            // case 5
            double kt_bar = 0.01*ii ;
            double km1_bar = 0.5 ;
            double km2_bar = 0.5 ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01 ;
            double stTilde_2 = 0.01 ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = kt_bar ;
        #endif

        #if case == 6
            // case 6
            double kt_bar = 0.5 ;
            double km1_bar = 0.01*ii ;
            double km2_bar = 0.5 ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01 ;
            double stTilde_2 = 0.01 ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = km1_bar ;
        #endif

        #if case == 7
            // case 7
            double kt_bar = 1. ;
            double km1_bar = 1. ;
            double km2_bar = 1. ;
            double stHat_1 = 0.01*ii ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01 ;
            double stTilde_2 = 0.01 ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = stHat_1 ;
        #endif

        #if case == 8
            // case 8
            double kt_bar = 1. ;
            double km1_bar = 1. ;
            double km2_bar = 1. ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01*ii ;
            double stTilde_2 = 0.01 ;
            double sm1Hat_1 = 0.1 ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = stTilde_1 ;
        #endif

        #if case == 9
            // case 9
            double kt_bar = 1. ;
            double km1_bar = 1. ;
            double km2_bar = 1. ;
            double stHat_1 = 0.1 ;
            double stHat_2 = 0.1 ;
            double stTilde_1 = 0.01 ;
            double stTilde_2 = 0.01 ;
            double sm1Hat_1 = 0.01*ii ;
            double sm1Hat_2 = 0.2 ;
            double sm2Hat_1 = 0.2 ;
            double sm2Hat_2 = 0.1 ;
            double target = sm1Hat_1 ;
        #endif

        double mean_inp_1 = 0. ;
        double mean_inp_2 = 0. ;
        double std_inp_1 = 0.01 ;
        double std_inp_2 = 0.01 ;

        vector<double>  mu_t_cas(sample);
        vector<double>  mu_t_est_cas(sample);
        vector<double>  std_t_cas(sample);
        vector<double>  std_t_est_cas(sample);

        double sum_a0 = 0. ;
        double sum_a1 = 0. ;

        for (int jj = 0; jj < set; jj++)
        {
            double sample_1[sample] ;
            for (int i = 0; i < sample; i++)
            {
                sample_1[i] = gauss(mean_inp_1,std_inp_1) ;
            }
            
            double sample_2[sample] ;
            for (int i = 0; i < sample; i++)
            {
                sample_2[i] = gauss(mean_inp_2,std_inp_2) ;
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
            double sum1_mm = 0. ;

            for (int i = 0; i < sample; i++)
            {
                double ds_s_1 = sample_1[i] ;
                double ds_s_2 = sample_2[i] ;

                double kt_i = kt_bar*(1.+stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2) ;
                double km1_i = km1_bar*(1.+sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2) ;
                double km2_i = km2_bar*(1.+sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2) ;

                #if dimensionless
                    kt_i = stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2 ;
                    km1_i = sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2 ;
                    km2_i = sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2 ;
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
            
            double mu_t = sum1_t/sample ;
            double mu_m_1 = sum1_m_1/sample ;
            double mu_m_2 = sum1_m_2/sample ;
            double mu_t2 = sum2_t/sample ;
            double mu_m2_1 = sum2_m_1/sample ;
            double mu_m2_2 = sum2_m_2/sample ;

            double var_t = (mu_t2-pow(mu_t,2))*sample/(sample-1) ;
            double var_m_1 = (mu_m2_1-pow(mu_m_1,2))*sample/(sample-1) ;
            double var_m_2 = (mu_m2_2-pow(mu_m_2,2))*sample/(sample-1) ;

            double var_t2 = (sum4_t/sample-pow(mu_t2,2))*sample/(sample-1);
            double var_m2_1 = (sum4_m_1/sample-pow(mu_m2_1,2))*sample/(sample-1);
            double var_m2_2 = (sum4_m_2/sample-pow(mu_m2_2,2))*sample/(sample-1);

            double cov_tm_1 = (sum1_tm_1/sample-mu_t*mu_m_1)*sample/(sample-1) ;
            double cov_tm_2 = (sum1_tm_2/sample-mu_t*mu_m_2)*sample/(sample-1) ;

            double cov_t2m2_1 = (sum2_tm_1/sample-mu_t2*mu_m2_1)*sample/(sample-1) ;
            double cov_t2m2_2 = (sum2_tm_2/sample-mu_t2*mu_m2_2)*sample/(sample-1) ;

            double cov_mm = (sum1_mm/sample-mu_m_1*mu_m_2)*sample/(sample-1) ;

            double M[2][2] ;
            M[0][0] = var_m_1 ;
            M[1][1] = var_m_2 ;
            M[0][1] = cov_mm ;
            M[1][0] = cov_mm ;

            // ofs<<var_m?<<" "<<var_m_2<<" "<<cov_mm<<"\n";

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
=======
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
>>>>>>> 2b7bfc688d139f55b0773e421cf6a0444ffe9f84
                }
                
            }
            
<<<<<<< HEAD
            double a_0 = a[0] ;
            double a_1 = a[1] ;

            sum_a0 += a_0 ;
            sum_a1 += a_1 ;

            // ofs<<a_0<<" "<<a_1<<"\n";
            
            for (int i = 0; i < sample; i++)
            {
                double ds_s_1 = sample_1[i] ;
                double ds_s_2 = sample_2[i] ;

                double kt_i = kt_bar*(1.+stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2) ;
                double km1_i = km1_bar*(1.+sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2) ;
                double km2_i = km2_bar*(1.+sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2) ;
                double km_i = pow(km1_i,a_0)*pow(km2_i,a_1) ;

                #if dimensionless
                    kt_i = stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2 ;
                    km_i = (sm1Hat_1*a_0+sm1Hat_2*a_1)*ds_s_1+(sm2Hat_1*a_0+sm2Hat_2*a_1)*ds_s_2 ;
                #endif

                // ofs<<kt_i<<" "<<km_i<<"\n" ;

                sum1_m += km_i ;
                sum2_m += pow(km_i,2);
                sum4_m += pow(km_i,4) ;
                sum1_tm += kt_i*km_i ;
                sum2_tm += pow(kt_i*km_i,2) ;
                sum1_tm_1 += kt_i*km1_i ;
                sum2_tm_1 += pow(kt_i*km1_i,2) ;
                sum1_tm_2 += kt_i*km2_i ;
                sum2_tm_2 += pow(kt_i*km2_i,2) ;
            }   
            
            double mu_m = sum1_m/sample ;
            double mu_m2 = sum2_m/sample ;

            double var_m = (mu_m2-pow(mu_m,2))*sample/(sample-1) ;

            double var_m2 = (sum4_m/sample-pow(mu_m2,2))*sample/(sample-1);

            double cov_tm = (sum1_tm/sample-mu_t*mu_m)*sample/(sample-1) ;
            double corr_tm = cov_tm/(sqrt(var_t)*sqrt(var_m)) ;

            double cov_t2m2 = (sum2_tm/sample-mu_t2*mu_m2)*sample/(sample-1) ;
            double corr_t2m2 = cov_t2m2/(sqrt(var_t2)*sqrt(var_m2)) ;

            double alpha = cov_tm/var_m; 
            double beta = cov_t2m2/var_m2;

            // ofs<<alpha<<" "<<beta<<"\n" ;

            double sample_H[sample] ;
            double sample_Hbar[sample] ;

            for (int i = 0; i < sample; i++)
            {
                double ds_s_1 = sample_1[i] ;
                double ds_s_2 = sample_2[i] ;

                double kt_i = kt_bar*(1.+stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2) ;
                double km1_i = km1_bar*(1.+sm1Hat_1*ds_s_1+sm1Hat_2*ds_s_2) ;
                double km2_i = km2_bar*(1.+sm2Hat_1*ds_s_1+sm2Hat_2*ds_s_2) ;
                double km_i = pow(km1_i,a_0)*pow(km2_i,a_1) ;

                double km_bar = pow(km1_bar,a_0)*pow(km2_bar,a_1) ;

                #if dimensionless
                    kt_i = stHat_1*ds_s_1+stHat_2*ds_s_2+stTilde_1*ds_s_1*ds_s_1+stTilde_2*ds_s_2*ds_s_2 ;
                    km_i = (sm1Hat_1*a_0+sm1Hat_2*a_1)*ds_s_1+(sm2Hat_1*a_0+sm2Hat_2*a_1)*ds_s_2 ;
                #endif

                sample_H[i] = kt_i-alpha*km_i;
                sample_Hbar[i] = pow(kt_i,2)-beta*pow(km_i,2) ;

                // ofs<<sample_H[i]<<"\n" ;
                // ofs<<kt_i<<" "<<km_i<<"\n" ;
                // ofs<<km1_i<<" "<<km2_i<<"\n" ;
                // ofs<<a_0<<" "<<a_1<<"\n" ;

            
            }

=======
            

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
            
>>>>>>> 2b7bfc688d139f55b0773e421cf6a0444ffe9f84
            double sum1_h = 0.;
            double sum2_h = 0.;
            double sum1_hbar = 0.;
            double sum2_hbar = 0.;

<<<<<<< HEAD
            for (int i = 0; i < sample; i++)
            {
                sum1_h += sample_H[i];
                sum2_h += pow(sample_H[i],2) ;

                sum1_hbar += sample_Hbar[i] ;
                sum2_hbar += pow(sample_Hbar[i],2) ;
            }
            
            double mu_h = sum1_h/sample;
            double mu_hbar = sum1_hbar/sample;

            // ofs<<mu_h<<" "<<mu_hbar<<"\n" ;

            double var_m2_1_rigorous = pow((sm1Hat_1*a_0+sm1Hat_2*a_1)*std_inp_1,2) ;
            double var_m2_2_rigorous = pow((sm2Hat_1*a_0+sm2Hat_2*a_1)*std_inp_2,2) ;

            #if conventional
                double mu_t_est = mu_h+alpha*pow(km1_bar,a_0)*pow(km2_bar,a_1);         
                double mu_m2_1_rigorous = km1_bar*km1_bar*var_m2_1_rigorous+km1_bar*km1_bar ;
                double mu_m2_2_rigorous = km2_bar*km2_bar*var_m2_2_rigorous+km2_bar*km2_bar ;
                double mu_m2_rigorous = pow(mu_m2_1_rigorous,a_0)*pow(mu_m2_2_rigorous,a_1) ;
            #endif     

            #if dimensionless
                double mu_t_est = mu_h ;         
                double mu_m2_rigorous = var_m2_1_rigorous+var_m2_2_rigorous ;
            #endif     

            // ofs<<mu_m2_rigorous<<" "<<mu_m2_1_rigorous<<" "<<mu_m2_2_rigorous<<"\n" ;
=======
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
>>>>>>> 2b7bfc688d139f55b0773e421cf6a0444ffe9f84

            double mu_t2_est = mu_hbar+beta*mu_m2_rigorous;
            double var_t_est = mu_t2_est-pow(mu_t_est,2);
            if(var_t_est<0)var_t_est = 0.;

<<<<<<< HEAD
            #if dimensionless
                mu_t = kt_bar*(1+mu_t) ;
                mu_t_est = kt_bar*(1+mu_t_est) ;
                var_t = kt_bar*kt_bar*(var_t) ;
                var_t_est = kt_bar*kt_bar*(var_t_est) ;
            #endif

            mu_t_cas[jj] = mu_t ;
            mu_t_est_cas[jj] = mu_t_est ;
            std_t_cas[jj] = var_t ;
            std_t_est_cas[jj] = var_t_est ;

            // ofs<<var_t_est<<" "<<mu_t2_est<<"\n" ;

        }
        
=======
            mu_t_cas[ii][i] = mu_t ;
            mu_t_est_cas[ii][i] = mu_t_est ;
            std_t_cas[ii][i] = sqrt(var_t) ;
            std_t_est_cas[ii][i] = sqrt(var_t_est) ;
        }
    
    }

    for (int ii = 0; ii < targetNum; ii++)
    {
>>>>>>> 2b7bfc688d139f55b0773e421cf6a0444ffe9f84
        double sum1_mut = 0.;
        double sum2_mut = 0.;
        double sum1_mutest = 0.;
        double sum2_mutest = 0.;
        double sum1_stdt = 0.;
        double sum2_stdt = 0.;
        double sum1_stdtest = 0.;
        double sum2_stdtest = 0.;
<<<<<<< HEAD
        for(int j = 0; j < set; j++){
            sum1_mut += mu_t_cas[j];
            sum2_mut += pow(mu_t_cas[j],2);
            sum1_mutest += mu_t_est_cas[j];
            sum2_mutest += pow(mu_t_est_cas[j],2);
            sum1_stdt += std_t_cas[j];
            sum2_stdt += pow(std_t_cas[j],2);
            sum1_stdtest += std_t_est_cas[j];
            sum2_stdtest += pow(std_t_est_cas[j],2);
        }

        sum1_mut /= set;
        sum1_mutest /= set;
        sum1_stdt /= set;
        sum1_stdtest /= set;

        double var_mut = (sum2_mut/set-pow(sum1_mut,2))*set/(set-1);
        double var_mutest = (sum2_mutest/set-pow(sum1_mutest,2))*set/(set-1);
        double var_stdt = (sum2_stdt/set-pow(sum1_stdt,2))*set/(set-1);
        double var_stdtest = (sum2_stdtest/set-pow(sum1_stdtest,2))*set/(set-1);
=======
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
>>>>>>> 2b7bfc688d139f55b0773e421cf6a0444ffe9f84
        if(var_mutest < 0.)var_mutest = 0.;
        if(var_stdtest < 0.)var_stdtest = 0.;    

        double ur_mu = sqrt(var_mutest)/sqrt(var_mut);    
        double ur_std = sqrt(var_stdtest)/sqrt(var_stdt); 

        #if mean
<<<<<<< HEAD
        ofs<<target<<" "<<ur_mu<<"\n";
        #endif

        #if deviation
        ofs<<target<<" "<<ur_std<<"\n";
        #endif

        #if test
        ofs<<target<<" "<<sum1_mut<<" "<<sum1_mutest<<" "<<sum1_stdt<<" "<<sum1_stdtest<<"\n" ;
        #endif

        // double a0_est = sum_a0/set ;
        // double a1_est = sum_a1/set ;

        // ofs<<target<<" "<<a0_est<<"\n";


    }
    
    
    
    
=======
        ofs<<target[ii]<<" "<<sum1_mut<<"\n";
        #endif

        #if deviation
        ofs<<target[ii]<<" "<<ur_std<<"\n";
        #endif
    }
    

>>>>>>> 2b7bfc688d139f55b0773e421cf6a0444ffe9f84
    return 0 ;
};