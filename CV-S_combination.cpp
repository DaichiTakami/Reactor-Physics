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

    return mu + sigma*z
}

double variance (double E_x2, double E2_x)
{
    return E2_x - E2_x ;
}

int main(){
    
    const char *filename = "a.dat" ;
    ofstream ofs(filename) ;
    
    srand(10) ;

    int samples = 100 ;
    int sets = 100 ;
    int loops = 101 ;

    int inps = 2 ;
    int oups = 1 ;
    int mocks = 2 ;

    for (int ii = 0; ii < loops; ii++)
    {
        #if case == 1
            // case 1
            double k_t_bar = ii * 0.01 ;
            double k_m_bar[mocks] = {0.5,0.5} ;
            double s_t1[inps] = {0.1,0.2} ;
            double s_t2[inps] = {0.,0.} ;
            double s_m1[mocks][inps] = {{0.1,0.1},
                                        {0.2,0.2}} ;
            double target = k_t_bar ;
        #endif
        
        double mean_inp[inps] = {0.,0.} ;
        double std_inp[inps] = {0.01,0.02} ;

        vector<double>  mu_t_cas(samples);
        vector<double>  mu_t_est_cas(samples);
        vector<double>  std_t_cas(samples);
        vector<double>  std_t_est_cas(samples);

        for (int jj = 0; jj < sets; jj++)
        {
            double sample[inps][samples] ;
            for (int i = 0; i < inps; i++)
            {
                for (int j = 0; j < samples; j++)
                {
                    sample[i][j] = gauss(mean_inp[i],std_inp[i]);
                }
                
            }

            
            
        }
        

    }
    
    
    
    
    return 0 ;
};