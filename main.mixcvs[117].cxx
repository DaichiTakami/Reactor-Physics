#include<vector>
#include<cmath>
#include<cstdlib>
#include<iostream>

using namespace std;

#define conventional 0
#define dimensionless 1



// To generate random number following the normal distribution

double ransu_gauss(double mu, double sigma)
{
  double fact=1./RAND_MAX;
  //double fact=1./2147483647;    
  double PI =3.1415926535897932384626433;
  double PI2=PI*2.;  
 
  double alpha=double(rand())*fact;
  double beta=double(rand())*fact*PI2;
  
  double BoxMuller1=sqrt(-2*log(alpha));
  double BoxMuller2=sin(beta);

  return sigma*(BoxMuller1*BoxMuller2)+mu;
};


  int main()
{
  double k_t_bar=1.0;
  double k_m1_bar=1.0;
  double k_m2_bar=1.0;
  double s_t_1=0.1;
  double s_t_2=0.1;
  double ss_t_1=0.;
  double ss_t_2=0.;
  double s_m1_1=0.1;
  double s_m1_2=0.2;  
  double s_m2_1=0.2;
  double s_m2_2=0.1;

  // weight for combined parameter
  double a1=10.;
  double a2=10.;  

  double std_inp1=0.01; // sigma for Z1
  double std_inp2=0.01; // sigma for Z2
 
  // +++++++++++++++++++++++

  int num=10000; // The number of samples

  int pnt_num=1;
  int pnt_array[]={num};

  int cas=1;
  //int cas=100;

  /*
  vector< vector<double> > mu_t_cas(cas);
  vector< vector<double> > mu_t_est_cas(cas);
  vector< vector<double> > std_t_cas(cas);
  vector< vector<double> > std_t_est_cas(cas);
  */
  
  for(int jj=0;jj<cas;jj++){

  vector<double> z1_list(num);
  vector<double> z2_list(num);  
  for(int i=0;i<num;i++){
    z1_list[i]=ransu_gauss(0, std_inp1);
    z2_list[i]=ransu_gauss(0, std_inp2);    
  };

    /*  
  vector<double> mu_t_store(num);
  vector<double> mu_t_est_store(num);
  vector<double> std_t_store(num);
  vector<double> std_t_est_store(num);
    */
    
  
  double sum1_t=0.;
  double sum2_t=0.;
  double sum4_t=0.;
  double sum1_m=0.;
  double sum2_m=0.;
  double sum4_m=0.;
  double sum1_tm=0.;
  double sum2_tm=0.;
  
  for(int i=0;i<num;i++){

    double z1=z1_list[i];
    double z2=z2_list[i];

    
    double k_m_bar=k_m1_bar; 
    double s_m_1=s_m1_1;
    double s_m_2=s_m1_2;    
    
#if conventional
    double k_t=k_t_bar*(1.+s_t_1*z1+s_t_2*z2+ss_t_1*z1*z1+ss_t_2*z2*z2);
    double k_m=k_m_bar*(1.+s_m_1*z1+s_m_2*z2);
#endif    

#if dimensionless    
    double k_t=s_t_1*z1+s_t_2*z2+ss_t_1*z1*z1+ss_t_2*z2*z2;
    //double k_m=s_m_1*z1+s_m_2*z2; // single mock-up
    double k_m=(s_m1_1*a1+s_m2_1*a2)*z1+(s_m1_2*a1+s_m2_2*a2)*z2; // combined mock-up
#endif

    //cout<<k_t<<" "<<k_m1<<" "<<k_m2<<"\n";

    sum1_t+=k_t;
    sum2_t+=pow(k_t,2);
    sum4_t+=pow(k_t,4);

    sum1_m+=k_m;
    sum2_m+=pow(k_m,2);
    sum4_m+=pow(k_m,4);

    sum1_tm+=k_t*k_m;
    sum2_tm+=pow(k_t*k_m,2);
    
    double mu_t=sum1_t/(i+1);
    double mu_m=sum1_m/(i+1);
    double mu2_t=sum2_t/(i+1);
    double mu2_m=sum2_m/(i+1);

    double var_t=(mu2_t-pow(mu_t,2))*(i+1)/i;
    double var_m=(mu2_m-pow(mu_m,2))*(i+1)/i;

    double var2_t=(sum4_t/(i+1)-pow(mu2_t,2))*(i+1)/i;
    double var2_m=(sum4_m/(i+1)-pow(mu2_m,2))*(i+1)/i;

    double cov_tm=(sum1_tm/(i+1)-mu_t*mu_m)*(i+1)/i;
    double cor_tm=cov_tm/(sqrt(var_t)*sqrt(var_m));
    double cov2_tm=(sum2_tm/(i+1)-mu2_t*mu2_m)*(i+1)/i;
    double cor2_tm=cov2_tm/(sqrt(var2_t)*sqrt(var2_m));

    // ... CV
    double alpha=cov_tm/var_m; 
    double beta=cov2_tm/var2_m;
      
    vector<double> sample_h(i+1);
    vector<double> sample_hbar(i+1);

    for(int ii=0;ii<=i;ii++){

      double zz1=z1_list[ii];
      double zz2=z2_list[ii];

#if conventional    
      double kk_t=k_t_bar*(1.+s_t_1*zz1+s_t_2*zz2+ss_t_1*zz1*zz1+ss_t_2*zz2*zz2);
      double kk_m=k_m_bar*(1.+s_m_1*zz1+s_m_2*zz2);
#endif

#if dimensionless    
      double kk_t=s_t_1*zz1+s_t_2*zz2+ss_t_1*zz1*zz1+ss_t_2*zz2*zz2;
      //double kk_m=s_m_1*zz1+s_m_2*zz2; // single mock-up
      double kk_m=(s_m1_1*a1+s_m2_1*a2)*zz1+(s_m1_2*a1+s_m2_2*a2)*zz2; // combined mock-up
#endif
   
      sample_h[ii]=kk_t-alpha*kk_m;
      sample_hbar[ii]=pow(kk_t,2)-beta*pow(kk_m,2);
    };

    double sum1_h=0.;
    double sum2_h=0.;
    double sum1_hbar=0.;
    double sum2_hbar=0.;
    for(int ii=0;ii<=i;ii++){
      sum1_h+=sample_h[ii];
      sum2_h+=pow(sample_h[ii],2);
      sum1_hbar+=sample_hbar[ii];
      sum2_hbar+=pow(sample_hbar[ii],2);
    };
    double mu_h=sum1_h/(i+1);
    double mu_hbar=sum1_hbar/(i+1);

#if conventional    
    double mu_t_est=mu_h+alpha*k_m_bar; // [MU] estimated by CV
    double mu_m2_rigorous=k_m_bar*k_m_bar;
    mu_m2_rigorous+=k_m_bar*k_m_bar*(s_m_1*s_m_1*std_inp1*std_inp1);
    mu_m2_rigorous+=k_m_bar*k_m_bar*(s_m_2*s_m_2*std_inp2*std_inp2);
#endif
    
#if dimensionless    
    double mu_t_est=mu_h; // [MU] estimated by CV
    /* 
    // ... single mock-up
    double mu_m2_rigorous=s_m_1*s_m_1*std_inp1*std_inp1;
    mu_m2_rigorous+=s_m_2*s_m_2*std_inp2*std_inp2;
    */
    // ... combined mock-up
    double mu_m2_rigorous=pow((s_m1_1*a1+s_m2_1*a2)*std_inp1,2);
    mu_m2_rigorous+=pow((s_m1_2*a1+s_m2_2*a2)*std_inp2,2);
#endif

    double mu_t2_est=mu_hbar+beta*mu_m2_rigorous;    
    double var_t_est=mu_t2_est-pow(mu_t_est,2); // [VAR] estimated by CV
    if(var_t_est<0)var_t_est=0.;
    
#if dimensionless    
    mu_t=k_t_bar*(1+mu_t);
    mu_t_est=k_t_bar*(1+mu_t_est);
    var_t=k_t_bar*k_t_bar*(var_t);
    var_t_est=k_t_bar*k_t_bar*(var_t_est);      
#endif    

    /*
      mu_t_store[i]=mu_t;
      mu_t_est_store[i]=mu_t_est;
      std_t_store[i]=sqrt(var_t);
      std_t_est_store[i]=sqrt(var_t_est);
      */
      /*
      cout<<i<<" "<<mu_x<<" "<<mu_y<<" ";
      cout<<sqrt(var_x)<<" "<<sqrt(var_y)<<" ";
      cout<<cov_xy/(sqrt(var_x)*sqrt(var_y))<<" ";
      cout<<alpha<<" "<<beta<<" ";
      cout<<mu_x_est<<" "<<var_x_est<<" ";
      cout<<"\n";
      */

    cout<<i+1<<" ";
    cout<<mu_t<<" "<<var_t<<" ";
    cout<<mu_t_est<<" "<<var_t_est<<" ";
    //cout<<mu_m2_rigorous-k_m_bar*k_m_bar<<" ";
    //cout<<mu_m1<<" "<<var_m1<<" ";    
    //cout<<mu_m2<<" "<<var_m2<<" ";
    //cout<<alp_m1<<" "<<bet_m1<<" ";
    //cout<<alp_m2<<" "<<bet_m2<<" ";    
      //cout<<sqrt(var_t)<<" "<<sqrt(var_t_est)<<" ";
    cout<<"\n";
    
  };


#if 0  
      double mu_t=sum1_t/num;
      double mu_m=sum1_m/num;
      double mu_t2=sum2_t/num;
      double mu_m2=sum2_m/num;

      double var_t=(sum2_t/num-pow(mu_t,2))*num/(num-1);
      double var_m=(sum2_m/num-pow(mu_m,2))*num/(num-1);

      double var_t2=(sum4_t/num-pow(mu_t2,2))*num/(num-1);
      double var_m2=(sum4_m/num-pow(mu_m2,2))*num/(num-1);

      double cov_tm=(sum1_tm/num-mu_t*mu_m)*num/(num-1);
      double corr_tm=cov_tm/(sqrt(var_t)*sqrt(var_m));

      double cov_t2m2=(sum2_tm/num-mu_t2*mu_m2)*num/(num-1);
      double corr_t2m2=cov_t2m2/(sqrt(var_t2)*sqrt(var_m2));

      double alpha=cov_tm/var_m; 
      double beta=cov_t2m2/var_m2;

      //cout<<i<<" "<<mu_t<<" "<<mu_m<<" "<<sqrt(var_t)<<" "<<sqrt(var_m)<<" "<<corr_tm<<"\n";
      //cout<<alpha<<" "<<beta<<"\n";

      // ... CV
      vector<double> sample_h(num);
      vector<double> sample_hbar(num);
      for(int ii=0;ii<num;ii++){
        double ds_s=sample[ii];
	double k_t_ii=k_t_bar*(1.+s_t1*ds_s+s_t2*ds_s*ds_s);
	double k_m_ii=k_m_bar*(1.+s_m1*ds_s);
        sample_h[ii]=k_t_ii-alpha*k_m_ii;
        sample_hbar[ii]=pow(k_t_ii,2)-beta*pow(k_m_ii,2);
      };
      double sum1_h=0.;
      double sum2_h=0.;
      double sum1_hbar=0.;
      double sum2_hbar=0.;
      for(int ii=0;ii<num;ii++){
        sum1_h+=sample_h[ii];
        sum2_h+=pow(sample_h[ii],2);
        sum1_hbar+=sample_hbar[ii];
        sum2_hbar+=pow(sample_hbar[ii],2);
      };
      double mu_h=sum1_h/(num);
      double mu_hbar=sum1_hbar/(num);

      double mu_t_est=mu_h+alpha*k_m_bar;         // estimated by CV
      double mu_m2_rigorous=k_m_bar*k_m_bar*s_m1*s_m1*std_inp*std_inp+k_m_bar*k_m_bar;
      double mu_t2_est=mu_hbar+beta*mu_m2_rigorous;
      double var_t_est=mu_t2_est-pow(mu_t_est,2); // estimated by CV
      if(var_t_est<0)var_t_est=0.;

      mu_t_store[i]=mu_t;
      mu_t_est_store[i]=mu_t_est;
      std_t_store[i]=sqrt(var_t);
      std_t_est_store[i]=sqrt(var_t_est);
      /*
      cout<<i<<" "<<mu_x<<" "<<mu_y<<" ";
      cout<<sqrt(var_x)<<" "<<sqrt(var_y)<<" ";
      cout<<cov_xy/(sqrt(var_x)*sqrt(var_y))<<" ";
      cout<<alpha<<" "<<beta<<" ";
      cout<<mu_x_est<<" "<<var_x_est<<" ";
      cout<<"\n";
      */
      /*
      cout<<num<<" ";
      cout<<mu_t<<" "<<mu_t_est<<" ";
      cout<<sqrt(var_t)<<" "<<sqrt(var_t_est)<<" ";
      cout<<"\n";
      */
#endif

#if 0      
  for(int i=0;i<pnt_num;i++){
    int tmp=pnt_array[i]-1;
    mu_t_cas[jj].push_back(mu_t_store[tmp]);
    mu_t_est_cas[jj].push_back(mu_t_est_store[tmp]);
    std_t_cas[jj].push_back(std_t_store[tmp]);
    std_t_est_cas[jj].push_back(std_t_est_store[tmp]);
  };
#endif

  }; // The end of the loop [for(int jj=0;jj<cas;jj++)]

#if 0

  for(int i=0;i<pnt_num;i++){

    double sum1_mut=0.;
    double sum2_mut=0.;
    double sum1_mutest=0.;
    double sum2_mutest=0.;
    double sum1_stdt=0.;
    double sum2_stdt=0.;
    double sum1_stdtest=0.;
    double sum2_stdtest=0.;
    for(int j=0;j<cas;j++){
      sum1_mut+=mu_t_cas[j][i];
      sum2_mut+=pow(mu_t_cas[j][i],2);
      sum1_mutest+=mu_t_est_cas[j][i];
      sum2_mutest+=pow(mu_t_est_cas[j][i],2);
      sum1_stdt+=std_t_cas[j][i];
      sum2_stdt+=pow(std_t_cas[j][i],2);
      sum1_stdtest+=std_t_est_cas[j][i];
      sum2_stdtest+=pow(std_t_est_cas[j][i],2);
    };

    sum1_mut/=cas;
    sum1_mutest/=cas;
    sum1_stdt/=cas;
    sum1_stdtest/=cas;

    double var_mut=(sum2_mut/cas-pow(sum1_mut,2))*cas/(cas-1);
    double var_mutest=(sum2_mutest/cas-pow(sum1_mutest,2))*cas/(cas-1);
    double var_stdt=(sum2_stdt/cas-pow(sum1_stdt,2))*cas/(cas-1);
    double var_stdtest=(sum2_stdtest/cas-pow(sum1_stdtest,2))*cas/(cas-1);
    if(var_mutest<0.)var_mutest=0.;
    if(var_stdtest<0.)var_stdtest=0.;    
    /*
    cout<<pnt_array[i]<<" ";
    cout<<sum1_mut<<" "<<sqrt(var_mut)<<" ";
    cout<<sum1_mutest<<" "<<sqrt(var_mutest)<<" ";
    cout<<sum1_stdt<<" "<<sqrt(var_stdt)<<" ";
    cout<<sum1_stdtest<<" "<<sqrt(var_stdtest)<<" ";
    cout<<sqrt(var_mutest)/sqrt(var_mut)<<" ";
    cout<<sqrt(var_stdtest)/sqrt(var_stdt)<<" ";
    cout<<"\n";
    */
    double ur_mu=sqrt(var_mutest)/sqrt(var_mut);
    double ur_std=sqrt(var_stdtest)/sqrt(var_stdt);
    //if(i==pnt_num-1)cout<<pnt_array[i]<<" "<<ur_mu<<" "<<ur_std<<"\n";
    if(i==pnt_num-1)cout<<target<<" "<<ur_mu<<" "<<ur_std<<"\n";
    //if(i==pnt_num-1)cout<<" "<<var_mutest<<" "<<var_mut<<" "<<ur_std<<"\n";        
  };

  };
#endif

  return 0;
}

