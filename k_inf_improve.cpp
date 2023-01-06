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

int main(){

    srand(10) ;

    ifstream ifs ;

    char openFile[256] ;
    cout<<"サンプルのあるフォルダ名を入力してください。"<<"\n" ;
    cin>>openFile ;

    const char *openFolder = openFile ;

    char output[256] ;
    cout<<"出力するファイルの名前を入力してください。"<<"\n" ;
    cout<<"例:a.dat,b"<<"\n" ;
    cin>>output;

    const char *outputfile = output ;
    ofstream ofs(outputfile) ;

    int fileNum = 0 ;

    int cal  ;

    char filename[256] ;
    char buf[256] ;

    int targetNum = 0 ;

    sprintf(filename, "%s/cal_point",openFolder) ;

    ifs.open(filename);

    if (ifs.fail())
    {
        cout<<"計算目標を指定するファイルが見つかりませんでした。"<<"\n";
        cout<<"目標パラメータの数を入力してください"<<"\n" ;
        cin>>targetNum ;

        double *cal_point ;
        cal_point = new double[targetNum] ;

        cout<<"出力データのxの値として出力する値を指定してください。"<<"\n" ;
        double cal_buf ;

        for (int jj = 0; jj < targetNum; jj++)
        {
            cout<<jj+1<<"番目 : ";
            cin>>cal_buf   ;
            cal_point[jj] = cal_buf ;
        }
        
    }

    // for (int i = 0; i < mockupNum; i++)
    // {
    //     cout<<cal_point[i]<<"\n" ;
    // }
    

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

    sprintf(filename,"%s/sample_0",openFolder) ;

    ifs.open(filename);

    int linenum = 0 ;

    while (ifs.getline(buf,sizeof(buf)))
    {
        linenum ++ ;
    }

    ifs.close() ;


    int mockupNum = linenum - targetNum ;

    int mockupfrag ;
    cout<<"類似パラメータはサンプル中の下から"<<mockupNum<<"個です。"<<"\n" ;
    cout<<"これらのの類似パラメータを使用しますか？"<<"\n" ;
    cout<<"1:はい　2:いいえ"<<"\n" ;
    cin>>mockupfrag ;

    if (mockupfrag == 2)
    {
        cout<<"いくつ使用しますか？"<<"\n" ;
        cin>>mockupNum ;
    }
    
    double ktSum ;
    vector<double> km_iSum(mockupNum) ;

    double ktSum2 ;
    double km_iSum2[mockupNum][mockupNum] ;

    double ktSum4 ;
    vector<double> km_iSum4(mockupNum) ;

    double ktkm_iSum[mockupNum] ;
    double ktkm_iSum2[mockupNum] ;

    int mockupdata[mockupNum] ;
    int mockupbuf ;

    if (mockupfrag == 2)
    {
        cout<<"どのデータを類似パラメータとして利用しますか？"<<"\n" ;
        cout<<"サンプルデータの上から何番目を利用するか指定してください。"<<"\n" ;
        for (int i = 0; i < mockupNum; i++)
        {
            cout<<i+1<<"つ目 : " ;
            cin>>mockupbuf ;
            mockupdata[i] = mockupbuf-1 ;
        }
        
    }

    // for (int i = 0; i < mockupNum; i++)
    // {
    //     cout<<mockupdata[i]<<"\n" ;
    // }
    
    

    for (int ii = 0; ; ii++)
    {
        
        sprintf(filename, "%s/sample_%d",openFolder,ii);

        ifs.open(filename) ;
        if (ifs.fail())
        {
            break ;
        }

        ifs.close() ;

        fileNum ++ ;
    }

    cout<<"サンプルファイルの数は"<<fileNum<<"です。"<<"\n" ;

    int recNum ;
    if (fileNum < 100)
    {
        recNum = 10*(fileNum/10) ;
    }else if (fileNum < 1000)
    {
        recNum = 100*(fileNum/100) ; ;
    }else{
        recNum = 1000*(fileNum/1000) ;
    }
    

    int sampleNum ;
    cout<<"計算で使用するサンプル数を入力してください"<<"\n" ;
    cout<<"推奨値は"<<recNum<<"です。"<<"\n" ;
    cin>>sampleNum ;

    int num[sampleNum] ;

    cout<<"計算回数を入力してください"<<"\n" ;
    cout<<"推奨値は100です。"<<"\n" ;
    cin>>cal ;

    int flag ;
    cout<<"出力形式を選択してください。"<<"\n";
    cout<<"1.UR 2.期待値と標準誤差"<<"\n" ;
    cin>>flag ;

    int outputFlag ;

    if (flag == 1)
    {
        cout<<"出力するパラメータを選択してください。"<<"\n" ;
        cout<<"1.平均 2.標準偏差 3.相関"<<"\n" ;
        cin>>outputFlag ;
        if (outputFlag != 1 && outputFlag != 2 && outputFlag != 3)
        {
            cout<<"無効な値が入力されました。"<<"\n" ;
            cout<<"プログラムを終了します。"<<"\n" ;
            return 0 ;
        }
        
    }else if (flag == 2)
    {
        cout<<"出力するパラメータを選択してください。"<<"\n" ;
        cout<<"1.平均 2.標準偏差"<<"\n" ;
        cin>>outputFlag ;
        if (outputFlag != 1 && outputFlag != 2)
        {
            cout<<"無効な値が入力されました。"<<"\n" ;
            cout<<"プログラムを終了します。"<<"\n" ;
            return 0 ;
        }
    } else {
        cout<<"無効な値が入力されました。"<<"\n" ;
        cout<<"プログラムを終了します。"<<"\n" ;
        return 0 ;
    }

    double mu_t_cas[cal] ;
    double mu_t_est_cas[cal] ;
    double std_t_cas[cal] ;
    double std_t_est_cas[cal] ;
    double corr_tm_cas[cal] ;

    double kt[targetNum][fileNum] ;
    double km_i[mockupNum][fileNum] ;

    double ref[linenum] ;
    for (int ii = 0; ii < linenum; ii++)
    {
        ref[ii] = 0. ;
    }
    

    for (int ii = 0; ii < fileNum; ii++)
    {
        
        sprintf(filename, "%s/sample_%d",openFolder, ii);

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
            ref[jj] += arr[jj] ;
        }
    }
    
    for (int ii = 0; ii < linenum; ii++)
    {
        ref[ii] /=fileNum ;
    }
    

    for (int ii = 0; ii < fileNum; ii++)
    {
        
        sprintf(filename, "%s/sample_%d",openFolder, ii);

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

        if (mockupfrag == 2)
        {
            for (int jj = 0; jj < mockupNum; jj++)
            {
                km_i[jj][ii] = arr[mockupdata[jj]] ;
            }
            
        }else{
            for (int jj = 0; jj < mockupNum; jj++)
            {
                km_i[jj][ii] = arr[jj+targetNum] ;
                
            }
        }
           
    }

    double cov_mm_rigorous[mockupNum][mockupNum] ;

    
    for (int jj = 0; jj < mockupNum; jj++)
    {
        km_iSum[jj] = 0. ;
        km_iSum4[jj] = 0. ;
        for (int kk = 0; kk < mockupNum; kk++)
        {
            km_iSum2[jj][kk] = 0. ;
        }
        
    }
    
    for (int ii = 0; ii < fileNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            km_iSum[jj] += km_i[jj][ii] ;
            km_iSum4[jj] += pow(km_i[jj][ii],4) ;
            for (int kk = 0; kk < mockupNum; kk++)
            {
                km_iSum2[jj][kk] += km_i[jj][ii]*km_i[kk][ii] ;
            }
            
        }
        
    }

    for (int ii = 0; ii < mockupNum; ii++)
    {
        for (int jj = 0; jj < mockupNum; jj++)
        {
            cov_mm_rigorous[ii][jj] = (km_iSum2[ii][jj]/fileNum - km_iSum[ii]/fileNum * km_iSum[jj]/fileNum)*fileNum/(fileNum-1) ;
        }
        
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
                    ktkm_iSum[kk] += kt[i][num[ii]]*km_i[kk][num[ii]] ;
                    ktkm_iSum2[kk] += pow(kt[i][num[ii]]*km_i[kk][num[ii]],2) ;
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

            // cout<<cov_mm[0][0]<<" "<<cov_mm[1][1]<<" "<<cov_mm[1][0]<<"\n" ;
            
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
            // cout<<a[0]<<" "<<a[1]<<"\n" ;

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
                km[jj] = 0. ;

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

            // ofs<<corr_tm_vi<<"\n" ;

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
            double mu_m2_rigorous  = 0. ;
            for (int ii = 0; ii < mockupNum; ii++)
            {
                for (int jj = 0; jj < mockupNum; jj++)
                {
                    mu_m2_rigorous += a[ii]*a[jj]*cov_mm_rigorous[ii][jj] ;
                }
                
            }
            


            double mu_t2_est = mu_hbar+beta*mu_m2_rigorous;
            double var_t_est = mu_t2_est-pow(mu_t_est,2);
            if(var_t_est<0)var_t_est = 0.;

            mu_t_cas[j] = mu_t ;
            mu_t_est_cas[j] = mu_t_est ;
            std_t_cas[j] = sqrt(var_t) ;
            std_t_est_cas[j] = sqrt(var_t_est) ;
            corr_tm_cas[j] = corr_tm_vi ;
            
        
        }

        
        double sum1_mut = 0.;
        double sum2_mut = 0.;
        double sum1_mutest = 0.;
        double sum2_mutest = 0.;
        double sum1_stdt = 0.;
        double sum2_stdt = 0.;
        double sum1_stdtest = 0.;
        double sum2_stdtest = 0.;
        double sum_corr = 0. ;
        for(int j = 0; j < cal; j++){
            sum1_mut += mu_t_cas[j];
            sum2_mut += pow(mu_t_cas[j],2);
            sum1_mutest += mu_t_est_cas[j];
            sum2_mutest += pow(mu_t_est_cas[j],2);
            sum1_stdt += std_t_cas[j];
            sum2_stdt += pow(std_t_cas[j],2);
            sum1_stdtest += std_t_est_cas[j];
            sum2_stdtest += pow(std_t_est_cas[j],2);
            sum_corr += corr_tm_cas[j] ;
        }

        sum1_mut /= cal;
        sum1_mutest /= cal;
        sum1_stdt /= cal;
        sum1_stdtest /= cal;
        sum_corr /= cal ;

        double var_mut = (sum2_mut/cal-pow(sum1_mut,2))*cal/(cal-1);
        double var_mutest = (sum2_mutest/cal-pow(sum1_mutest,2))*cal/(cal-1);
        double var_stdt = (sum2_stdt/cal-pow(sum1_stdt,2))*cal/(cal-1);
        double var_stdtest = (sum2_stdtest/cal-pow(sum1_stdtest,2))*cal/(cal-1);
        if(var_mutest < 0.)var_mutest = 0.;
        if(var_stdtest < 0.)var_stdtest = 0.;    

        double ur_mu = sqrt(var_mutest)/sqrt(var_mut);    
        double ur_std = sqrt(var_stdtest)/sqrt(var_stdt); 

        if (flag == 1)
        {
            if (outputFlag == 1)
            {
                ofs<<cal_point[i]<<" "<<ur_mu<<"\n" ;
            }

            if (outputFlag == 2)
            {
                ofs<<cal_point[i]<<" "<<ur_std<<"\n" ;
            }
            
            if (outputFlag == 3)
            {
                ofs<<cal_point[i]<<" "<<sum_corr<<"\n" ;
            }
        }
        
        if (flag == 2)
        {
            if (outputFlag == 1)
            {
                ofs<<cal_point[i]<<" "<<sum1_mutest<<" "<<var_mutest<<"\n" ;
            }

            if (outputFlag == 2)
            {
                ofs<<cal_point[i]<<" "<<sum1_stdtest<<" "<<var_stdtest<<"\n" ;
            }
            
        }
        
        
        
        
    }
    
    // for (int i = 0; i < mockupNum; i++)
    // {
    //     for (int j = 0; j < mockupNum; j++)
    //     {
    //         cout<<cov_mm_rigorous[i][j]<<" ";
    //     }
    //     cout<<"\n" ;
    // }
    
    

        
    

    return 0 ;
};