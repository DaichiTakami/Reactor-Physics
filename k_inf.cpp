#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>
#include<string>

using namespace std ;


int main(){

    const char *filename = "a.dat" ;
    ofstream ofs(filename) ;

    FILE *fp ;

    const int sampleNum = 10000 ;

    char openfile[100] ;

    for (int ii = 0; ii < sampleNum; ii++)
    {
        
        sprintf(openfile,"sample_%d",ii);
        fp = fopen(openfile,"r") ;

        



        fclose(fp);
    }
    

    

    
    

    return 0 ;
};