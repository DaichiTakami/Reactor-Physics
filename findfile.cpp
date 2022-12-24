#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>

using namespace std ;

int main(){

    ifstream ifs ;

    const char *inputfile = "a.dat" ;
    ofstream ofs(inputfile) ;

    int fileNum = 10000 ;

    char filename[256] ;
    char buf[256] ;

    for (int ii = 0; ii < fileNum; ii++)
    {
        sprintf(filename, "/Users/takamidaichi/Documents/GitHub/Reactor-Physics/SAMPLE/sample_%d", ii);

        // ofs<<filename<<"\n" ;

        ifs.open(filename) ;
        if (ifs.fail())
        {
            ofs<<"can not open file."<<"\n" ;
            exit(1);
        }
        
        int linenum = 0 ;

        while (ifs.getline(buf,sizeof(buf)))
        {
            linenum ++ ;
        }
        
        ofs<<"line number = "<<linenum<<"\n" ;

        ifs.clear() ;
        ifs.seekg(0, std::ios::beg) ;

        double *arr ;
        arr = new double[linenum] ;

        for (int jj = 0; jj < linenum; jj++)
        {
            ifs.getline(buf,sizeof(buf)) ;
            arr[jj] = atof(buf) ;
        }

        ifs.close() ;
        
    }
    
    

    return 0 ;
};