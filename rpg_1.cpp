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
using System.Collections.Generic


vector<string> split(string& input, char delimiter)
{
    istringstream stream(input) ;
    string field ;
    vector<string> result ;
    while (getline(stream,field,delimiter))
    {
        result.push_back(field) ;
    }
    return result ;
}


int main(){

    ifstream ifs("sigmaData.csv") ;

    string line ;

    int energyGroupNumber ;

     while (getline(ifs, line)) {
        
        vector<string> strvec = split(line, ',');

        energyGroupNumber = getline(ifs, line) ;
        
        for (int i=0; i<strvec.size();i++){
            printf("%5d\n", stoi(strvec.at(i)));
        }
        
    }

    double sigma[energyGroupNumber][energyGroupNumber] ;


    


    return 0 ;
};