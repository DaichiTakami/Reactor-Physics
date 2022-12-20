#define _USE_MATH_DEFINES
#include<cmath>
#include<iostream>
#include<vector>
#include<stdlib.h>
#include<stdio.h>
#include<fstream>
#include<math.h>

#include<Windows.h>
#include<string.h>



using namespace std ;



int main(){

    WIN32_FIND_DATA ffd ;
    HANDLE hNextFile = FindFirstFileW("C:\\Users\\rokoroko10\\OneDrive\\デスクトップ\\プログラミング\\Reactor-Physics\\SAMPLE\\*.*" , &ffd) ;
    if (hNextFile != INVALID_HANDLE_VALUE) 
    {
        do
        {
            puts(ffd.cFileName) ;
        } while (FindNextFile(hNextFile, &ffd));
        FindClose(hNextFile) ;
        
    }
    

    return 0 ;
};