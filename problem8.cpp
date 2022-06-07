#include<iostream>
#include<cmath>
#include<vector>
#include<stdlib.h>
#include<string>
#include<cstdlib>

using namespace std;

class student{
    private: 
    int intelligence ;
    int health ;
    bool stdnt ;

    public:
    student(bool stdntimp = true){intelligence = 0;health = 10; stdnt = stdntimp;};
    ~student(){};
    int GetInt(){return intelligence;};
    int GetH(){return health;};
    
    void PutInt(int inp){intelligence = inp;};
    void Seminar(){intelligence++;};
    void DrinkingParty(){health--;};

    void Makestudent(){stdnt=true;};
    void MakeTeacher(){stdnt=false;};

    bool Student(){return stdnt;};
    
    bool WritingADissertation(student &other){
        if (intelligence<10||(stdnt&&other.Student()))
        {
            return false ;
        }
      return true ;  
    };

};

int main(){

    student boku(true);
    boku.PutInt(8);
    boku.DrinkingParty();
    boku.Seminar();
    boku.Seminar();

    student watashi(false);
    

    cout<<"#Intelligence of boku   : "<<boku.GetInt()<<"\n";
    cout<<"#Health of boku   : "<<boku.GetH()<<"\n";

    bool judge = boku.WritingADissertation(watashi);
    

    if (judge)
    {
        cout<<"Writing a dissertation."<<"\n";
    }else{
            cout<<"# No dissertation.\n";
    }
    

    
        
    

    return 0 ;
};