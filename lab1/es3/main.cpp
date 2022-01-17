#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include "lib.h"
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/block_stat.h"
#include <string>
#include <cmath>


using namespace std;

int main (int argc, char *argv[]){

        //INITIALISING SOME VARIABLES                                                   

        Random rnd;          //random generator     
        rnd.Initialize();    //initializing random generator
        const string file_out="stima_pi.csv";   //output file

        int N_blocks=100;         //# of blocks
        int N_throws=10000;       //# of throws
        int n=N_throws/N_blocks;  //# of throws in each block
        BlockStat block_stat=BlockStat(N_blocks, N_throws);     //initializing class for block statistics

        int M=13000;        //number of needles 
        double d=1.;        //distance between lines
        double L=0.83;      //lenght of the needle
        
        double ave[N_blocks]{};   //array of mean in each block
        double av2[N_blocks]{};   //array of mean squared in each block


        //EVALUATE AVERAGE IN EACH BLOCK
        cout<<"Evaluating pi with the buffon experiment"<<endl;

        for(int i=0; i<N_blocks; i++){
                double sum=0;
                for(int j=0;j<n;j++){ 
                        sum+=pi(M,L,d,rnd);         //cumulative sum in each block
                }
                ave[i]=sum/n;                       //estimating mean for each block
                av2[i]=pow(ave[i],2);               //estimating mean squared for each block
        }


        //BLOCK STATISTICS
        cout<<"Doing block statistics"<<endl;

        block_stat.progressive_mean(ave);           //progressive mean
        block_stat.progressive_mean_squared(av2);   //progressive mean squared
        block_stat.progressive_error();             //progressive standard deviation

        //PRINTING OUPUT
        cout<<"Printing output in file "<<file_out<<endl<<endl;
        block_stat.print_output(file_out);          //printing the outpu: progressive mean and progressive error 

return 0;
}


