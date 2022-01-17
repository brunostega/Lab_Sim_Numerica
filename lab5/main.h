#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <cmath>
#include "../random.h"
#include "../block_stat.h"
#include "lib.h"
#include <complex>

//INITIALIZING SOME VARIBALES

const int N_throws=5000000;       //# throws
const int N_blocks=100;           //# blocks
const int L=N_throws/N_blocks;    //# throws per block
const int N_equilib=250;         //# equilibration steps      
double X=20.,Y=3.,Z=1.1;         //starting position
posizione x(X,Y,Z);               //starting position class

Hydrogen100 H100;          //defining the probability function of 1s hydrogen
Hydrogen210 H210;          //defining the probability function of 2p hydrogen

double l1=1.2, l2=2.8, l3=0.7, l4=1.6; //leghts of the metropolis steps 
T_uniform T_unif1(l1);    //Transition probabiity functions with a defined step to get the 50% rule
T_uniform T_unif2(l2);    //""
T_gaussian T_gauss1(l3);  //""
T_gaussian T_gauss2(l4);  //""

Metropolis MR100(N_equilib,& H100, & T_unif1, x, "MR100");          //defining the metropolis class with the probability density to sample and the transition probability function
Metropolis MR210(N_equilib,& H210, & T_unif2, x,"MR210");          //""
Metropolis MR100_gauss(N_equilib,& H100, & T_gauss1, x,"MR100_gauss");   //"" 
Metropolis MR210_gauss(N_equilib,& H210, & T_gauss2, x,"MR210_gauss");   //""

BlockStat blockstat(N_blocks, N_throws);       //initializing class for block statistics

const std::string file_out_x_100="output/x_out_100.out";            		 //output file: position of sampled points of 1s hydrogen density probability 
const std::string file_out_x_210="output/x_out_210.out";            		 //output file: position of sampled points of 2p hydrogen density probability
const std::string file_out_r_mean_100="output/r_mean_100.out";      		 //output file: block statistics for mean radius of 1s orbital
const std::string file_out_r_mean_210="output/r_mean_210.out";      		 //output file: block statistics for mean radius of 2p orbital
const std::string file_out_r_mean_100_gauss="output/r_mean_100_gauss.out";  //output file: block statistics for mean radius of 1s orbital
const std::string file_out_r_mean_210_gauss="output/r_mean_210_gauss.out";  //output file: block statistics for mean radius of 2p orbital

double ave_r_100[N_blocks]={};        //arrays for average 
double ave_r_210[N_blocks]={};        //""
double ave_r_100_gauss[N_blocks]={};  //""
double ave_r_210_gauss[N_blocks]={};  //""

