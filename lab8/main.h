#include "../random.h"
#include "../block_stat.h"
#include "lib.h"
#include <string>
#include <iostream>
#include <fstream>



//block statistics
int N_throws=1000000;               //# throws
const int N_blocks=100;           //# blocks
int L=N_throws/N_blocks;          //# throws per block
double ave_H[N_blocks]={};        //arrays for average 
const std::string file_out_H="H_mean.out";            		//output file:  
const std::string file_out_H_Opt="H_mean_Opt.out";            		//output file:  


//metropolis
int N=1;                 //# of temporal steps
double X=0.1;            //starting position
posizione x(X);          //starting position class
posizione x_appo;

//variables
double sigma=0.2,mu=0.5;
double sigma_start=0.05,mu_start=0.05;
double delta_grid=0.05;
double step_grid1=20;
double step_grid2=40;
bool BoolOptimization;
int nstep_opt=100000;



Wave_func WF(mu,sigma);
T_uniform T_unif1(5);                        //Transition probabiity functions with a defined step to get the 50% rule
BlockStat blockstat(N_blocks, N_throws);     //initializing class for block statistics
Metropolis metro(N, & WF, & T_unif1);        //defining the metropolis class with the probability density to sample and the transition probability functions

std::string file_out_x;    //output file: position of sampled points of 1s hydrogen density probability 


//functions
void Initialization();
double Hamiltonian(double, double, double);
void Optimization();
void BlockPrint();