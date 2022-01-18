#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../../random.h"
#include "../../block_stat.h"
#include "lib.h"
#include <string>

using namespace std;

int main (int argc, char *argv[]){

//INITIALISING SOME VARIABLES

int M=1000000;        //# of pseudo-random number
int N=100;            //# of blocks
int L=M/N;            //# for block
Random rand;          //random generator 
rand.Initialize();    //initializing the random variable

BlockStat block_stat(N,M);    //initializing class for block statistics

double ave[N]{};        //array of mean in each block
double av2[N]{};        //array of mean squared in each block


///////////////////////////////////////////////////////////////////////////////////////////////

//1 STEP: MEAN
//testing the pseudo random generator estimating (with uncertainties)
//the mean value using block statistic


//EVALUATE AVERAGE IN EACH BLOCK
cout<<"Step 1:Evaluating mean value"<<endl;

double sum=0;    
for(int i=0; i<N; i++){
	for(int j=0;j<L;j++){
		sum+=rand.Rannyu();            //cumulative sum in each block
	}
	ave[i]=sum/L;                      //estimating mean for each block
	av2[i]=pow(ave[i],2);              //estimating mean squared for each block
	sum=0;	
}


//BLOCK STATISTICS
cout<<"Doing block statistics"<<endl;

block_stat.progressive_mean(ave);           //progressive mean
block_stat.progressive_mean_squared(av2);   //progressive mean squared
block_stat.progressive_error();             //progressive standard deviation

//PRINTING THE RESULT
const string file_out="stima_r.csv";        //output file
cout<<"printint the result in"<<file_out<<endl;
block_stat.print_output(file_out);          //printing the outpu: progressive mean and progressive error 


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//2 STEP: VARIANCE
//testing the pseudo random generator estimating (with uncertainties)
//the variance using block statistic

//EVALUATE AVERAGE IN EACH BLOCK
cout<<"Step 2:Evaluating the variance"<<endl;

sum=0;   //this is going to be the variance
for(int i=0; i<N; i++){
	for(int j=0;j<L;j++){
		sum+=pow((rand.Rannyu()-0.5),2);             //cumulative sum in each block	
	}
	ave[i]=sum/L;                                    //estimating mean for each block
	av2[i]=pow(ave[i],2);                            //estimating mean squared for each block
	sum=0;
}

//BLOCK STATISTICS
cout<<"Doing block statistics"<<endl;

block_stat.progressive_mean(ave);           //progressive mean
block_stat.progressive_mean_squared(av2);   //progressive mean squared
block_stat.progressive_error();             //progressive standard deviation

//PRINTING THE RESULT
const string file_out2="stima_sigma2.csv";        //output file
cout<<"printint the result in"<<file_out2<<endl;
block_stat.print_output(file_out2);          //printing the outpu: progressive mean and progressive error 



/////////////////////////////////////////////////////////////////////////////////////////

//3 STEP: CHI-SQUARED DISTRIBUTION
//chi squared test
//n=10000 throws in [0,1] divided in 100 intervals
//we expect the variance of the number of throws to be n/M (equal to the mean value since the distribution is a poisson one)
//so the chi-squared should be equal to circa 100

cout<<"Step 3: chi-squared test"<<endl;

M=100;           //# of subintervals in [0,1]
int n=10000;         //# of throws
double subint[M]{};   //vector of subintervals
double chi2[M]{};     //chi squared array
int count[M]{};      //counter of number of throws in each subinterval already initialized to 0

for(int i=0;i<M;i++){      //check on the counter
	if(count[i]!=0){
		cout<<"errore, il contatore non è inizializzato a zero"<<endl;
	}
}

for(int i=0;i<M;i++){       //divide [0,1] in M subintervals
	subint[i]=i*double(1)/M;
}

double appo=0;    //support variable

for(int i=0;i<M;i++){  //running the "experiment" (10000 throws) M times
	for(int j=0;j<n;j++){
		appo=rand.Rannyu();
		for(int k=0;k<M;k++){      //count throws in each subinterval
			if((appo>=subint[k])&(appo<=subint[k]+double(1)/M)){
				count[k]+=1;
			}
		}
	}
	for(int l=0;l<M;l++){      //calculate chi for each run
		chi2[i]+=pow((count[l]-n/M),2)/(n/M);
	}
	for(int s=0;s<M;s++){      //ri-setting count=0  
		count[s]=0;
	}
}

//PRINTING OUTPUT

const string file_out3="chi2.csv"; //output file  
cout<<"printing the result in: "<<file_out3<<endl;
ofstream output;
output.open("chi2.csv");           //outputting chi2
if(output.is_open()){
	for(int i=0;i<N;i++){
		output<<chi2[i]<<endl;
	}
}
else{
	cerr<<"PROBEL, couldn't open the file"<<endl;
}
output.close();


return 0;
}
