
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "block_stat.h"
#include "random.h"
#include <string>

using namespace std;

///////
//CLASSE STATISTICA A BLOCCHI

BlockStat :: BlockStat(int N_blocks,  int N_throws){     //costruttore (posso eventualmente dargli una costruzione predefinita ma generalemte dovrò all'interno del codice settare la classe)
	
	_N_blocks=N_blocks;                    //#of blocks
	_N_throws=N_throws;                    //#of throws
	if(_N_blocks==0 or _N_throws==0){cout<<endl<<"Error initialization of BlockStat went wrong: N_blocks or N_throws is 0"<<endl<<endl;}
	if(N_throws%N_blocks==0){              //CHECK: N_throws must be an integer multiple of N_blocks 
		_n=N_throws/N_blocks;              //#of throws in each block
	}
	else{cout<<"error number of throws must be an exact multiple of number of blocks"<<endl;}
	_final_mean=0;
	_final_error=0;

	_x=new int[N_blocks]{};                //progressive number of throws per block  
	_mean_prog=new double[N_blocks]{};     //progressive sum of mean
	_mean2_prog=new double[N_blocks]{};    //progressive sum of squared mean
	_err_prog=new double[N_blocks]{};      //progressive standard deviation

	for(int i=0; i<N_blocks;i++){
        _x[i]=(i+1)*_n;                        //number of throws per block 
    }
}        
 
BlockStat :: ~BlockStat(){   //distruttore (necessario quando nella classe ci sono vettori dinamici che vanno deallocati dalla memoria)
	delete [] _x;            //deleting all dynamical arrys in the class
	delete [] _mean_prog;
	delete [] _mean2_prog;
	delete [] _err_prog;
}        

void BlockStat :: reset(){
	for(int i=0;i<_N_blocks; i++){
		_mean_prog[i]=0;
		_mean2_prog[i]=0;
		_err_prog[i]=0;
	}
}


void BlockStat :: progressive_mean(double* ave){   //construct an array of cumulative mean 

	//check before starting the cumulative statistic that mean array is correctly initialized to zero
	int check=0;
	for(int i=0; i<_N_blocks; i++){
		if(_mean_prog[i]!=0){
			check++;
			_mean_prog[i]=0;
		}
	}
	if(check>0){
		//cout<<"Progressive mean array was not initialized to zero. Already resetted to zero before starting block statistic"<<endl;
		_final_mean=0;	     //resetting to zero olso the final value of the mean
	}


	for(int i=0; i<_N_blocks; i++){
		for(int j=0;j<i+1;j++){
			_mean_prog[i]+=ave[j];    //cumulative sum of the means in each block
		}
		_mean_prog[i]/=(i+1);             //progressive mean
	}

	_final_mean=_mean_prog[_N_blocks-1];

}


void BlockStat :: progressive_mean_squared(double* av2){   //construct an array of cumulative mean squared 

	int check=0;
	for(int i=0; i<_N_blocks; i++){
		if(_mean2_prog[i]!=0){
			check++;
			_mean2_prog[i]=0;
		}
	}
	//if(check>0){cout<<"Progressive mean squared array was not initialized to zero. Already resetted to zero before starting block statistic"<<endl;}


	for(int i=0; i<_N_blocks; i++){
		for(int j=0;j<i+1;j++){
			_mean2_prog[i]+=av2[j];  //cumulative sum of the means squaredin each block
		}
	_mean2_prog[i]/=(i+1);           //progressive mean squared
	}
}


void BlockStat :: progressive_error(){

	//check before starting the cumulative statistic that error array is correctly initialized to zero
	int check=0;
	for(int i=0; i<_N_blocks; i++){
		if(_err_prog[i]!=0){
			check++;
			_err_prog[i]=0;
		}
	}
	if(check>0){
		//cout<<"Progressive error array was not initialized to zero. Already resetted to zero before starting block statistic"<<endl;
		_final_error=0;        //resetting to zero olso the final value of the mean
	}

	for(int i=0; i<_N_blocks; i++){
		_err_prog[i]=this->error(_mean_prog[i],_mean2_prog[i],i);   //statistical uncertenty
	}
	_final_error=_err_prog[_N_blocks-1];
	//cout<<_final_error<<endl;
}


void BlockStat :: print_output(const string file_out){
	    ofstream output;
        output.open(file_out);
        if(output.is_open()){
                for(int i=0;i<_N_blocks;i++){
                        output<<_x[i]<<","<<_mean_prog[i]<<","<<_err_prog[i]<<endl;
                }
        }
        else{
                cerr<<"PROBEL, couldn't open the file"<<endl;
        }
        output.close();
}



double BlockStat :: error(double mean, double mean2, int n){   //returning the statistical uncertenty as 
	if(n==0){
		return 0;
	}
	else if(mean2-mean*mean<0){
		cout<<"error media quadro: "<<mean2-mean*mean<<endl;
		return 0;
	}
	else{
		return pow((mean2-mean*mean)/n,0.5);
	}
}


double BlockStat :: get_final_mean(){

	if(_final_mean==0){       //check: if is 0, it means block statistics wasn't done already
		cout<<"Problem: final mean was not estimated already"<<endl;
	}
	return _final_mean;
}

double BlockStat :: get_final_error(){

	if(_final_error==0){       //check: if is 0, it means block statistics wasn't done already
		cout<<"Problem: final error was not estimated already"<<endl;
	}
	return _final_error;
}


void BlockStat :: block_stat_all(const string file_out, double * ave){    //makes all block statistics: estimate progressive mean, mean squared, error and print in file
	
	double av2[_N_blocks]{};   //construct the array of mean squared
	cout<<endl;
	for(int i=0;i<_N_blocks;i++){
		av2[i]=ave[i]*ave[i];
	}

	this->progressive_mean(ave);
	this->progressive_mean_squared(av2);
	this->progressive_error();
	this->print_output(file_out);

}


/////////////////////////
