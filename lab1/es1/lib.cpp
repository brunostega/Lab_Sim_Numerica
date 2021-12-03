
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "lib.h"
#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include <string>

using namespace std;

float error(float mean, float mean2, int n){  	//Calculate the standard deviation given the mean and the mean squared

	if(n==0){
		return 0;
	}
	else if(mean2-mean*mean<0){
		cout<<"error"<<endl;
		return 0;
	}
	else{
		return pow((mean2-mean*mean)/n,0.5);
	}
};



void cumulative_mean(float* values, float* mean, int N){   	//construct an array of cumulative mean (mean) 
															//from a starting vector of values (values)

	for(int i=0; i<N; i++){
		for(int j=0;j<i+1;j++){
			mean[i]+=values[j];     //cumulative sum of the means in each block
		}
		mean[i]=mean[i]/(i+1);          //cumulative mean
	}
}


int count(const string& file){       //counts the number of elements in a file

        int N=0;
        double appo;
        ifstream in;
        in.open(file);

	if(in.is_open()){
        	in>>appo;
        	while(!in.eof()){
        		N++;
        		in>>appo;
			//cout<<appo<<endl;
		}
	in.close();
        } else cerr << "PROBLEM: Unable to open file" << endl;
        
return N;
};


////////////////////
