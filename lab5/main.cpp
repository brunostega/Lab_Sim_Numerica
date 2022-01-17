#include "main.h"



using namespace std;

int main (int argc, char *argv[]){

	//Experiment
	cout<<"Hydrogen wave function sampling using metropolis alghoritm and block statistics with an equilibration phase"<<endl;
	cout<<"of "<<N_equilib<<"steps"<<endl<<endl;
	if(L<1000){cout<<"the number of steps per block is too low! possible problem when saving position markov chain"<<endl;}

	//Acceptance percetage control
	cout<<endl<<"controlling for every probability function to be sampled if the Metropolis alghoritm satisfy the 50% rule"<<endl<<endl;
	ofstream out;
	out.open("output/50rule.dat");
	out<<l1<<","<<MR100.Find_delta()<<endl;   //checks if the acceptance in metropolis is more or less 50%
	out<<l2<<","<<MR210.Find_delta()<<endl;   
	out<<l3<<","<<MR100_gauss.Find_delta()<<endl;
	out<<l4<<","<<MR210_gauss.Find_delta()<<endl;
	out.close();

	ofstream x_out_100;             //output files for all positions
	ofstream x_out_210;
	x_out_100.open(file_out_x_100);
	x_out_210.open(file_out_x_210);

	cout<<endl<<"Metropolis alghoritm for all probability distributions"<<endl;
	posizione x_appo,x_appo2;
	for(int i=0;i<N_blocks;i++){
		if(i%10==0){
			cout<<"numero di blocco="<<i<<endl;
		}
		for(int j=0;j<L;j++){
			x_appo=MR100.MR_alghoritm();

			ave_r_100[i]+=x_appo.getR();
			ave_r_100_gauss[i]+=MR100_gauss.MR_alghoritm().getR();

			x_appo2=MR210.MR_alghoritm();
			ave_r_210[i]+=x_appo2.getR();
			ave_r_210_gauss[i]+=MR210_gauss.MR_alghoritm().getR();
			if(j%100==0){
				x_out_100<<x_appo.getX()<<","<<x_appo.getY()<<","<<x_appo.getZ()<<endl;
				x_out_210<<x_appo2.getX()<<","<<x_appo2.getY()<<","<<x_appo2.getZ()<<endl;
			}
		}
		ave_r_100[i]/=L;
		ave_r_210[i]/=L;
		ave_r_100_gauss[i]/=L;
		ave_r_210_gauss[i]/=L;
	}
	x_out_100.close();
	x_out_210.close();

	//BLOCK STATISTICS AND PRINTING OUTPUTS
	cout<<endl<<"Block Statistic using"<<endl;
	cout<<"N_blocks="<<N_blocks<<endl;
	cout<<"N_steps per block="<<L<<endl;
    blockstat.block_stat_all(file_out_r_mean_100, ave_r_100);
	blockstat.block_stat_all(file_out_r_mean_210, ave_r_210);
   	blockstat.block_stat_all(file_out_r_mean_100_gauss, ave_r_100_gauss);
	blockstat.block_stat_all(file_out_r_mean_210_gauss, ave_r_210_gauss);

return 0;
}