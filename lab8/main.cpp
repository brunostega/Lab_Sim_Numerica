#include <iostream>
#include <math.h>
#include <cmath>
#include "main.h"
#include <complex>

using namespace std;


int main (int argc, char *argv[]){

	Initialization();

	metro.Find_delta();   //checks if the acceptance in metropolis is more or less 50%
	
	ofstream x_out;       //output files for all positions
	x_out.open(file_out_x);
	cout<<file_out_x<<endl;
	cout<<endl<<"Metropolis alghoritm for all probability distributions using"<<endl;
	cout<<"n_step="<<N<<endl<<endl;

	metro.setX_start(x);
	for(int i=0;i<N_blocks;i++){
		if(i%10==0){
			cout<<"numero di blocco="<<i<<endl;
		}
		for(int j=0;j<L;j++){
			x_appo=metro.MR_alghoritm();
			x_out<<x_appo.getX()<<endl;
			ave_H[i]+=Hamiltonian(x_appo.getX(), mu, sigma);
		}
		ave_H[i]/=L;
	}
	x_out.close();

	//BLOCK STATISTICS AND PRINTING OUTPUTS
	cout<<endl<<"Block Statistic using"<<endl;
	cout<<"N_blocks="<<N_blocks<<endl;
	cout<<"N_steps per block="<<L<<endl;
	BlockPrint();

return 0;
}




void BlockPrint(){
	if(BoolOptimization==1){ blockstat.block_stat_all(file_out_H_Opt, ave_H);}
	else{blockstat.block_stat_all(file_out_H, ave_H);
}
}

void Initialization(){
	ifstream ReadInput;
	ReadInput.open("input.dat");
	
  	ReadInput >> BoolOptimization;
	x=posizione(X);


	//brutal equilibration
	for(int j=0;j<1000;j++)
	{
		x_appo=metro.MR_alghoritm();
	}

	if(BoolOptimization==1)
	{
		cout<<"The program searches for the best values of mu and sigma with a grid search"<<endl;
		file_out_x="x_out_Optimized.out";
		
		Optimization();
		cout<<"Optimization of parameters. Best values are:"<<endl;
		cout<<"mu="<<mu<< "  sigma="<<sigma<<endl;
		WF.reset(mu,sigma);               //reset wave function class with corrected parameters
		metro.SetProbabillity( & WF);     //changing probability function to sample

	}
	else
	{
		cout<<"The program uses standard values for mu and sigma: mu="<<mu<<"  sigma="<<sigma<<endl;
		file_out_x="x_out.out";	
	}
}

double Hamiltonian(double pos, double m_mu, double m_sigma){

	double exp2=exp( - pow(pos-m_mu,2) / (2* m_sigma*m_sigma ) );
	double exp1=exp( - pow(pos+m_mu,2) / (2* m_sigma*m_sigma ) );
	double V=pow(pos,4)-5./2* pow(pos,2);
	double T= (1./(2. * pow(m_sigma,2))) * (1. - (pow(pos,2) + pow(m_mu,2))/pow(m_sigma,2) - (2 * pos * m_mu/pow(m_sigma,2)) * (exp1 - exp2)/(exp1+ exp2));

	return  T+V;

}

void Optimization(){
	ofstream Firstgrid, Secondgrid;
	Firstgrid.open("output.Firstgrid.dat");
	double  mu_appo,sigma_appo, mu_min=0, sigma_min=0;
	double H_appo_min=10000;
	Wave_func WF_appo(mu_min, sigma_min);
	Metropolis metro_appo(N, & WF_appo, & T_unif1);          //defining the metropolis class with the probability density to sample and the transition probability functions
	posizione x_test;

	for(int j=0;j<1000;j++)
	{
		x_test=metro.MR_alghoritm();
	}

	for(int i=0;i<step_grid1;i++)
	{
		mu_appo=mu_start+(i+1)*delta_grid;            //redefine mu
		for(int j=0;j<step_grid1;j++)
		{
			sigma_appo=sigma_start+(j+1)*delta_grid;  //redefine sigma
			WF_appo.reset(mu_appo,sigma_appo);        //reset Wave Function
			posizione x_begin=posizione(0.01);        //reset starting position
			metro_appo.setX_start(x_begin);            

			metro_appo.SetProbabillity(& WF_appo);    //reset metropolis

			for(int j=0;j<1000;j++)
			{
				x_test=metro.MR_alghoritm();
			}

			double H_appo_new=0;                      //reset energy
			for(int k=0;k<nstep_opt;k++)
			{
				posizione x_test=metro_appo.MR_alghoritm();   //metropolis alghoritm
				x_begin=x_test;
				H_appo_new+=Hamiltonian(x_test.getX(), mu_appo, sigma_appo);        //estimate energy
			}
			H_appo_new/=nstep_opt;
			Firstgrid<< mu_appo << ',' << sigma_appo << ',' << H_appo_new <<endl;

			if(H_appo_new<H_appo_min)
			{
				cout<<H_appo_new<< "  mu="<<mu_appo<< "  sigma="<<sigma_appo<<endl;
				mu_min=mu_appo;
				sigma_min=sigma_appo;
				H_appo_min=H_appo_new;
			}
		}
	}

	Firstgrid.close();

	//second refined search grid
	double mu_start2=mu_min-0.04;
	double sigma_start2=sigma_min-0.04;
	double delta_grid2=0.002;
	H_appo_min=10000;

	Secondgrid.open("output.Secondgrid.dat");
	for(int i=0;i<step_grid2;i++)
	{
		mu_appo=mu_start2+(i+1)*delta_grid2;            //redefine mu
		for(int j=0;j<step_grid2;j++)
		{
			sigma_appo=sigma_start2+(j+1)*delta_grid2;  //redefine sigma
			WF_appo.reset(mu_appo,sigma_appo);          //reset Wave Function
			posizione x_begin=posizione(0.01);          //reset starting position

			metro_appo.setX_start(x_begin);
			metro_appo.SetProbabillity(& WF_appo);      //reset metropolis

			for(int j=0;j<1000;j++)
			{
				x_test=metro.MR_alghoritm();
				x_begin=x_test;
			}
			double H_appo_new=0;                        //reset energy
			for(int k=0;k<nstep_opt;k++)
			{
				posizione x_test=metro_appo.MR_alghoritm();   //metropolis alghoritm
				x_begin=x_test;
				H_appo_new+=Hamiltonian(x_test.getX(), mu_appo, sigma_appo);        //estimate energy
			}
			H_appo_new/=nstep_opt;

			Secondgrid<< mu_appo << ',' << sigma_appo << ',' << H_appo_new <<endl;

			if(H_appo_new<H_appo_min)
			{
				cout<<H_appo_new<< "  mu="<<mu_appo<< "  sigma="<<sigma_appo<<endl;
				mu_min=mu_appo;
				sigma_min=sigma_appo;
				H_appo_min=H_appo_new;
			}
		}
	}
	Secondgrid.close();

	mu=mu_min;        //assign the best values of mu and sigma
	sigma=sigma_min;
}


