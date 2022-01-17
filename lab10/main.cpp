#include "main.h"

using namespace std;


int main (int argc, char *argv[]){

	string shape=argv[1];
	string beta_law=argv[2];
	Travel_Salesman TS(shape, beta_law);
	TS.Sim_Annealing();

return 0;
}


Travel_Salesman::Travel_Salesman(string sh, string bl){   //inizialize class with every variable needed
	
	shape=sh;
	beta_law=bl;
	ifstream ReadInput;
	ReadInput.open("input.dat");
	ReadInput >> N_city;      //# cities
	ReadInput >> N_beta;      //# beta
	ReadInput >> N_step;      //# steps per beta
	ReadInput >> Beta_step;   //beta step
	ReadInput.close();

	beta=vector<double>(N_beta);    //array of temperature
	rnd.Initialize();              //Random generator
	this->City_Initialize();       //defining all cities
	this->Beta_Initialize();

}

Travel_Salesman::~Travel_Salesman(){}

void Travel_Salesman::Sim_Annealing(){
	ofstream out_order, out_L, out_L_equil;
	out_order.open("output.order_"+shape+beta_law+".dat");
	out_L.open("output.L_"+shape+beta_law+".dat");
	out_L_equil.open("output.L_"+shape+beta_law+"equil.dat");

	double L_check=0;
	int N_check=0;
	Individuo Order_new;
	cout<<beta_law<<endl;
	cout<<"N_step="<<N_step<<endl;
	for(int j=0;j<N_beta;j++)
	{
		for(int i=0; i<N_step;i++)
		{
			Order_new=this->Mutation();     //mutate cities order
			Order_new.Fitness();            //estimate loss function
			Order.Fitness();                //estimate loss function
			Boltz_weight=exp( - ( Order_new.GetLoss() - Order.GetLoss() ) * beta[j]);    
			double r=rnd.Rannyu();

			if(1<Boltz_weight)              //metropolis: Boltzmann weight is more then one, change order 
			{
				Order=Order_new;
			}
			else                            //otherwise change with probability r=boltzmann weight
			{
				if(r>Boltz_weight)
				{
					Order=Order;
				}
				else
				{
					Order=Order_new;
				}
			}
			out_L_equil<<j*N_step+i<<","<<Order.GetLoss()<<endl;   //printing out fitness as a function of beta and steps
		}

		if(L_check==Order.GetLoss())    //convergence alghoritm
		{
			N_check++;
			//cout<<N_check<<",  "<<Order.GetLoss()<<endl;
		}
		else
		{
			L_check=Order.GetLoss();
			N_check=0;
		}
		if(N_check==50)
		{
			cout<<"reached convergence at step: "<<j<<" and at temperature: "<<beta[j]<<endl;
			break;
		}
		out_L<<beta[j]<<","<<j<<","<<Order.GetLoss()<<endl;   //printing out fitness as a function of beta and steps
	}

	for(int i=0; i<N_city;i++){out_order<<Order.GetElement(i)<<endl;}  //print city order of the final result
	out_order<<Order.GetElement(0);

	out_order.close();
	out_L.close();
}

void Travel_Salesman::Beta_Initialize(){
	double A=0.08;
	double B=0.1;
	double C=2.3;
	beta[0]=beta0;
	if(beta_law=="A")   //linear change of beta 
	{
		for(int i=1; i<N_beta; i++)
		{
			if(i<=(int)(N_beta*0.2)){beta[i]=beta[i-1]+A;}
			else if(i>(int)(N_beta*0.2) and i<=(int)(N_beta*0.5)){beta[i]=beta[i-1]+B;}
			else{beta[i]=beta[i-1]+C;}
		}
	}
	else if(beta_law=="B")   //quadratic change of beta
	{
		for(int i=1; i<N_beta; i++)
		{
			Beta_step=0.001;
			beta[i]=beta[0]+Beta_step*i*i;
		}
	}
	else if(beta_law=="C")   //cubic change of beta
	{
		for(int i=1; i<N_beta; i++)
		{
			Beta_step=0.001;
			beta[i]=beta[0]+Beta_step*i*i*i;
		}
	}
	else if(beta_law=="D")   //exponential change of beta
	{
		for(int i=1; i<N_beta; i++)
		{
			Beta_step=0.6*log(10)/100;
			beta[i]=beta[0]*exp(Beta_step*i);
		}
	}
	else
	{
		cout<<endl<<"Problem, didn't choose any of the possible values A,B,C,D"<<endl<<endl;
	}
}

void Travel_Salesman::City_Initialize(){    //initialization of random distribution of cities
	
	ofstream cities_pos;
	vector<vector<double> > pos(N_city, vector<double>(2));
	double theta=0,x=0,y=0;
	double x0=0,y0=0;

	//generate random position in a square or on a circle
	cities_pos.open("output."+shape+"_pos.dat");
	rnd.Rannyu();
	rnd.Rannyu();
	for(unsigned int i=0;i<N_city;i++)
	{
		if(shape=="circle")
		{
			theta=rnd.Rannyu()*2*pi_greco;
			x=cos(theta);
			y=sin(theta);
		}
		else if(shape=="square")
		{
			x=rnd.Rannyu();
			y=rnd.Rannyu();
		}
		else
		{
			cout<<"problem! the string must be square or circle. (x,y) is set to zero"<<endl;
			x=0;
			y=0;
		}

		if(i==0)
		{
			x0=x;
			y0=y;
		}

		pos[i][0]=x;
		pos[i][1]=y;

		cities_pos<<x<<","<<y<<endl;
	}
	Order.SetPosition(pos);
	cities_pos<<x0<<","<<y0<<endl;
	cities_pos.close();
}

Individuo Travel_Salesman::Mutation(){
	Individuo Order_new=Order;
	double r=rnd.Rannyu();
	if(r<1./3){	Order_new.Mutation1(& rnd);};
	if(r<2./3 and r>1./3){Order_new.Mutation2(& rnd);}
	if(r>2./3){	Order_new.Mutation3(& rnd);}

	Order_new.Check();
	return Order_new;
}


//CLASS INDIVIDUO
Individuo::Individuo(){
	ifstream ReadInput;
	ReadInput.open("input.dat");
	ReadInput >> N_city;      //# cities
	ReadInput.close();

	Order=vector<int>(N_city);
	pos=vector< vector <double> >(N_city, vector<double>(2));
	L=0;
	exp=1;
	for(unsigned int i=0;i<N_city;i++){Order[i]=i+1;}
}
Individuo::Individuo(unsigned int N){

	N_city=N;      //# cities
	Order=vector<int>(N_city);
	L=0;
	exp=1;
	for(unsigned int i=0;i<N_city;i++){Order[i]=i+1;}
}
Individuo::~Individuo(){}
double Individuo::GetLoss(){return L;}
void Individuo::SetLoss(double loss){L=loss;}
void Individuo::SetElement(int appo, int pos){Order[pos]=appo;}
void Individuo::SetPosition(vector< vector<double> > position){pos=position;}
int Individuo::GetElement(int pos){return Order[pos];}

void Individuo::Fitness(){  //estimate fitness 
	double loss=0;
	for(unsigned int j=0;j<N_city;j++)
	{
		int j2=j+1;
		if(j==N_city-1){j2=0;}
		loss+=pow(  abs(pos[ this->GetElement(j)-1 ][0]-pos[ this->GetElement(j2)-1 ][0])  ,  exp) + pow(  abs(pos[this->GetElement(j)-1 ][1]-pos[this->GetElement(j2)-1 ][1]) , exp);
		L=loss;
	}
}

void Individuo::Check(){  //checks that the city order is correct (no city is reapeated more then once and that the first city is 1)
	int chk[N_city]={};
		int control=0;
		for(unsigned int i=0;i<N_city;i++)
		{
			for(unsigned int j=0;j<N_city;j++)
			{
				if(this->GetElement(j)==i+1){chk[i]+=1;}
			}
		}

		for(unsigned int i=0;i<N_city;i++)
		{
			if(chk[i]!=1){control+=1;}
		}

		if(control!=0 or this->GetElement(0)!=1){
			cout<<"problema!!"<<endl;
			for( unsigned int i=0;i<N_city;i++)
			{
				cout<<this->GetElement(i)<<" ";
			}
			cout<<endl;
		}
}
void Individuo::PrintOrder(){
	cout<< "L="<<left<<setw(10)<<this->GetLoss()<<"Path: ";
	for(unsigned int j=0;j<N_city;j++)
	{
		cout<<this->GetElement(j)<<" ";
	}
	cout<<endl;
}
void Individuo::Mutation1(Random * R){   //swap to numbers

	int appo_i=(int)(R->Rannyu(1,N_city));
	int appo_j=(int)(R->Rannyu(1,N_city));
	while(appo_j==appo_i){appo_j=(int)(R->Rannyu(1,N_city));}
	double appo=0;
	appo=this->GetElement(appo_i);
	this->SetElement(this->GetElement(appo_j),appo_i);
	this->SetElement(appo,appo_j);
}
void Individuo::Mutation2(Random * R){   //tralsate all numbers by some tot
	Individuo copy(*this);
	int pos=0;
	int n=(int)R->Rannyu(1,N_city);
	for(unsigned int i=1;i<N_city;i++)
	{
		if(i+n>=N_city){pos=-N_city+n+i+1;}
		else{pos=i+n;}
		this->SetElement(copy.GetElement(pos),i);
	}
}
void Individuo::Mutation3(Random * R){   //inverte l'ordine di comparsa tra m numeri a partire da start (escluso il primo)
	Individuo copy(*this);
	int pos_iniz=0, pos_fin=0;
	int m=(int)(R->Rannyu(1,N_city));
	int start=(int)(R->Rannyu(1,N_city));

	//cout<<"invert "<<m<<" positions, starting from "<<start<<endl;
	for(int i=0;i<=m;i++)
	{
		if(start+i>=N_city){pos_iniz=-N_city+start+i+1;}
		else{pos_iniz=i+start;}
		if(start+m-i>=N_city){pos_fin=-N_city+start+m-i+1;}
		else{pos_fin=start+m-i;}
		this->SetElement(copy.GetElement(pos_fin),pos_iniz);
	}

}
