#include "main.h"

using namespace std;


int main (int argc, char *argv[]){

	string shape=argv[1];
	ofstream order, out_L;
	order.open("order_"+shape+".dat");
	out_L.open("L_"+shape+".dat");
	Gene_Alghoritm gene(shape);
	int generation=gene.Get_NGeneration();
	int N_city=gene.Get_NCity();
	int N_indiv=gene.Get_NIndiv();



	for(int i=0;i<generation;i++)
	{
		gene.New_Generation();
		if(i%10==0)
		{
			cout<<"generation "<<i<<", best fitness="<<gene.GetIndividuo(0).GetLoss()<<endl;
			for(int j=0;j<(int)(N_indiv/2);j++)
			{
				ave_L+=gene.GetIndividuo(j).GetLoss();
			}
			ave_L/=(int)(N_indiv/2);
			out_L<<i<<","<<gene.GetIndividuo(0).GetLoss()<<","<<ave_L<<endl;
			ave_L=0;
		}
	}	
	

	for(int i=0;i<N_city;i++){
		order<<gene.GetIndividuo(0).GetElement(i)<<endl;
	}
	order<<gene.GetIndividuo(0).GetElement(0)<<endl;
	order.close();
	out_L.close();

return 0;
}



//CLASS iNDIVIDUO
Individuo::Individuo(){
	ifstream ReadInput;
	ReadInput.open("input.dat");
	ReadInput >> N_city;      //# cities
	ReadInput.close();

	indiv=vector<int>(N_city);
	L=0;
	exp=1;
	for(unsigned int i=0;i<N_city;i++){indiv[i]=i+1;}
	rnd.Initialize();
}
Individuo::Individuo(unsigned int N){

	N_city=N;      //# cities
	indiv=vector<int>(N_city);
	L=0;
	exp=1;
	for(unsigned int i=0;i<N_city;i++){indiv[i]=i+1;}
	rnd.Initialize();
}
Individuo::~Individuo(){}
void Individuo::PrintIndiv(){
	cout<< "L="<<left<<setw(10)<<this->GetLoss()<<"Path: ";
	for(unsigned int j=0;j<N_city;j++)
	{
		cout<<this->GetElement(j)<<" ";
	}
	cout<<endl;
}
void Individuo::Check(){
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
		for(int i=0;i<N_city;i++)
		{
			cout<<this->GetElement(i)<<" ";
		}
		cout<<endl;
	}
}
double Individuo::GetLoss(){return L;}
void Individuo::SetLoss(double loss){L=loss;}
void Individuo::SetElement(int appo, int pos){indiv[pos]=appo;}
int Individuo::GetElement(int pos){return indiv[pos];}

//CLASS Gene ALghoritm
Gene_Alghoritm::Gene_Alghoritm(string shape){     //constructor

	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> N_city;        //# cities
	ReadInput >> N_indiv;       //#individuals
	ReadInput >> N_generation;  //#generations
	ReadInput >> P_c;           //crossing probability
	ReadInput >> P_m;           //mutation probability
	ReadInput >> fit_exp;       //selection exponent
	ReadInput >> exp;           //exponent of loss function
	ReadInput >> cross;          //1 or 2. uses different crossover functions
	ReadInput >> grid;          //1 or 2. uses different crossover functions


	cout<<"genetic alghoritm to solve salesman problem"<<endl;
	cout<<"# cities="<<N_city<<", casually distributed on a "<<shape<<endl;
	cout<<"# individuals="<<N_indiv<<endl;
	cout<<"usng crossover "<<cross<<endl;

	rnd.Initialize();      //initialize Random class

	//construct all vectors
	population.resize(N_indiv);
	population_start.resize(N_indiv);

	pos=vector< vector <double> >(N_city, vector<double>(2));
	pos_start=vector< vector <double> >(N_city, vector<double>(2));

	this->City_initialize(shape);
	if(grid==1)
	{
	this->GridSearch(shape);

	this->Reset_Cities();
	this->fitness();
	this->Sort();
	}

	else if(grid==0)
	{
		this->fitness();
		this->Sort();
	}

}

Gene_Alghoritm::~Gene_Alghoritm(){}               //deconstructor


void Gene_Alghoritm::GridSearch(string shape){
	ofstream out;
	out.open("L_"+shape+"_grid.dat");
	double dpm=0.003;
	double dpc=0.02;
	int n_grid_step=10;
	double pm_start=0.07;
	double pc_start=0.5;
	double pm_test, pc_test;


	for(int j=0;j<n_grid_step;j++)
	{
		pm_test=pm_start+j*dpm;   //changes P_m
		this->SetPm(pm_test);     //Set P_m
		for(int k=0;k<n_grid_step;k++)
		{
			pc_test=pc_start+k*dpc;  //Changes P_c
			this->SetPc(pc_test);    //Set P_c

			//reset cities
			this->Reset_Cities();  
			//this->City_initialize(shape);  
			this->fitness();
			this->Sort();
			//this->PrintPopulation();

			for(int i=0;i<this->Get_NGeneration();i++)
			{
				this->New_Generation();
			}	
			out<<pm_test<<","<<pc_test<<","<<this->GetIndividuo(0).GetLoss()<<endl;
			cout<<pm_test<<","<<pc_test<<","<<this->GetIndividuo(0).GetLoss()<<endl;
		}

	}
	out.close();
}
void Gene_Alghoritm::New_Generation(){            //create a new generation. Estimate loss function and reorder the individuals

	if(this->GetCross()==1){	this->crossover1();}
	else if(this->GetCross()==2){	this->crossover2();}
	else
	{
		this->crossover1();	
	}
	this->mutation();
	this->fitness();
	this->Sort();
}
Individuo Gene_Alghoritm::GetIndividuo(int i){return population[i];}  
int Gene_Alghoritm::Get_NGeneration(){return N_generation;}
int Gene_Alghoritm::Get_NCity(){return N_city;}
int Gene_Alghoritm::GetExp(){return exp;}
int Gene_Alghoritm::GetCross(){return cross;}
int Gene_Alghoritm::Get_NIndiv(){return N_indiv;}
int Gene_Alghoritm::GetPos(int i, int j){return pos[i][j];}
void Gene_Alghoritm::SetIndividuo(Individuo copy, int i){
	population[i]=copy;
	population[i].SetLoss(copy.GetLoss());
}


void Gene_Alghoritm::SetPm(double pm){P_m=pm;}
void Gene_Alghoritm::SetPc(double pc){P_c=pc;}


void Gene_Alghoritm::PrintPopulation(){
	for(int i=0;i<N_indiv;i++){population[i].PrintIndiv();}	
	cout<<endl;
}
void Gene_Alghoritm::City_initialize(string shape){    //initialization of random distribution of cities
	
	ofstream cities_pos;
	double theta=0,x=0,y=0;
	double x0=0,y0=0;

	//generate random position in a square or on a circle
	cities_pos.open(shape+"_pos.dat");
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
		pos_start[i][0]=x;
		pos_start[i][1]=y;

		cities_pos<<x<<","<<y<<endl;
	}
	cities_pos<<x0<<","<<y0<<endl;
	cities_pos.close();
	
	//randomize all vectors
	for(unsigned int i=0;i<N_indiv;i++)
	{
		for(unsigned int j=0;j<N_city;j++)
		{
			int appo_i=(int)(rnd.Rannyu(1,N_city));
			int appo_j=(int)(rnd.Rannyu(1,N_city));
			int  appo=0;
			appo=population[i].GetElement(appo_j);
			population[i].SetElement(population[i].GetElement(appo_i),appo_j);
			population[i].SetElement(appo,appo_i);
		}
		population_start[i]=population[i];

		this->Check(population[i]);
		this->Check(population_start[i]);
	}
}

void Gene_Alghoritm::Reset_Cities(){

	for(unsigned int i=0;i<N_city;i++)
	{
		pos[i][0]=pos_start[i][0];
		pos[i][1]=pos_start[i][1];

	}

	for(unsigned int i=0;i<N_indiv;i++)
	{
		population[i]=population_start[i];
		this->Check(population[i]);
	}
}

void Gene_Alghoritm::Check(Individuo V){  //controls one individual is in correct in form
	int chk[N_city]={};
	int control=0;
	for(unsigned int i=0;i<N_city;i++)
	{
		for(unsigned int j=0;j<N_city;j++)
		{
			if(V.GetElement(j)==i+1){chk[i]+=1;}
		}
	}

	for(unsigned int i=0;i<N_city;i++)
	{
		if(chk[i]!=1){control+=1;}
	}

	if(control!=0 or V.GetElement(0)!=1){
		cout<<"problema!!"<<endl;
		for(int i=0;i<N_city;i++)
		{
			cout<<V.GetElement(i)<<" ";
		}
		cout<<endl;
	}
}
void Gene_Alghoritm::Check_all(){for(unsigned int i=0;i<N_indiv;i++){population[i].Check();} //Check if all individuals are correct
}
void Gene_Alghoritm::Sort(){      //sort individuals based on loss function
	double min=0;
	int indice=0;
	double appo=0, appo_int=0;
	Individuo appo_copy(N_city);

	for(unsigned int i=0;i<N_indiv;i++)
	{
		indice=i;
		min=population[i].GetLoss();
		for(int j=i;j<N_indiv;j++)
		{
			if(population[j].GetLoss()<min)
			{
				indice=j;
				min=population[j].GetLoss();
			}
		}
		appo_copy=population[indice];
		this->SetIndividuo(population[i], indice);
		this->SetIndividuo(appo_copy, i);
	}
	//Sortcontrol
	for(unsigned int i=0;i<N_indiv-1;i++)
	{
		if(population[i].GetLoss()>population[i+1].GetLoss()){cout<<"Sorting problem"<<endl;}
	}
}
void Gene_Alghoritm::fitness(){   //estimate loss function for every individual

	for(unsigned int i=0;i<N_indiv;i++)
	{
		double loss=0;
		for(unsigned int j=0;j<N_city;j++)
		{
			int j2=j+1;
			if(j==N_city-1){j2=0;}
			loss+=pow(  abs(pos[ population[i].GetElement(j)-1 ][0]-pos[ population[i].GetElement(j2)-1 ][0])  ,  exp) + pow(  abs(pos[population[i].GetElement(j)-1 ][1]-pos[population[i].GetElement(j2)-1 ][1]) , exp);
			population[i].SetLoss(loss);
		}
		loss=0;
	}
}
int Gene_Alghoritm::select(){     //choosing an individual with weighted probability based on loss function
	return (int)( (N_indiv) * pow(rnd.Rannyu(),fit_exp)  );
}
void Gene_Alghoritm::crossover1(){

	int mother=0, father=0,start=0;
	vector<int> appo_M(N_city);
	vector<int> appo_F(N_city);
	int min_M=N_city;
	int min_F=N_city;
	vector<int> order_mother(N_city);
	vector<int> order_father(N_city);
	vector<Individuo> copy_population(population);

	for(unsigned int i=0;i<N_indiv;i+=2)
	{
		mother=this->select();       //chooses two individuals based on their goodness
		father=this->select();
		while(mother==father){father=this->select();}

		if(rnd.Rannyu()<P_c)
		{
			start=(int)(rnd.Rannyu(2, N_city-2));        //set the number of genes that will remain in the two individuals

			int i_m=start;
			int i_f=start;
			for(int j=0;j<N_city;j++)             //search the position in order 
			{
				for(int k=start;k<N_city;k++)
				{
					if(population[mother].GetElement(j)==population[father].GetElement(k))
					{
						order_mother[i_m]=population[mother].GetElement(j);    //order in wich fathers elements are dispoused in mother
						i_m++;
					}
					if(population[father].GetElement(j)==population[mother].GetElement(k))
					{
						order_father[i_f]=population[father].GetElement(j);    //order in wich mothers elements are dispoused in father
						i_f++;
					}
				}
			}
			
			for(int j=0;j<N_city;j++)
			{
				if(j<start)
				{
					copy_population[i].SetElement(population[mother].GetElement(j), j);
					copy_population[i+1].SetElement(population[father].GetElement(j), j);
				}
				else
				{
					copy_population[i].SetElement(order_father[j], j);
					copy_population[i+1].SetElement(order_mother[j], j);
				}
			}
		}

		else
		{
			copy_population[i]=population[mother];
			copy_population[i+1]=population[father];
		}
	}
	population=copy_population;
	this->Check_all();
}
void Gene_Alghoritm::crossover2(){

	int mother=0, father=0,gene=0;
	vector<int> appo_M(N_city);
	vector<int> appo_F(N_city);
	int min_M=N_city;
	int min_F=N_city;
	vector<Individuo> copy_population(population);

	for(unsigned int i=0;i<N_indiv;i+=2)
	{
		mother=this->select();       //chooses two individuals based on their goodness
		father=this->select();       
		while(mother==father){father=this->select();}
		//cout<<"mother="<<mother<<", father="<<father<<endl;

		//population[mother].PrintIndiv();
		//population[father].PrintIndiv();

		if(rnd.Rannyu()<P_c)
		{
			gene=(int)(rnd.Rannyu(2, N_city-2));   //set the number of genes that will remain in the two individuals
			for(int j=gene;j<N_city;j++)             //search the position in order 
			{
				int pos_F=0,pos_M=0;
				for(int k=gene;k<N_city;k++)
				{
					if(population[mother].GetElement(k)<=min_M and population[mother].GetElement(k)>population[mother].GetElement(appo_M[j-1]))
					{
						pos_M=k;
						min_M=population[mother].GetElement(k);
					}
					if(population[father].GetElement(k)<=min_F and population[father].GetElement(k)>population[father].GetElement(appo_F[j-1]))
					{
						pos_F=k;
						min_F=population[father].GetElement(k);
					}
				}	
				appo_M[j]=pos_M;	//saves the position of ascendent numbers
				appo_F[j]=pos_F;	
				if(pos_F==0 or pos_M==0){
					cout<<"prob in cicle "<<i+1<<endl;
					population[mother].PrintIndiv();
				}
				min_M=N_city;
				min_F=N_city;
			}

			for(unsigned int j=0;j<gene;j++)   //copy first genes
			{
				copy_population[i].SetElement(population[mother].GetElement(j),j);
				copy_population[i+1].SetElement(population[father].GetElement(j),j);
			}

			for(unsigned int j=gene;j<N_city;j++)   //copy first genes
			{
				copy_population[i].SetElement(population[mother].GetElement(appo_M[j]),appo_F[j]);
				copy_population[i+1].SetElement(population[father].GetElement(appo_F[j]),appo_M[j]);
			}

			for(int j=0;j<N_city;j++) //reset the arrays to zero
			{
				appo_F[j]=0;
				appo_M[j]=0;
			}
		}

		else
		{
			copy_population[i]=population[mother];
			copy_population[i+1]=population[father];
		}
	//cout<<gene<<endl;
	//copy_population[i].PrintIndiv();
	//copy_population[i+1].PrintIndiv();
	//cout<<endl;

	}

	population=copy_population;
	this->Check_all();

}
void Gene_Alghoritm::mutation(){
	double r;
	for(unsigned int i=0;i<N_indiv;i++)
	{
		/*
		if(rnd.Rannyu()<P_m)
		{
			r=rnd.Rannyu();
			if(r<1./3){	mutation1(population[i]);};
			if(r<2./3 and r>1./3){mutation2(population[i]);}
			if(r>2./3){	mutation2(population[i]);}
		}
		*/
		if(rnd.Rannyu()<P_m){mutation1(population[i]);}
		if(rnd.Rannyu()<P_m){mutation2(population[i]);}
		if(rnd.Rannyu()<P_m){mutation3(population[i]);}
		if(rnd.Rannyu()<P_m){mutation4(population[i]);}

		
	}
	this->Check_all();
}
void Gene_Alghoritm::mutation1(Individuo &V){   //swap to numbers
	
	int appo_i=(int)(rnd.Rannyu(1,N_city));
	int appo_j=(int)(rnd.Rannyu(1,N_city));
	while(appo_j==appo_i){appo_j=(int)(rnd.Rannyu(1,N_city));}
	double appo=0;
	appo=V.GetElement(appo_i);
	V.SetElement(V.GetElement(appo_j),appo_i);
	V.SetElement(appo,appo_j);
}
void Gene_Alghoritm::mutation2(Individuo &V){   //tralsate all numbers by some tot
	Individuo copy(V);
	int pos=0;

	int n=(int)rnd.Rannyu(1,N_city);
	for(unsigned int i=1;i<N_city;i++)
	{
		if(i+n>=N_city){pos=-N_city+n+i+1;}
		else{pos=i+n;}
		V.SetElement(copy.GetElement(pos),i);
	}
}
void Gene_Alghoritm::mutation3(Individuo & V){   //inverte l'ordine di comparsa tra m numeri a partire da start (escluso il primo)
	Individuo copy(V);
	int pos_iniz=0, pos_fin=0;
	int m=(int)(rnd.Rannyu(1,N_city));
	int start=(int)(rnd.Rannyu(1,N_city));

	//cout<<"invert "<<m<<" positions, starting from "<<start<<endl;
	for(int i=0;i<=m;i++)
	{
		if(start+i>=N_city){pos_iniz=-N_city+start+i+1;}
		else{pos_iniz=i+start;}
		if(start+m-i>=N_city){pos_fin=-N_city+start+m-i+1;}
		else{pos_fin=start+m-i;}
		V.SetElement(copy.GetElement(pos_fin),pos_iniz);
	}

}

void Gene_Alghoritm::mutation4(Individuo &V){   //swap m (m<N/4) numbers in the first half with m numbers in the second half
	
	int appo_i=(int)(rnd.Rannyu(1,N_city/4));
	int appo_j=(int)(rnd.Rannyu(N_city/2,N_city*3/4));
	int m=(int)(rnd.Rannyu(1, N_city/4));
	double appo=0;
	for(int i=0; i<m; i++)
	{
		appo=V.GetElement(appo_i);
		V.SetElement(V.GetElement(appo_j),appo_i);
		V.SetElement(appo,appo_j);
		appo_i++;
		appo_j++;
	}
}

void Gene_Alghoritm::reset_L(){
	for(int i=0;i<N_indiv;i++){population[i].SetLoss(0);}
}


