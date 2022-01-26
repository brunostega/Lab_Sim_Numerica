#include "main.h"

using namespace std;

string path="../../";


int main (int argc, char *argv[]){

	int size, rank;
	int A=0; 
	int N_migration=30;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double tstart = MPI_Wtime();
	Gene_Alghoritm gene(rank);
	int generation=gene.Get_NGeneration();
	int N_city=gene.Get_NCity();
	int N_indiv=gene.Get_NIndiv();

	
//make sure every node has the same cities starting position/////
///////////////

	double * cities_pos_x=new double[N_city]{};
	double * cities_pos_y=new double[N_city]{};
	if(rank==0)
	{
		for(int i=0;i<N_city;i++)
		{
			cities_pos_x[i]=gene.GetPos_Start(i,0);
			cities_pos_y[i]=gene.GetPos_Start(i,1);
		}
	}

	MPI_Bcast(cities_pos_x,N_city,MPI_DOUBLE,0, MPI_COMM_WORLD);	 //Bcasting the starting position to every node
	MPI_Bcast(cities_pos_y,N_city,MPI_DOUBLE,0, MPI_COMM_WORLD);	


	for(int i=0;i<N_city;i++)   //setting position to every gene class of the different node
	{
		gene.SetPos_Start(i,0,cities_pos_x[i]);
		gene.SetPos_Start(i,1,cities_pos_y[i]);
	}
	gene.Print_Cities();   //printing in output every starting position for different ranks 


	gene.fitness();
	gene.Sort();

	delete [] cities_pos_x;
	delete [] cities_pos_y;

/////////////////

	for(int i=0;i<generation;i++)   //genetic alghoritm
	{
		gene.New_Generation();  //new generation: crossover+ mutation+ fitness+ order

		if(i%N_migration==0)    //mpi exchange of best paths
		{

			if(rank==0){
				A=int(gene.rnd.Rannyu(1,size));    //choose a integer between 1 and size-1 (number of cores)
			}

			MPI_Bcast(&A,1,MPI_INTEGER,0, MPI_COMM_WORLD);	 


			MPI_Status stato;
			int BestOrderA [N_city]{};	
			int BestOrderB [N_city]{};	


			for(int j=0; j<size; j++)
			{
				int isend=j+A;    //send indice
				int irecive=j-A;  //recive indice

				if(isend>size-1){isend=isend-size;}   //periodic boudary condition
				if(irecive<0){irecive=irecive+size;}  //periodic boudary condition

				vector <MPI_Request> reqs_send=vector <MPI_Request> (size);
				vector <MPI_Request> reqs_recv=vector <MPI_Request> (size);
				vector <MPI_Status> stats=vector <MPI_Status> (size);
				vector <int> itags= vector <int>(size);

				MPI_Status statoA,statoB;

				for(int i=0;i<size;i++){
					itags[i]=i+1;
				}

				if(rank==irecive)  //irecive sends to j
				{

					for(int i=0; i<N_city;i++){ BestOrderB[i]=gene.GetIndividuo(0).GetElement(i);}
					MPI_Send(&BestOrderB, N_city, MPI_INTEGER, j, itags[irecive], MPI_COMM_WORLD);

				}
				if(rank==j)     //j sends A to isend and recives B from irecive
				{

					MPI_Recv(&BestOrderB, N_city, MPI_INTEGER, irecive, itags[irecive], MPI_COMM_WORLD, &statoB);

					Individuo appo(N_city);
					for(int s=0; s<N_city;s++){appo.SetElement(BestOrderB[s],s);}
					gene.SetIndividuo(appo,0); //setting first individuo
				}

				if(rank==j)     //j sends A to isend and recives B from irecive
				{
					
					for(int i=0; i<N_city;i++){ BestOrderA[i]=gene.GetIndividuo(0).GetElement(i);}
					MPI_Send(&BestOrderA, N_city, MPI_INTEGER, isend, itags[isend], MPI_COMM_WORLD);	
				}

				if(rank==isend) //isend recives from j
				{
					MPI_Recv(&BestOrderA, N_city, MPI_INTEGER, j, itags[isend], MPI_COMM_WORLD, &statoA);
					Individuo appo(N_city);
					for(int s=0; s<N_city;s++){appo.SetElement(BestOrderA[s],s);}
					gene.SetIndividuo(appo,0);
				}

			}
			gene.fitness();
			gene.Sort();
		}


		if(i%10==0)
		{	
			gene.Out_L(i);

		}

	}	
	gene.Print_Order();

	double tend = MPI_Wtime();
	double dt = tend - tstart;
	ofstream  time;
	time.open("output."+gene.shape+".time.dat", ios::app);

	if (time.is_open()){time<<rank<<","<<dt<<endl;} 
	else cerr << "PROBLEM: Unable to open random.out" << endl;
	time.close();

	MPI_Finalize();

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
	rnd.Initialize(path);
}
Individuo::Individuo(unsigned int N){

	N_city=N;      //# cities
	indiv=vector<int>(N_city);
	L=0;
	exp=1;
	for(unsigned int i=0;i<N_city;i++){indiv[i]=i+1;}
	rnd.Initialize(path);
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
void Individuo::SetElement(int appo, int pos){
	indiv[pos]=appo;
}
int Individuo::GetElement(int pos){return indiv[pos];}

//CLASS Gene ALghoritm
Gene_Alghoritm::Gene_Alghoritm(int Rank){     //constructor

	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> N_city;        //#cities
	ReadInput >> N_indiv;       //#individuals
	ReadInput >> N_generation;  //#generations
	ReadInput >> P_c;           //crossing probability
	ReadInput >> P_m;           //mutation probability
	ReadInput >> fit_exp;       //selection exponent
	ReadInput >> exp;           //exponent of loss function
	ReadInput >> cross;          //1 or 2. uses different crossover functions
	ReadInput >>shape;          //shape for cities distribution


	cout<<"genetic alghoritm to solve salesman problem"<<endl;
	cout<<"# cities="<<N_city<<", casually distributed on a "<<shape<<endl;
	cout<<"# individuals="<<N_indiv<<endl;
	cout<<"usng crossover "<<cross<<endl;

	
	rank=Rank;   //mpi rank
	rnd.Initialize_chooseSeed(path, rank+1);      //initialize Random class


	//construct all vectors
	population.resize(N_indiv);
	population_start.resize(N_indiv);

	pos=vector< vector <double> >(N_city, vector<double>(2));
	pos_start=vector< vector <double> >(N_city, vector<double>(2));

	this->City_initialize();
	this->fitness();
	this->Sort();

}

Gene_Alghoritm::~Gene_Alghoritm(){}               //deconstructor



void Gene_Alghoritm::Out_L(int i){
	ofstream  out_L;
	double ave_L=0;
	string whichrank=to_string(rank+1);
	out_L.open("output."+shape+".rank"+whichrank+".L.dat",ios::app);

	if (out_L.is_open()){} 
	else cerr << "PROBLEM: Unable to open random.out" << endl;
	for(int j=0;j<(int)(N_indiv/2);j++)
	{
		ave_L+=this->GetIndividuo(j).GetLoss();
	}
	ave_L/=(int)(N_indiv/2);
	out_L<<i<<","<<this->GetIndividuo(0).GetLoss()<<","<<ave_L<<endl;

	out_L.close();
}

void Gene_Alghoritm::Print_Order(){
	ofstream order;

	string whichrank=to_string(rank+1);
	order.open("output."+shape+".rank"+whichrank+".order.dat");

	if (order.is_open()){} 
	else cerr << "PROBLEM: Unable to open order.out" << endl;
	//this->GetIndividuo(0).PrintIndiv();
	for(int i=0;i<N_city;i++){
		order<<this->GetIndividuo(0).GetElement(i)<<endl;
	}
	order<<this->GetIndividuo(0).GetElement(0)<<endl;
	order.close();
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

double Gene_Alghoritm::GetPos(int i, int j)
{
	if(i>=N_city){cout<<"Error. indice is bigger then the number of cities"<<endl;}
	if(j>=2){cout<<"Error, there are only two dimentions"<<endl;}
	return pos[i][j];
}

double Gene_Alghoritm::GetPos_Start(int i, int j)
{
	if(i>=N_city){cout<<"Error. indice is bigger then the number of cities"<<endl;}
	if(j>=2){cout<<"Error, there are only two dimentions"<<endl;}
	return pos_start[i][j];
}

void Gene_Alghoritm::SetPos_Start(int i, int j, double position)
{
	if(i>=N_city){cout<<"Error. indice is bigger then the number of cities"<<endl;}
	if(j>=2){cout<<"Error, there are only two dimentions"<<endl;}
	pos_start[i][j]=position;	
	pos[i][j]=position;	
}



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
void Gene_Alghoritm::City_initialize(){    //initialization of random distribution of cities
	
	ofstream cities_pos;
	double theta=0,x=0,y=0;
	double x0=0,y0=0;
	
	string whichrank=to_string(rank+1);

	//generate random position in a square or on a circle
	cities_pos.open("output."+shape+".rank"+whichrank+".pos");
	if (cities_pos.is_open()){} 
	else cerr << "PROBLEM: Unable to open cities_pos" << endl;
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


void Gene_Alghoritm::Print_Cities(){
	ofstream cities_pos;

	string whichrank=to_string(rank+1);

	//generate random position in a square or on a circle
	cities_pos.open("output."+shape+".rank"+whichrank+".pos");
	if (cities_pos.is_open()){cout <<"output."+shape+".rank"+whichrank+".pos"<<endl;} 
	else cerr << "PROBLEM: Unable to open random.out" << endl;

	for(unsigned int i=0;i<N_city;i++)
	{
		cities_pos<<this->GetPos_Start(i,0)<<","<<this->GetPos_Start(i,1)<<endl;
	}
		cities_pos<<this->GetPos_Start(0,0)<<","<<this->GetPos_Start(0,1)<<endl;
	cities_pos.close();
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
void Gene_Alghoritm::Check_all(){
for(unsigned int i=0;i<N_indiv;i++){population[i].Check();} //Check if all individuals are correct
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


