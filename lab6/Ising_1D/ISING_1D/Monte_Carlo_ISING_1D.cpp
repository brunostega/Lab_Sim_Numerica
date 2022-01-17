#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "lib.h"
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  Equilibration();
  cout<<delta_T<<endl;
  for(int t=0;t<=step_temp;t++){
    temp=temp_iniz+t*delta_T;
    cout<<temp<<endl;
    beta=1./temp;
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move(metro);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
    }
    ConfFinal(); //Write final configuration
  }

  return 0;
  
}


void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  rnd.Initialize();

//Read input informations

  ReadInput.open("input.dat");

  ReadInput >> temp_iniz;
  cout << "Temperature = " << temp << endl;
  ReadInput >> temp_final;
  ReadInput >> step_temp;
  delta_T=(double)(temp_final-temp_iniz)/(double)step_temp;

  beta = 1.0/temp_iniz;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> BoolStart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Initialize markov chain
Boltzmann_prob Boltz_prob;
int N=1; //making only one filp per atom with markov chain
//Metropolis Metro(N, Boltz_prob, T)

//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

//initial configurationifference in energy between two states
  if(BoolStart==0){
    cout<<"sono qui"<<endl;
    for (int i=0; i<nspin; ++i)
    {
      if(rnd.Rannyu() >= 0.5) s[i] = 1;
      else s[i] = 1;
    }
  }
  else{
    
    ifstream config;
    config.open("config.final");
    for(int i=0;i<nspin;i++){
      config >> s[i];
      }
  }
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}

void Equilibration(){
  ofstream Ene;
  Ene.open("Equilib.ene");
  for(int i=0; i<100000; i++){
    Move(metro);
    Measure();
    Ene<<walker[iu]/(double)nspin<<endl;
  }
}


void Move(int metro)
{
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;
  int appo=0;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) //Metropolis
    {

      //genero un numero casuale tra 1 e -1(nuovo spin)
      int r_int=rnd.UpDown();
      double r=rnd.Rannyu();
      //calcolo la differenza di energia tra vecchia e nuova configurazione (usando le periodic boundary condition)
      double DeltaE=Boltzmann(r_int,o)-Boltzmann(s[o],o);
      //metropolis con un solo passo usando come distribuzione di prob quella di boltzmann
      if(1<=exp(-beta*DeltaE)){
        //cout<<exp(-beta*DeltaE)<<endl;
        appo++;
        //cout<<endl<<s[o]<<"    "<<r_int<<endl;
		    s[o]=r_int;
	    }
	    else{
        if(exp(-beta*DeltaE)<r){
          //cout<<"non è cambiato lo spin"<<endl;
        }
        else{
        appo++;
		    s[o]=r_int;          
        }
      }
    } 
    else //Gibbs sampling
    {
    double K=2*beta*J*(s[Pbc(o-1)] + s[Pbc(o+1)])+2*h;
    double r=rnd.Rannyu();
  
    if(r<=1.0/(1+exp(-K))){s[o]=1;}
    else{s[o]=-1;}

    }
  }
}


double Boltzmann(int sm, int ip)      //energy of one spin (sm)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}


void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m+=s[i];
// INCLUDE YOUR CODE HERE
  }
  walker[iu] = u;
  walker[ic]=u*u;
  walker[im]=m;
  walker[ix]=m*m;
// INCLUDE YOUR CODE HERE
}




void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}




void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}




void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   ofstream E_final, Heat_final, M_final,X_final;
   const int wd=12;
    
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u=Error(glob_av[iu],glob_av2[iu],iblk);
    Ene << iblk << ","<< stima_u <<"," << glob_av[iu]/(double)iblk <<","<< err_u << endl;
    Ene.close();

    
    stima_c=(blk_av[ic]/blk_norm/(double)nspin/(double)nspin-stima_u*stima_u)*beta*beta;
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c=Error(glob_av[ic],glob_av2[ic],iblk);

    stima_m = blk_av[im]/blk_norm; //Energy
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m=Error(glob_av[im],glob_av2[im],iblk);
    
    stima_x = beta*blk_av[ix]/blk_norm; //Energy
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x=Error(glob_av[ix],glob_av2[ix],iblk);


    E_final.open("E_final.dat", ios::app);
    if(iblk==nblk){
      E_final<<temp<<","<<glob_av[iu]/(double)iblk<<","<<err_u<<endl;
    }
    E_final.close();

    Heat_final.open("Heat_final.dat", ios::app);
    if(iblk==nblk){
      Heat_final<<temp<<","<<glob_av[ic]/(double)iblk<<","<<err_c<<endl;
    }
    Heat_final.close();

    M_final.open("M_final.dat", ios::app);
    if(iblk==nblk){
      M_final<<temp<<","<<glob_av[im]/(double)iblk<<","<<err_m<<endl;
    }
    M_final.close();

    X_final.open("X_final.dat", ios::app);
    if(iblk==nblk){
      X_final<<temp<<","<<glob_av[ix]/(double)iblk<<","<<err_x<<endl;
    }
    X_final.close();


// INCLUDE YOUR CODE HERE

    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}



int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

