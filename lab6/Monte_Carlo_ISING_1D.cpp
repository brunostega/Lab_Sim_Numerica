#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "lib.h"
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char *argv[])
{ 
  Metro_Gibbs =argv[1];
  Input(BoolEquilib); //Inizialization
  Measure();
  cout << "energy at beginning of the simulation = " << walker[iu]/(double)nspin << endl;

  for(int t=0;t<=step_temp;t++){
    temp=temp_iniz+t*delta_T;
    cout<<"temperature="<<temp<<endl;
    beta=1./temp;
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep)
      {
        Move(Metro_Gibbs);
        Measure();
        Accumulate(); //Update block averages
      }
      Averages(iblk);   //Print results for current block
      accepted=0;
      attempted=0;
    }
  }

  return 0;
}



void Input(bool BoolEquilib)
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
  temp=temp_iniz;
  cout << "Initial Temperature = " << temp << endl;
  ReadInput >> temp_final;
  ReadInput >> step_temp;
  delta_T=(double)(temp_final-temp_iniz)/(double)step_temp;
  cout<<"dT="<<delta_T<<endl<<endl;
  beta = 1.0/temp_iniz;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> BoolStart;

  ReadInput >> BoolEquilib;

  ReadInput >> Nequilib_step;

  if(Metro_Gibbs=="metro"){
    cout << "The program perform Metropolis moves" << endl;
    cartella="outputMetropolis";
    cout<<cartella<<endl;
  }
  else if(Metro_Gibbs=="gibbs")
  {
    cout << "The program perform Gibbs moves" << endl;
    cartella="outputGibbs";
    cout<<cartella<<endl;
  }
  else{cout<<"problem! argv[1] was not either metro or gibbs"<<endl;}
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


  //Initialize markov chain
  Boltzmann_prob Boltz_prob;

  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

  //Equilibration
  if(BoolEquilib==1)
  {
    cout<<"Equilibrating the system starting from 2 different configuration. Then saving last final configuration"<<endl;
    Equilibration_LowT(Nequilib_step);  //equilibrate the system starting from low temperature configuration
    Equilibration_HighT(Nequilib_step);  //equilibrate the system starting from high temperature configuration
    ConfFinal();                        //Write final configuration (high temperature one)
  }
  else if(BoolEquilib==0)
  {
    if(BoolStart==0)
    {
      cout<<endl<<"Starting the system from a random configuration. Attention! the system is not equilibrated"<<endl;
      for (int i=0; i<nspin; ++i)
      {
        if(rnd.Rannyu() >= 0.5) s[i] = 1;
        else s[i] = -1;
      }
    }

    else
    {
      cout<<"Starting from a previus (already equilibrated) configuration"<<endl;
      ifstream config;
      config.open(cartella+"/config.final");
      if(config.is_open())
      {
        for(int i=0;i<nspin;i++){
          config >> s[i];
        }
      }      
      else{cout<<"problem with file config.final"<<endl;}
      config.close();
    }
  }
  else{cout<<"Problem! argv[2] was not either 0 or 1"<<endl;}
  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Equilibration_LowT(int Nequilib_step){    //equilibrate the system starting from 2 configuration: infinite T and low T
  ofstream Ene_LowT;

  //Equilibration from Low temperature-->all spin are parallel
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = 1;
  }
  Ene_LowT.open(cartella+"/Equilib_LowT.ene");
  for(int i=0; i<Nequilib_step; i++){
    Move(Metro_Gibbs);
    Measure();
    Ene_LowT<<walker[iu]/(double)nspin<<endl;
    //cout<<walker[iu]/(double)nspin<<endl;
  }
  Ene_LowT.close();
}


void Equilibration_HighT(int Nequilib_step){    //equilibrate the system starting from configuration of infinite temperature
  ofstream  Ene_HighT;

  //Equilibration from High temperature-->all spin are random
  for (int i=0; i<nspin; ++i)
  {
    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }
  Ene_HighT.open(cartella+"/Equilib_HighT.ene");
  for(int i=0; i<Nequilib_step; i++){
    Move(Metro_Gibbs);
    Measure();
    Ene_HighT<<walker[iu]/(double)nspin<<endl;
  }
  Ene_HighT.close();
}



void Move(string Metro_Gibbs)   //if metro==1 uses metropolis alghoritm, otherwise uses Gibbs
{
  int o;
  double p, sm, DeltaE;
  double energy_up, energy_down;
  int appo=0;  //to eventually is all spins maneged to be flipped
  double K, r;
  int r_int;
  for(int i=0; i<nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);
    attempted+=1;
    if(Metro_Gibbs=="metro") //Metropolis
    {
      //genero un numero casuale tra 1 e -1(nuovo spin)
      r_int=rnd.UpDown();
      r=rnd.Rannyu();
      //calcolo la differenza di energia tra vecchia e nuova configurazione (usando le periodic boundary condition)
      DeltaE=Boltzmann(r_int,o)-Boltzmann(s[o],o);
      //metropolis con un solo passo usando come distribuzione di prob quella di boltzmann
      if(1<=exp(-beta*DeltaE)){
        appo++;
		    s[o]=r_int;
        accepted+=1;
	    }
	    else{
        //cout<<exp(-beta*DeltaE)<<", "<<r<<endl;
        if(exp(-beta*DeltaE)>r){
          appo++;
          s[o]=r_int;  
          accepted+=1;
        }
        else{appo++;}
      }
    } 
    else if(Metro_Gibbs=="gibbs")//Gibbs sampling
    {
    K=2*beta*J*(s[Pbc(o-1)] + s[Pbc(o+1)])+2*beta*h;
    r=rnd.Rannyu();

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
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
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
    
  //cout << "Block number " << iblk << endl;
  //cout << "Acceptance rate " << accepted/attempted << endl << endl;
  
  //Energy
  stima_u = blk_av[iu]/blk_norm; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);

  //Heat capacity
  stima_c=(blk_av[ic]/blk_norm-stima_u*stima_u)*beta*beta; 
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);

  //Magnetization
  stima_m = blk_av[im]/blk_norm; 
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);
  
  //Susceptibility
  stima_x = beta*blk_av[ix]/blk_norm; 
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);


  E_final.open(cartella+"/output.E_final.dat", ios::app);
  Heat_final.open(cartella+"/output.Heat_final.dat", ios::app);
  M_final.open(cartella+"/output.M_final.dat", ios::app);
  X_final.open(cartella+"/output.X_final.dat", ios::app);

  //print last values of block averages for observable per spin
  if(iblk==nblk){
    E_final<<temp<<","<<glob_av[iu]/(double)iblk/nspin<<","<<err_u/nspin<<endl;  
    Heat_final<<temp<<","<<glob_av[ic]/(double)iblk/nspin<<","<<err_c/nspin<<endl;
    M_final<<temp<<","<<glob_av[im]/(double)iblk/nspin<<","<<err_m/nspin<<endl;
    X_final<<temp<<","<<glob_av[ix]/(double)iblk/nspin<<","<<err_x/nspin<<endl;
  }

  E_final.close();
  Heat_final.close();
  M_final.close();
  X_final.close();

    //cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open(cartella+"/config.final");
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

