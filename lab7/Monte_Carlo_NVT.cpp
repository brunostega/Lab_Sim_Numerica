
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    appo=0;
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
      appo++;

      if(istep%10 == 0){
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  OldFinal();

  return 0;
}


void Input(void)
{
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> Bool;

  ReadInput >> Eq_steps;

  ReadInput >> TermBool;

  cout<< "if TermBool variable is 1, the program makes the equilibration. Term bool= "<<TermBool<<endl;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();

//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
 
  n_props = 2; //Number of observables

//measurement of g(r)
  igofr = 2;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

  //equilibration
  if(TermBool==1)
  {
    Equilibration();
  }
//Read initial configuration

  if(Bool==false){
    cout << "Read initial configuration from file config.0 " << endl << endl;
    ReadConf.open("config.0");
    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
    }
    ReadConf.close();
    //Prepare initial velocities and old position xold
  }

  if(Bool==true){
    ifstream ReadOld;
    cout << "Read initial configuration from file config.final " << endl << endl;
    ReadConf.open("config.final");
    ReadOld.open("old.final");

    for (int i=0; i<npart; ++i){
      ReadConf >> x[i] >> y[i] >> z[i];
      x[i] = x[i] * box;
      y[i] = y[i] * box;
      z[i] = z[i] * box;
      ReadOld >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
    }
    ReadConf.close();
    //Prepare initial velocities and old position xold
  }
}


void Equilibration(void){   //equilibrate the system 

  ifstream ReadConf;
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

  cout<<"Equilibration of the system"<<endl;
  for(int i=0;i<Eq_steps;i++){
      if(i%(Eq_steps/10)==0){cout<<"step = "<<i<<endl;}
      Move();
      Measure_Equilib();
  }
  ConfFinal();
  OldFinal();
}


void Move(void)
{
  int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;

  for(int i=0; i<npart; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (int)(rnd.Rannyu()*npart);

    //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

    //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

    //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu())  
    {
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;
    
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
}


double Boltzmann(double xx, double yy, double zz, int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (int i=0; i<npart; ++i)
  {
    if(i != ip)
    {
      //distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}


void Measure_Equilib()  //used to measure and save output in a different file (equilibration output)
{
  ofstream en,vir;
  en.open("output.equilib_energy.dat",ios::app);
  vir.open("output.equilib_press.dat",ios::app);
  int bin, position;
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);  

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      //update of the histogram of g(r)
      position=int(dr*nbins/(box/2));
      if(igofr+position>m_props){cout<<igofr+position<<endl;}
      walker[igofr+position]+=2;

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

        //contribution to energy and virial
        v += vij;
        w += wij;
      }
    }          
  }
  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;

  en <<walker[iv]/(double)npart+vtail<<endl;
  vir << rho * temp + (walker[iw] + ptail * (double)npart) / vol<<endl;
  en.close();
  vir.close();
}



void Measure()   //use to measure variables and save output in files
{
  ofstream en,vir,g;
  en.open("output.energia_istant.dat",ios::app);
  vir.open("output.virial_istant.dat",ios::app);
  g.open("output.gofr.0",ios::app);

  int bin, position;
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i)
  {
    for (int j=i+1; j<npart; ++j)
    {
      // distance i-j in pbc
      dx = Pbc(x[i] - x[j]);
      dy = Pbc(y[i] - y[j]);
      dz = Pbc(z[i] - z[j]);  

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

       position=int(dr*nbins/(box/2));
      if(igofr+position>m_props){cout<<igofr+position<<endl;}
      walker[igofr+position]+=2;

      if(dr < rcut)
      {
        vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
        wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

        //contribution to energy and virial
        v += vij;
        w += wij;
      }
    }          
  }
  walker[iv] = 4.0 * v;
  walker[iw] = 48.0 * w / 3.0;

  en <<walker[iv]/(double)npart+vtail<<endl;
  vir << rho * temp + (walker[iw] + ptail * (double)npart) / vol<<endl;
  en.close();
  vir.close();
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
    //cout<<i<<endl;
    //cout<<blk_av[i]<<endl;
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{    
  double r, gdir;
  ofstream Gofr, Gave, Epot, Pres;
  const int wd=12;
    
  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted/attempted << endl << endl;
  
  Epot.open("output.epot.0",ios::app);
  Pres.open("output.pres.0",ios::app);
  Gofr.open("output.gofr.0",ios::app);
  Gave.open("output.gave.0",ios::app);
  
  stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
  glob_av[iv] += stima_pot;
  glob_av2[iv] += stima_pot*stima_pot;
  err_pot=Error(glob_av[iv],glob_av2[iv],iblk);
  
  stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
  glob_av[iw] += stima_pres;
  glob_av2[iw] += stima_pres*stima_pres;
  err_press=Error(glob_av[iw],glob_av2[iw],iblk);

  for(int i=igofr; i<igofr+nbins;i++){
    double r=box/2*(i-igofr)/nbins;
    double dr=box/2/nbins;
    double stimag=blk_av[i]/blk_norm/(rho*npart*pi*4/3*(pow(r+dr,3)-pow(r,3)));
    glob_av[i]+=stimag;
    glob_av2[i]+=stimag*stimag;
    err_gdir=Error(glob_av[i],glob_av2[i],iblk);
  }
//Potential energy per particle
  cout <<  iblk <<  "," << glob_av[iv]/(double)iblk << "," << err_pot << endl;
//Pressure
  cout << iblk <<  ","  << glob_av[iw]/(double)iblk << "," << err_press << endl;

//Potential energy per particle
  Epot <<  iblk <<  ","  << glob_av[iv]/(double)iblk << "," << err_pot << endl;
//Pressure
  Pres << iblk <<  ","  << glob_av[iw]/(double)iblk << "," << err_press << endl;
//g(r)
  if(iblk==nblk){
    cout<<iblk<<". Sono qui dentro"<<endl;
    for(int i=igofr;i<igofr+nbins;i++){
      double r=box/2*(i-igofr)/nbins;
      Gave << r << "," << glob_av[i]/(double)iblk << "," <<Error(glob_av[i],glob_av2[i],iblk) <<endl;
    }
  }

  cout << "----------------------------" << endl << endl;

  Epot.close();
  Pres.close();
  Gofr.close();
  Gave.close();
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<npart; ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void OldFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file old.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}


void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}


