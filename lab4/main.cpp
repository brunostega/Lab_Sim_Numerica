//#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>       // cin, cout: Standard Input/Output Streams Library
#include <fstream>        // Stream class to both read and write from/to files.
#include <cmath>          // rint, pow
#include "main.h"
#include "lib.h"
#include "../random.h"
#include "../block_stat.h"



using namespace std;


int main(int argc, char *argv[]){ 
  string percorso;
  cout<<argc<<endl;
  if(argc==2){
    percorso=argv[1];
    cout<<percorso<<endl;
  }
  else{
    percorso="../esercizi/lab4";
    cout<<percorso<<endl;
  }


  //INITIALIZING SOME VARIABLES

  const string  file_epot=percorso+"/ave_epot.out";      //output file for potential energy
  const string  file_ekin=percorso+"/ave_ekin.out";      //output file for kinetic energy
  const string  file_etot=percorso+"/ave_etot.out";      //output file for total energy
  const string  file_temp=percorso+"/ave_temp.out";      //output file for temperature
  const string  file_pressure=percorso+"/ave_pressure.out";      //output file for temperature


	double ave_epot[N_blocks]={};                //block mean for the different properties
	double ave_ekin[N_blocks]={};                //
	double ave_etot[N_blocks]={};                //
	double ave_temp[N_blocks]={};                //
	double ave_pressure[N_blocks]={};                //


  MolDyn Mdin=MolDyn(percorso);      //Inizialization with equilibration of molecular dynamics class
  int L=nstep/N_blocks;              //Number of steps per block

  BlockStat blockstat(N_blocks,nstep);         //Initializing class for block statistics


  //EXPERIMENT

  cout<<endl<<"Actual experiment"<<endl;
  for(int istep=0; istep < N_blocks; ++istep){
    if((istep+1)%iprint == 0){cout << "Number of time-steps: " << istep+1 << endl;}
    for(int j=0; j<L;j++){
      Mdin.Move();            //Move particles with Verlet algorithm
      Mdin.Measure();         //Properties measurement 
      ave_epot[istep]+=Mdin.getEpot();   //upgrade potential energy 
      ave_ekin[istep]+=Mdin.getEkin();   //upgrade kinetic energy
      ave_etot[istep]+=Mdin.getEtot();   //upgrade total energy
      ave_temp[istep]+=Mdin.getTemp();   //upgrade temperature 
      ave_pressure[istep]+=Mdin.getPressure();   //upgrade temperature 

    }

    //mean estimate per block
    ave_epot[istep]/=L;          
    ave_ekin[istep]/=L;
    ave_etot[istep]/=L;
    ave_temp[istep]/=L;
    ave_pressure[istep]/=L;
  }

  //block statistics and printing average and error(per block) in files
  blockstat.block_stat_all(file_epot, ave_epot);
  blockstat.block_stat_all(file_ekin, ave_ekin);
  blockstat.block_stat_all(file_etot, ave_etot);
  blockstat.block_stat_all(file_temp, ave_temp);
  cout<<"media finale a blocchi = "<<blockstat.get_final_mean()<<" +- "<<blockstat.get_final_error()<<endl<<endl<<endl;
  blockstat.block_stat_all(file_pressure, ave_pressure);

  Mdin.OldFinal();
  Mdin.ConfFinal();         //Write final configuration to restart

  return 0;
}


MolDyn :: MolDyn(const string path){         //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  Random rand_appo;
  rand=rand_appo;
  rand.Initialize();
  //double ep, ek, pr, et, vir;
  percorso=path;


  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //ReadInput.open(percorso+"/input.dat"); //Read input
  ReadInput.open("input.dat"); //Read input

    ReadInput >> temp;
    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;
    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;
    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> iprint;
    ReadInput >> TFbool;
    ReadInput >> termBool;
    ReadInput >> N_cicle;
    ReadInput >> N_equil_step;
    cout <<"Boolean variable (true=1, false=0)="<<TFbool<<". Chooses to start from config0 or from configuration from last run"<<endl;
    cout <<"Boolean variable (true=1, false=0)="<<termBool<<". Chooses to equilibrate of the system"<<endl;
    cout <<"If requested the program equilibrate the system (before making the experiment) with "<<N_equil_step<<" steps for every cicle ("<<N_cicle<<" cicles)"<<endl;
    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;

  ReadInput.close();

  //Prepare array for measurements

  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables


  //Equilibrating the system

  if(termBool==true){
    this->Equilibration();
  }


  //Read initial configuration

  if(TFbool==false){
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
    this->RandomVelocity();
  }

  if(TFbool==true){
    ifstream ReadOld;
    cout << "Read initial configuration from file config.final " << endl << endl;
    ReadConf.open(percorso+"/config.final");
    ReadOld.open(percorso+"/old.final");

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
    //this->ScaleVelocity();
  }
}

MolDyn::~MolDyn(){}

void MolDyn::ScaleVelocity(void){       //rinormalize velocity to equilibrate the system to fixed temperature
  double sumv2 = 0.0, fs;

  for(int i=0; i<npart; i++){
    vx[i] = this->Pbc(x[i] - xold[i])/delta;         //le velocità sono tali da non fare passi più lunghi di box/2
    vy[i] = this->Pbc(y[i] - yold[i])/delta;
    vz[i] = this->Pbc(z[i] - zold[i])/delta;
    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;

  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor (making sure that v_tot=sqrt(3T) )
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = this->Pbc(x[i] - vx[i] * delta);
    yold[i] = this->Pbc(y[i] - vy[i] * delta);
    zold[i] = this->Pbc(z[i] - vz[i] * delta);
  }
  return;
}

void MolDyn::RandomVelocity(void){    //uses random velocities to generate the old configuration

  //Prepare initial velocities
  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){
    vx[i] = rand.Rannyu() - 0.5;      
    vy[i] = rand.Rannyu() - 0.5;
    vz[i] = rand.Rannyu() - 0.5; 

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
    double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){  //choosing v such that <v>=0 (not so many particle to assure mean 0)
      vx[i] = vx[i] - sumv[0];
      vy[i] = vy[i] - sumv[1];
      vz[i] = vz[i] - sumv[2];

      sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
    }
    sumv2 /= (double)npart;

    fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor (making sure that v_tot=sqrt(3T) )
    for (int i=0; i<npart; ++i){
      vx[i] *= fs;
      vy[i] *= fs;
      vz[i] *= fs;

      xold[i] = this->Pbc(x[i] - vx[i] * delta);
      yold[i] = this->Pbc(y[i] - vy[i] * delta);
      zold[i] = this->Pbc(z[i] - vz[i] * delta);
    }
  return;
}

void MolDyn::Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = this->Force(i,0);
    fy[i] = this->Force(i,1);
    fz[i] = this->Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = this->Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = this->Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = this->Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = this->Pbc(xnew - xold[i])/(2.0 * delta);         //le velocità sono tali da non fare passi più lunghi di box/2
    vy[i] = this->Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = this->Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double MolDyn::Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = this->Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = this->Pbc( y[ip] - y[i] );
      dvec[2] = this->Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
} 

double MolDyn::getEpot(){
  return stima_pot;
}

double MolDyn::getEkin(){
  return stima_kin;
}

double MolDyn::getEtot(){
  return stima_etot;
}

double MolDyn::getTemp(){
  return stima_temp;
}

double MolDyn::getPressure(){
  return stima_P;
}

void MolDyn::Measure(void){ //Properties measurement

  double v, t, vij, P, Pij;
  double dx, dy, dz, dr;

  ofstream Epot, Ekin, Etot, Temp;

  v = 0.0; //reset observables
  P=0.0;
  t = 0.0;

  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

      dx = this->Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
      dy = this->Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
      dz = this->Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut){
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        Pij=1.0/pow(dr,12) - 0.5/pow(dr,6);
  //Potential energy
       v += vij;
  //Pressure
       P += 48*Pij*1/vol;

      }
    }          
  }

  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart;                  //Potential energy per particle
    stima_kin = t/(double)npart;                  //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart;   //Temperature
    stima_etot = (t+v)/(double)npart;             //Total energy per particle
    stima_P=rho*temp+P;

    return;
}

void MolDyn::Measure_equilib(void){ //Properties measurement
  //int bin;
  double v, t, vij, P, Pij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Temp.open(percorso+"/output_temp.dat",ios::app);

  //Epot.open(percorso+"/output_epot.dat",ios::app);
  //Ekin.open(percorso+"/output_ekin.dat",ios::app);
  //Etot.open(percorso+"/output_etot.dat",ios::app);

  v = 0.0; //reset observables
  P=0.0;
  t = 0.0;
/*
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

      dx = this->Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
      dy = this->Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
      dz = this->Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut){
        vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
        Pij=1.0/pow(dr,12) - 0.5/pow(dr,6)
  //Potential energy
       v += vij;
       P += 48*Pij;
      }
    }          
  }
*/
  //Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  stima_pot = v/(double)npart;                  //Potential energy per particle
  stima_kin = t/(double)npart;                  //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart;   //Temperature
  stima_etot = (t+v)/(double)npart;             //Total energy per particle
  stima_P=rho*temp+P;


  Temp << stima_temp << endl;

  //Epot << stima_pot  << endl;
  //Ekin << stima_kin  << endl;
  //Etot << stima_etot << endl;

  //Etot.close();
  //Epot.close();
  //Ekin.close();
  Temp.close();

  return;
}

void MolDyn::ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  //cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open(percorso+"/config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void MolDyn::OldFinal(void){ //Write final configuration
  ofstream WriteConf;

  //cout << "Print final configuration to file old.final " << endl << endl;
  WriteConf.open(percorso+"/old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void MolDyn::ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double MolDyn::Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);        
}

void MolDyn::Equilibration(){

  ifstream ReadConf;

  cout<<"Thermodynamical equilibration of the system"<<endl;
  int count=0;

  for(int j=1;j<=N_cicle;j++){
    
    if(j==1){       //Read initial configuration
      cout <<"Cicle="<<j <<". Read initial configuration from file config.0 " << endl;
      ReadConf.open("config.0");
      for (int i=0; i<npart; ++i){
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
      }
      ReadConf.close();

      //Prepare initial velocities and old position xold
      this->RandomVelocity();
    }

    else{
      cout<<"Cicle="<<j<<". Read initial configuration from file last configuration (x,x_old) " << endl;
      ifstream ReadOld;
      //cout << "Read initial configuration from file config.final " << endl << endl;
      ReadConf.open(percorso+"/config.final");
      ReadOld.open(percorso+"/old.final");

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
      this->ScaleVelocity();
    }

    for(int istep=1; istep <= N_equil_step; ++istep){
      this->Move();                //Move particles with Verlet algorithm
      this->Measure_equilib();     //Properties measurement
      count++;
    }
    this->OldFinal();
    this->ConfFinal();         //Write final configuration to restart
  }

  cout<<endl<<"Equilibration ended. Num of measures: "<<count<<endl<<endl;
}

