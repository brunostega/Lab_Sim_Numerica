
#ifndef __NVT__
#define __NVT__

//Random numbers
#include "../random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000;
int n_props, iv, iw, igofr;
double vtail,ptail,bin_size,nbins,sd;
double walker[m_props];
bool Bool,TermBool,VelvetBool;
int appo=0;

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_pot,stima_pres,err_pot,err_press,err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];
double xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];


// thermodynamical state
int npart;
double beta,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, Eq_steps;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Equilibration();
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void MoveVelvet(void);
void ScaleVelocity();
void RandomVelocity();
double Force(int, int);
void ConfFinal(void);
void OldFinal(void);
void ConfXYZ(int);
void Measure(void);
void Measure_Equilib(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);



#endif

