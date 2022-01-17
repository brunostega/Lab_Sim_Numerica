#include "../random.h"
#include "../block_stat.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <complex>
#include <bits/stdc++.h>

double pi_greco=3.14159;
double beta0=0.001;

class Individuo{
    
    public:
    Individuo();
    Individuo(unsigned int);
    ~Individuo();
    int GetElement(int);            //get element of position
    void SetElement(int,int);       //set position element
    void SetPosition(std::vector< std::vector<double> >);
    void PrintOrder();              //print cities order
    void Fitness();                 //estimate fitness
    double GetLoss();               //get fitness
    void SetLoss(double);           //set fitness
    void Check();                   //checks if order of cities is correct
    void Mutation1(Random *);
    void Mutation2(Random *);
    void Mutation3(Random *);

    private:
    unsigned int N_city;
    double L, exp;
    std::vector<int> Order;
    std::vector< std::vector<double> > pos;
};

class Travel_Salesman{

    public:
    Travel_Salesman(std::string, std::string);     //inizialize positions of cities, beta 
    ~Travel_Salesman();
    void Sim_Annealing();  //makes montecarlo N_steps at fixed termperature
    void City_Initialize();
    void Beta_Initialize();
    Individuo Mutation();

    private:
    double Boltz_weight, Beta_step;
    Individuo Order;
    unsigned int N_city, N_beta, N_step;
    std::vector<double> beta;  
    std::string shape, beta_law;
    Random rnd;
};

