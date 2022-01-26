#include "../../random.h"
#include "../../block_stat.h"
#include <string>
#include "mpi.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <complex>
#include <bits/stdc++.h>

double pi_greco=3.14159;


class Individuo{

    public:
    Individuo();
    Individuo(unsigned int);
    ~Individuo();
    void Shuffle_Indiv();
    void Check();
    double GetLoss();
    void SetLoss(double);
    int GetElement(int);
    void SetElement(int,int);
    void PrintIndiv();


    private:
    unsigned int N_city;
    double L;
    int exp;
    Random rnd;
    std::vector<int> indiv;
};


class Gene_Alghoritm{

    public:
    Gene_Alghoritm(int);
    ~Gene_Alghoritm();
    void Out_L(int);
    void Print_Order();
    void GridSearch(std::string);
    void New_Generation();
    void City_initialize();
    void Reset_Cities();
    void Print_Cities();
    void PrintPopulation();
    void Sort();
    void fitness();
    int select();
    void reset_L();
    void crossover1();
    void crossover2();
    void mutation();
    void Check(Individuo);
    void Check_all();
    void mutation1(Individuo &);
    void mutation2(Individuo &);
    void mutation3(Individuo &);
    void mutation4(Individuo &);
    double GetPos(int, int);
    double GetPos_Start(int, int);
    void SetPos_Start(int, int, double);
    int Get_NGeneration();
    int Get_NCity();
    int Get_NIndiv();
    int GetExp();
    int GetCross();
    Individuo GetIndividuo(int);
    void SetIndividuo(Individuo, int);
    void SetPm(double);
    void SetPc(double);

    Random rnd;
    std::string shape;

    private:

    //cities
    unsigned int N_city;


    //genetic alghoritm data
    unsigned int N_indiv, N_generation;
    std::vector<Individuo> population;
    std::vector<Individuo> population_start;
    std::vector<std::vector<double> > pos;
    std::vector<std::vector<double> > pos_start;
    std::vector<double> L;
    double P_c, P_m;

    int exp, fit_exp, cross;

    //mpi
    int rank;
};



