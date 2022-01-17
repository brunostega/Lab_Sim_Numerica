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
double ave_L=0;

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
    Gene_Alghoritm(std::string);
    ~Gene_Alghoritm();
    void GridSearch(std::string);
    void New_Generation();
    void City_initialize(std::string);
    void Reset_Cities();
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
    int GetPos(int, int);
    int Get_NGeneration();
    int Get_NCity();
    int Get_NIndiv();
    int GetExp();
    int GetCross();
    Individuo GetIndividuo(int);
    void SetIndividuo(Individuo, int);
    void SetPm(double);
    void SetPc(double);



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
    Random rnd;
    int exp, fit_exp, cross, grid;
};



