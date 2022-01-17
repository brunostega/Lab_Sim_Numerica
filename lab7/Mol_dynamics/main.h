#include "../../random.h"
#include "../../block_stat.h"
#include <string>


//functions
void ChangeConfig(void);

const int m_props=1000;
const int m_part=108;

// simulation
int nstep, iprint, seed;
double delta;

//block statistics
int N_blocks=100;        //#blocks
const double pi_greco=3.1415927;



//class
class MolDyn{
    public:
        MolDyn(std::string);
        ~MolDyn();
        void Input(void);
        void Move(void);
        void ConfFinal(void);
        void OldFinal(void);
        void ConfXYZ(int);
        void RandomVelocity(void);
        void ScaleVelocity(void);
        void Reset(int);
        void Accumulate(void);
        void Averages(int);
        void Measure(void);
        void Measure_equilib(void);
        void Equilibration();
        double getEpot();
        double getEkin();
        double getEtot();
        double getTemp();
        double getPressure();
        double Force(int, int);
        double Pbc(double);
        double Error(double,double,int);


    private:
        //parameters, observables
        int n_props, igofr;
        int iv,ik,it,ie;
        double stima_pot, stima_kin, stima_etot, stima_temp, stima_P;
        double bin_size,nbins;
        bool TFbool;
        Random rand;
        std::string percorso;
        double walker[m_props];


        //equilibration
        bool termBool;
        int N_cicle, N_equil_step;

        // averages
        double acc,att;
        double blk_av[m_props],blk_norm,accepted,attempted;
        double glob_av[m_props],glob_av2[m_props];

        //configuration
        double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
        double vx[m_part],vy[m_part],vz[m_part];

        // thermodynamical state
        int npart;
        double energy,temp,vol,rho,box,rcut;

};

