#include "random.h"
#include <string>
#include <fstream>

#ifndef __BlockStat__
#define __BlockStat__

class BlockStat {

    private:
    
    public:
        int _N_blocks=0, _N_throws=0, _n=0;    //# of blocks, # of throws, # of throws in each block
        int * _x;             //progressive number of throws per block      
        double * _mean_prog;  //progressive sum of mean
        double * _mean2_prog; //progressive sum of squared mean
        double * _err_prog;   //progressive standard deviation

        double _final_mean,_final_error;

        double error(double, double, int);
        double get_final_mean();
        double get_final_error();
        void progressive_mean(double *);
        void progressive_mean_squared(double *);
        void progressive_error();
        void print_output(const std::string);
        void block_stat_all(const std::string, double*);
        void reset();

        //constructor
        BlockStat(int, int);
        //Destructir
        ~BlockStat();
};

#endif // __Random__



double pi(int,double, double, Random &);
