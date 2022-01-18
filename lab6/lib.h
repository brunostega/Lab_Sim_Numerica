#include "../random.h"
#include <string>




//CLASS PROBABILITY (mother class to define functions)

class probability{

	public:
	virtual double eval(double)=0;
};


//HYDROGEN BOLTZMANN FUNCTIONS
class Boltzmann_prob: public probability{
	virtual double eval(double);
};

//CLASS ACCEPTANCE (mother class to define functions): TRANSITION PROBABILITY FUNCTION

class acceptance{
	public:
    acceptance();
    ~acceptance();
	virtual double eval(double)=0;

    protected:
    double _delta;         //lenght of the step for the random walk
    Random _rnd;

};

class newspin: public acceptance{
	
};



//CLASS METROPOLIS for metropolis algoritm

class Metropolis{

	public:
	Metropolis(int, probability *, acceptance *);  //costruttore
	~Metropolis();                                 //distruttore 
	double Acceptance(double, double);         //Acceptance function
	double MR_alghoritm(double);             //Metropolis alghoritm with "_passi" steps with the transition probability distrib given in constructor

	private:

	int _passi;                    //number of steps of Metropolis algoritm
	Random _rnd;                   //Random class variable to generate pseudo-casual numbers
	probability * _probability;    //probability density function to sample
    acceptance * _acceptance;      //acceptance function (transition probability distribution)
};




