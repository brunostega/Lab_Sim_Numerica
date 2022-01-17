#include "/home/bruno/Scrivania/universita/magistrale1/lab_simulaz/esercizi/random.h"
#include <string>


//CLASSE POSIZIONE (useful to pass from cartesian coordinates to polar coordinates)

class posizione{

	public:

	posizione();	//costruttori
	posizione(double x); 	//costruttori
	~posizione(); 	//distruttore


	double getX() const;

	//void setR(double )const;

	double getDistanza(const posizione&) const;

	protected:
	
	double m_x, m_y, m_z;

};

//CLASS PROBABILITY (mother class to define functions)

class probability{

	public:
	probability();
	~probability();
	virtual double eval(posizione)=0;
};


//HYDROGEN PROBABILITY FUNCTIONS
class Wave_func: public probability{
	
	public:
	Wave_func(double, double);
	~Wave_func();
	void reset(double, double);
	double getMu();
	double getSigma();
	virtual double eval(posizione);

	private:
	double _sigma,_mu;
};





//CLASS ACCEPTANCE (mother class to define functions): TRANSITION PROBABILITY FUNCTION

class acceptance{
	public:
    acceptance();
    ~acceptance();
	virtual posizione eval(posizione)=0;

    protected:
    double _delta;         //lenght of the step for the random walk
    Random _rnd;

};

class T_uniform: public acceptance{      //Transition probability function
    public:
    T_uniform(double);     
    T_uniform();
    ~T_uniform();
	virtual posizione eval(posizione);
};

class T_gaussian: public acceptance{     //Transition probability function
    public:  
    T_gaussian(double);
    T_gaussian();
    ~T_gaussian();
	virtual posizione eval(posizione);
};



//CLASS METROPOLIS for metropolis algoritm

class Metropolis{

	public:
	Metropolis(int, probability *, acceptance *);  //costruttore
	~Metropolis();                                 //distruttore 
	void Find_delta();                             //finding scaling factor with 50% rule
	posizione Acceptance(posizione, posizione);    //Acceptance function
	posizione MR_alghoritm();             //Metropolis alghoritm with "_passi" steps with the transition probability distrib given in constructor
	void setX_start(posizione);
	void SetProbabillity(probability *);

	private:

	int _passi;                    //number of steps of Metropolis algoritm
	Random _rnd;                   //Random class variable to generate pseudo-casual numbers
	probability * _probability;    //probability density function to sample
    acceptance * _acceptance;      //acceptance function (transition probability distribution)
	posizione _x_start;
};




