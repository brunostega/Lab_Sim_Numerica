#include "../random.h"
#include "../function.h"
#include <string>

Random initialize_random();

class BlackScholes{
	public:
		BlackScholes(double, double, double, double, double, double, Random &);
		~BlackScholes();
		double Eur_call_direct(double);
		double Eur_put_direct(double);
		double Eur_call_discrete(double, int);
		double Eur_put_discrete(double, int);
		double W(double);
		double S(double);
		double S_discrete(double, int);
		double d1(double, double);
		double d2(double, double);
		double N(double);

	private:
		double _T,_K,_r,_mu,_sigma,_S0;
		Random _rnd;
};

