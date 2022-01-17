#ifndef _cercaZeri_h_
#define _cercaZeri_h_

class funzioneBase{
	public:

	virtual double eval(double) const=0;

};

class seno: public funzioneBase{
	
	virtual double eval(double) const;
};

class coseno: public funzioneBase{
	
	virtual double eval(double) const;
};


class parabola: public funzioneBase{

	public:
	
	parabola(double, double, double);
	~parabola();

	double getA();
	double getB();
	double getC();

	virtual double eval(double) const;

	private:
	double m_a, m_b, m_c;

};

class gaussiana: public funzioneBase{

	public:
		gaussiana(double, double);
		~gaussiana();
	
		double getMean();
		double getSigma();
		virtual double eval(double) const;

	private:
		double _mean,_sigma;
};



#endif
