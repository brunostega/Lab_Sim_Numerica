#include "function.h"

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();
  // destructor
  ~Random();
  // methods
  void Initialize();
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Gauss(double mean, double sigma);
  double Exponential(double lambda);
  double Cauchy(double mean, double gamma);
  double Dado(int N);
  double acceptReject(funzioneBase *, double, double, double);
  void RW_discrete(double, double *,int,  int);
  void RW_continum3D(double, double *,int,  int);
  int UpDown();


};

#endif // __Random__
