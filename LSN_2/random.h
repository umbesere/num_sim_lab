/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructor
  Random();
  // destructor
  ~Random();
  
  /** methods **/
  void SetRandom(int * , int, int);
  void SaveSeed();
	// generating uniform double or int rnd numbers
  double Rannyu(void);
  double Rannyu(double min, double max);
  int RannyuDiscr(void);
  int RannyuDiscr(int min, int max);
	// generating rnd numbers that follow a precise pdf
  double Gauss(double mean, double sigma);
  double Exponential(double gamma); 
  double Cauchy(double mean, double gamma);
	double P_x(void);
  // generating rnd angles and rnd solid angles
  double Angle(void);
  void SolidAngle(double &theta, double &phi);
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
