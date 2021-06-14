#ifndef __Random__
#define __Random__

#include<string>
#include<cstdlib>
#include<vector>

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
  std::string m_name; // seed out file name
  
protected:

public:
  
  // revisited constructor
  Random( std::string name );
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
	// generating numbers following a pdf
  double Gauss();
  double Gauss(double mean, double sigma);
  double Exponential(double gamma); 
  double Cauchy(double mean, double gamma);
	double P_x(void);
  // generating rnd angles and rnd solid angles
  double Angle(void);
  void SolidAngle(double &theta, double &phi);
};

#endif // __Random__

