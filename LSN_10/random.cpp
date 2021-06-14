#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "random.h"

using namespace std;

// revisited constructor
Random :: Random( string seedin, string name ){
	m_name = name; 
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input(seedin);
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	for(int i=0; i<20; i++){
		Rannyu();
	}	
}

// destructor with saving seed 
Random :: ~Random(){
	SaveSeed();
}

// SaveSeed
void Random :: SaveSeed(){ 
   ofstream WriteSeed;
   WriteSeed.open(m_name);
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

// SetRandom
void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0];
  l2 = s[1];
  l3 = s[2];
  l4 = s[3];
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

// uniform generator between two double
double Random :: Rannyu(void){

  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));
  return r;
}

// uniform generator in (min, max)
double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

// discrete generator in (0,1)
int Random :: RannyuDiscr(void){
	return int(Rannyu(0.,2.));
};

// discrete generator in (min, max) (int version)
int Random :: RannyuDiscr(int min, int max){
	return int(Rannyu(min, max));
};

// gaussian generators
double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

// stadard normal generator
double Random :: Gauss() {
   return Gauss(0.,1.);
}

// exponential generator
double Random :: Exponential(double gamma){
	return -1/gamma*log(1-Rannyu());
}

// Cauchy-Lorentz generator
double Random :: Cauchy(double mean, double gamma){
	return mean + gamma*tan( M_PI*( Rannyu() - 0.5 ) );
}

// angle generator (without using pi)
double Random :: Angle(void){
	double x,y;
	for(;;){
		x = Rannyu(-1.,1.);
		y = Rannyu(0.,1.);
		if(x*x+y*y <= 1)
			break;
	}
	return 2*acos(x/sqrt(x*x+y*y));
}

// solid angle generator
void Random :: SolidAngle(double &theta, double &phi){
	phi = Rannyu(0., 2*M_PI);
	theta = acos(1 - 2*Rannyu());
};

// p(x) = 2(1-x) generator
double Random :: P_x(void){
	return 1 - sqrt(1-Rannyu());
};