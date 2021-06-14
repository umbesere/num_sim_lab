#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "random.h"

using namespace std;

// error function to calculate standard deviation from mean and mean squared
double error(double ave, double av2, int n){
	if (n == 0)
		return 0;
	else
		return sqrt((av2 - pow(ave,2))/n);
}
 
int main (int argc, char *argv[]){
	
	// initializing random
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
	Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
	while ( !input.eof() ){
	input >> property;
	if( property == "RANDOMSEED" ){
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	}
	}
	input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	for(int i=0; i<20; i++){
	cout << rnd.Rannyu() << endl;
	 	}

	rnd.SaveSeed();
	cout << endl;
   	
	//** INTEGRAL COMPUTATION **//   	

	// initializing variables
	int M = 100000;							// #points
	int N = 100;								// #blocks
	int L = int(M/N);						// bock's length
	vector<double> I;						// vectors to fill with block's mean and squared block's mean
	vector<double> I2;
	vector<double> I_new;
	vector<double> I2_new;
	double sum1, sum2, a, sum_prog=0, su2_prog=0, sum_new=0, su2_new=0, err_prog=0;		// useful vars

	// opening output file   	
	ofstream WriteResults;
	WriteResults.open("es02_1.txt");
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open es02_1.txt!" << endl;
		return -1;
	}
  
  // filling vectors
	for (int i=0; i<N; i++){
		sum1 = 0;
		sum2 = 0;
		for (int k=0; k<L; k++){
			sum1 += M_PI/2*cos(M_PI/2*rnd.Rannyu());			// uniform sampling
			a = rnd.P_x();										// generating rnd numbers with new pdf
			sum2 += M_PI/2*cos(M_PI/2*a) / (2*(1-a));			// importance sampling
		}
		sum1 /= L;
		I.push_back(sum1);										// loading vectors
		I2.push_back(sum1*sum1);
		sum2 /= L;
		I_new.push_back(sum2);
		I2_new.push_back(sum2*sum2);
	}

	// progressive means and errors
	for(int i=0; i<N; i++){
		// uniform sampling
		sum_prog += I[i];																	// accumulating means
		su2_prog += I2[i];
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);								// and errors
		cout << (i+1)*L << " " << sum_prog/(i+1) << " " << err_prog;
		WriteResults << (i+1)*L << " " << sum_prog/(i+1) - 1 << " " << err_prog;			// wirting output
		// importance sampling
		sum_new += I_new[i];																// accumulating means 
		su2_new += I2_new[i];
		err_prog = error(sum_new/(i+1), su2_new/(i+1), i);									// and error
		cout << " " << sum_new/(i+1) << " " << err_prog << endl;
		WriteResults << " " << sum_new/(i+1) - 1 << " " << err_prog << endl;				// writting output
	}
	
	WriteResults.close();
	
	return 0;
}

