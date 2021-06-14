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
		return sqrt((av2 - pow(ave,2))/(n+1));
}
 
int main (int argc, char *argv[]){

	// initializing random class
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
	
	//** SIMULATION OF BUFFFON'S EXPERIMENT **//
	
	int M = 1000000, N = 100;																				// #throws and #blocks
	int N_hit = 0;																									// #hitting neelds
	double L = 4, d = 5;																						// experiment's paramenters
	int n = int(M/N);																								// block's length
	double z, theta, a, sum_prog = 0, su2_prog = 0,err_prog;				
	vector<double> pi_est;																					// pi and pi squared vectors
	vector<double> pi2_est;
	
	// experiment's simulation
	for (int k = 0; k < N; k++){
		N_hit = 0;
		for (int i = 0; i < n; i++){
			z = rnd.Rannyu(0., d/2);				// throwing a needle
			theta = rnd.Angle()/4;					
			if( z - L/2*cos(theta) < 0 )		// condition for hitting a line
				N_hit++;	
		}
		a = 2*L*n/(N_hit*d);
		pi_est.push_back( a );						// loading pi and pi squared vectors
		pi2_est.push_back( a*a );
	}
	
	// opening output file
	ofstream WriteResults;
	WriteResults.open("es01_3.txt");
	if(!WriteResults){
		cerr << "PROBLEM: Unable to open es01_3.txt!" << endl;
		return -1;
	}
	
	// computing pi mean and wiritng output
	for(int i=0; i<N; i++){
		sum_prog += pi_est[i];																							// progressive means
		su2_prog += pi2_est[i];
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);								// and errors
		cout << (i+1)*n << " " << sum_prog/(i+1) << " " << err_prog << endl;
		WriteResults << (i+1)*n << " " << sum_prog/(i+1) - M_PI << " " << err_prog << endl;		// writing output
	}
		
	WriteResults.close();
	return 0;
}
