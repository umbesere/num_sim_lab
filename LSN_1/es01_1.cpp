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

	//** PROGRESSIVE MEANS AND STANDARD DEVIATIONS **//
	
	// initializing variables
	int M = 100000; 				// #throws
	int N = 100;						// #blocks
	int L = int(M/N);				// block's length
	vector<double> ave;			// mean, mean squared, sigma , sigma squared vectors
	vector<double> av2;
	vector<double> var;
	vector<double> va2;
 	double sum1, sum2, a, sum_prog=0, su2_prog=0, sum_var_prog = 0, su2_var_prog = 0, err_prog=0; 		// useful variables
   
	ofstream WriteResults;	
	
  // loading vectors
	for (int i=0; i<N; i++){
		sum1 = 0;
		sum2 = 0;
		for (int k=0; k<L; k++){			// computing block's mean
			a = rnd.Rannyu();
			sum1 += a;
			sum2 += pow( (a-0.5),2 );
		}
		sum1 /= L;
		ave.push_back(sum1);					// loading averages
		av2.push_back(sum1*sum1);			// loading squared averages
		sum2 /= L;										
		var.push_back(sum2);					// loading dev std
		va2.push_back(sum2*sum2);			// loading squared dev std
	}
	
	// opening file
	WriteResults.open("ES01_1.txt");
  	if(!WriteResults){
  		cerr << "PROBLEM: Unable to open ES01_1.txt!" << endl;
  		return -1;
	}
	
	// computing progressive averages with uncertainty and output
	for(int i=0; i<N; i++){
		sum_prog += ave[i];																					// progressive sums
		su2_prog += av2[i];
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);												// and errors
		cout << (i+1)*L << " " << sum_prog/(i+1) << " " << err_prog;
		WriteResults << (i+1)*L << " " << sum_prog/(i+1) - 0.5 << " " << err_prog;							//	writing progressive means with errors
		sum_var_prog += var[i];																				// progressive std dev sums
		su2_var_prog += va2[i];
		err_prog = error(sum_var_prog/(i+1), su2_var_prog/(i+1), i);										// and errors
		cout << " " << sum_var_prog/(i+1) << " " << err_prog << endl;																
		WriteResults << " " << sum_var_prog/(i+1) - 1/(double)12 << " " << err_prog << endl;				// writing progressive std devs with errors
	}
	
	cout << endl;
	
	//** CHI^2 **//
	
	// initializing variables
	int rep = 100;											// #chi2
	double chi = 0;											// initializie chi	
	M = 100;														// #intervals
	N = 10000;													// #points
	double expect = int(N/M);						// expected #points in an interval
	vector<int> obs(M, 0.);							// useful vector to register interval's occupation	
	
	// computing chi2
	for(int j=0; j<rep; j++){																
		fill( obs.begin(), obs.end(), 0 );		
		chi = 0;
		for(int i=0; i<N; i++){				// filling obs
			a = rnd.Rannyu();
			for(int k=0;k<M;k++){
				if( int(a*M) == k )
					obs[k]++;				
			}
		}
		for(int k=0; k<M; k++)														
			chi += pow( (obs[k]-expect), 2) / expect;						// computing chi2
		cout << j << " " << chi << endl;
		WriteResults << j << " " << chi - 100 << endl;				// writing chi2
	}
		

	WriteResults.close();
	return 0;
}

