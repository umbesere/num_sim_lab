#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "random.h"
#include "random_walk.h"

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
   	
	//** SIMULATION OF DESCRETE AND CONTINUOUS RANDOM WALKS **//   	

	// initializing variables
	int M = 10000;											// #RWs
	int N = 100;												// #blocks
	int L = int(M/N);										// block's length
	int N_step = 100;										// total number of steps
	double sum_discr, mean_discr, mea2_discr, err_discr, sum_cont, mean_cont, mea2_cont, err_cont;	
	vector<double> r2_discr(N);					// vectors to be filled by block's mean of |r|^2 for descrete/continuous RWs
	vector<double> r22_discr(N);
	vector<double> r2_cont(N);
	vector<double> r22_cont(N);
	
	// initializing M random walks
	RandomWalk RW_discr[M];
	RandomWalk RW_cont[M];
  
  	// setting initial points to origin
  	for(int i=0; i<M; i++){
  		RW_discr[i].SetCoord();
  		RW_cont[i].SetCoord();
  	}
  
  // opening output file
	ofstream WriteResults;
	WriteResults.open("es02_2.txt");
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open es02_2.txt!" << endl;
		return -1;
	}

	// simulating RWs	
	for (int j=0; j<N_step; j++){
		for(int i=0; i<N; i++){
			sum_discr = 0;
			sum_cont = 0;
			// simulating a step for each RW of the considered block
			for(int k=L*i; k<L*(i+1); k++){
				RW_discr[k].DiscrStep(rnd);
				RW_cont[k].ContStep(rnd);	
				sum_discr += pow( RW_discr[k].getR(), 2);
				sum_cont += pow( RW_cont[k].getR(), 2);
			}
			// loading vectors with means and squared means
			sum_discr /= L;
			sum_cont /= L;
			r2_discr[i] = sum_discr;
			r22_discr[i] = sum_discr*sum_discr;
			r2_cont[i] = sum_cont;
			r22_cont[i] = sum_cont*sum_cont;
		}
		// computing progressive means and errors
		mean_discr = accumulate(r2_discr.begin(), r2_discr.end(), 0.0) / N;
		mea2_discr = accumulate(r22_discr.begin(), r22_discr.end(), 0.0) / N;
		mean_cont = accumulate(r2_cont.begin(), r2_cont.end(), 0.0) / N;
		mea2_cont = accumulate(r22_cont.begin(), r22_cont.end(), 0.0) / N;
		err_discr = error(mean_discr, mea2_discr, N);
		err_cont = error(mean_cont, mea2_cont, N);
		// writing results	
		cout << j << sqrt(mean_discr) << " " << err_discr/(2*sqrt(mean_discr)) << " " << sqrt(mean_cont) << " " << err_cont/(2*sqrt(mean_cont)) << endl;
		WriteResults << j << " " << sqrt(mean_discr) << " " << err_discr/(2*sqrt(mean_discr)) << " " << sqrt(mean_cont) << " " << err_cont/(2*sqrt(mean_cont)) << endl;
	}
	
	WriteResults.close();	
	return 0;
}

