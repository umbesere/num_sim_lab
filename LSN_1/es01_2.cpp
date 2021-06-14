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
	
	//** TESTING GAUSSIAN, EXPONENTIAL AND CAUCHY GENERATORS **//
	
	// initializing variables
	int M = 10000;														// #points for each histogram
	int dim[4] = {1, 2, 10, 100};							// #points used for the mean
	double sum1 = 0, sum2 = 0, sum3 = 0;			// useful sums
	
	// opening output file
	ofstream WriteResults;
	WriteResults.open("es01_2.txt");
	if(!WriteResults){
		cerr << "PROBLEM: Unable to open es01_2.txt!" << endl;
		return -1;
	}
	
	// generating points
	for(int i = 0; i < 4; i++){
		for(int k = 0; k < M; k++){
			sum1 = 0; sum2 = 0; sum3 = 0;
			for(int j = 0; j < dim[i]; j++){
				// computing means using variious generators
				sum1 += rnd.Rannyu() / dim[i];
				sum2 += rnd.Exponential(1.) / dim[i];
				sum3 += rnd.Cauchy(0.,1.) / dim[i];
			}
			WriteResults << sum1 << " " << sum2 << " " << sum3 << endl;			// writing output
		}
	}
	
	WriteResults.close();
	return 0;
}
