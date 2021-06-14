#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <numeric>
#include "random.h"

using namespace std;

double error(double ave, double av2, int n){
	if (n == 0)
		return 0;
	else
		return sqrt((av2 - pow(ave,2))/n);
}
 
int main (int argc, char *argv[]){

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


	int M = 100000;
	int N = 100;
	int L = int(M/N);
	vector<double> ave;
	vector<double> av2;
	vector<double> var;
	vector<double> va2;
   	double sum1, sum2, a, sum_prog=0, su2_prog=0, sum_var_prog = 0, su2_var_prog = 0, err_prog=0;
   	
   	ofstream WriteResults;
   
   	for (int i=0; i<N; i++){
		sum1 = 0;
		sum2 = 0;
		for (int k=0; k<L; k++){
			a = rnd.Rannyu();
			sum1 += a;
			sum2 += pow( (a-0.5),2 );
		}
		sum1 /= L;
		ave.push_back(sum1);
		av2.push_back(sum1*sum1);
		sum2 /= L;
		var.push_back(sum2);
		va2.push_back(sum2*sum2);
	}
	
	WriteResults.open("ES01_1.txt");
   	if(!WriteResults){
   		cerr << "PROBLEM: Unable to open ES01_1.txt!" << endl;
   		return -1;
   	}
	for(int i=0; i<N; i++){
		sum_prog += ave[i];
		su2_prog += av2[i];
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);
		cout << (i+1)*L << " " << sum_prog/(i+1) << " " << err_prog;
		WriteResults << (i+1)*L << " " << sum_prog/(i+1) - 0.5 << " " << err_prog;
		sum_var_prog += var[i];
		su2_var_prog += va2[i];
		err_prog = error(sum_var_prog/(i+1), su2_var_prog/(i+1), i);
		cout << " " << sum_var_prog/(i+1) << " " << err_prog << endl;
		WriteResults << " " << sum_var_prog/(i+1) - 1/(double)12 << " " << err_prog << endl;	
	}
	
	int dim[M];
	vector<double> chi;
	vector<double> chi2;
	double chi_prog = 0;
	int j_max = 100;	
	M = 100;
	N = 10000;
	ex = int(N/M);
	sum1=0;
	
	for(int j=0; j<j_max; j++){
		sum = 0;
		for(int i=0; i<N; i++){
			a = rnd.Rannyu();
			for(int k=0;k<M;k++){
				if( int(a*M) == k )
					dim[k]++;				
			}
		}
		for(int k=0; k<M; k++)
			sum1 += pow( (dim[k]-ex), 2) / ex;		
		chi.push_bak(sum1);
		chi2.push_back(sum1*sum1);
	}
	for( int i=0; i<chi.size(); i++){
		chi_prog += chi[i];
		chi2_prog += chi2[i];
		err_prog = error( chi2/(i+1), chi/(i+1), i) 
		cout << (i+1)*N << " " << chi_prog/(i+1) << " " << err_prog << endl;
		WriteResults << (i+1)*N << " " << chi_prog/(i+1) - 100 << " " << err_prog << endl;
	} 

	WriteResults.close();
	return 0;
}

