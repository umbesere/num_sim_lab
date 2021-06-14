#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<memory>
#include "blocking.h"
#include "brownian_motion.h"
#include "random.h"

using namespace std;

// error function to compute error from mean and mean squared
double error(double ave, double av2, int n){
	if (n == 0)
		return 0;
	else
		return sqrt((av2 - pow(ave,2))/n);
}

// destructor
BasicBlockMean::~BasicBlockMean(){}

/* METHODS TO SIMULATE THE PRRICE OF CALL/PUT VANILLA EU OPTIONS */

// Call with single step
double StepCall :: BlockMean( int L ){
	double sum = 0;
	for (int k=0; k<L; k++){
			m_BM->SetPos();								// to reset GBM position
			m_BM->Gstep( m_rnd, m_time );				// geometric step
			sum += max( 0., exp(-m_BM->GetDrift()*m_time)*(m_BM->GetPos()-m_strike_price));		// single evaluation of call's price
		}
	return sum/L;		// mean of a single block
}

// Call with multiple steps
double PathCall :: BlockMean( int L ){
	double sum = 0;
	for (int k=0; k<L; k++){
			m_BM->SetPos();								// to reset GBM position
			m_BM->Gpath( m_rnd, m_time, m_nstep );		// #nstep of geometric steps
			sum += max( 0., exp(-m_BM->GetDrift()*m_time)*(m_BM->GetPos()-m_strike_price) );	// single evaluation of call's price
		}
	return sum/L;		 // mean of a single block
}

// Put with single step
double StepPut :: BlockMean( int L ){
	double sum = 0;
	for (int k=0; k<L; k++){
			m_BM->SetPos();								// to reset GBM position
			m_BM->Gstep( m_rnd, m_time );				// geometric step
			sum += max( 0., exp(-m_BM->GetDrift()*m_time)*(-m_BM->GetPos()+m_strike_price) );	// single evaluation of put's price
		}
	return sum/L;		// mean of a single block
}

// Put with multiple steps
double PathPut :: BlockMean( int L ){
	double sum = 0;
	for (int k=0; k<L; k++){
			m_BM->SetPos();								// to reset GBM position
			m_BM->Gpath( m_rnd, m_time, m_nstep );		// #nstep of geometric steps
			sum += max( 0., exp(-m_BM->GetDrift()*m_time)*(-m_BM->GetPos()+m_strike_price) );	// single evaluation of put's price
		}
	return sum/L;		// mean of a single block
}



/* FUNCTIONS TO COMPUTE AND PRINT MEANS VIA BLOCKING METHODS */

// NB ... TO TEST!!!!
// computing and printing final mean and error on an already open file (multiple process version)
void Print( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults ){
	int size = f.size();
	int L = int(nthrows/nblock);
	double mean;
	vector<double> sum(size, 0.), su2(size, 0.);
	for (int i=0; i<nblock; i++){
		for(int j=0; j<size; j++){
			mean = f[j]->BlockMean( L );
			sum[j] += mean;
			su2[j] += mean*mean; 
		}
	}
	WriteResults << nthrows << " ";
	for(int j=0; j<size; j++)
		WriteResults << sum[j]/nblock << " " << error(sum[j]/nblock, su2[j]/nblock, nblock) << " ";
	WriteResults << endl;
}

// computing and printing progressive mean and error on a new file (single process version)
void PrintProg( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int L = int(nthrows/nblock);
	double mean = 0., sum_prog = 0., su2_prog = 0., err_prog = 0.;
	for (int i=0; i<nblock; i++){
		mean = f->BlockMean( L );
		sum_prog += mean;
		su2_prog += mean*mean;
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);
		WriteResults << sum_prog/(i+1) << " " << err_prog;
	}
	WriteResults.close();
}

// computing and printing progressive means and errors on a new file (multiple processes version)
void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int size = f.size();
	double mean;
	int L = int(nthrows/nblock);
	vector<double> sum_prog(size, 0.), su2_prog(size, 0.), err_prog(size, 0.);
	for (int i=0; i<nblock; i++){
		WriteResults << (i+1)*L;
		for(int j=0; j<size; j++){
			mean = f[j]->BlockMean( L );
			sum_prog[j] += mean;
			su2_prog[j] += mean*mean;
			err_prog[j] = error(sum_prog[j]/(i+1), su2_prog[j]/(i+1), i);
			WriteResults << " " << sum_prog[j]/(i+1) << " " << err_prog[j];
		}
		WriteResults << endl;
	}
	WriteResults.close();
};

/*  USEFUL ??
void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults ){
	int size = f.size();
	int L = int(nthrows/nblock);
	vector<double> sum(size, 0.), sum_prog(size, 0.), su2_prog(size, 0.), err_prog(size, 0.);
	for (int i=0; i<nblock; i++){
		WriteResults << (i+1)*L;
		for(int j=0; j<size; j++){
			sum[j] = f[j]->BlockMean( L );
			sum_prog[j] += sum[j];
			su2_prog[j] += sum[j]*sum[j];
			err_prog[j] = error(sum_prog[j]/(i+1), su2_prog[j]/(i+1), i);
			WriteResults << " " << sum_prog[j]/(i+1) << " " << err_prog[j];
		}
		WriteResults << endl;
	}
}
*/

