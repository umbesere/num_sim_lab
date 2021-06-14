#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<memory>
#include <iomanip>
#include "random.h"
#include "MCMC.h"
#include "blocking.h"

using namespace std;

// error function to compute error from mean and mean squared
double error(double ave, double av2, int n){
	if (n == 0)
		return 0;
	else
		return sqrt((av2 - pow(ave,2))/(n+1)); //OCCHIO, DEVI CORREGGERE N -> N+1 IN TUTTI I FILE !!!!!
}

// destructor
BasicBlockMean::~BasicBlockMean(){};

double MCMCIntegral :: BlockMean( int L ){			//Uncomment to have a code that save the sampled distribution
	ofstream out(m_MCMC->m_name, ios::app);			//<---
	int wd = 18;										//<--
	double sum = 0;
	for (int i =0; i<L; i++){
		m_MCMC->UnifStep();
		sum += m_integrand->Eval(m_MCMC->getPos());
		if (i % 5 == 0)								//<--
			out << setw(wd) << m_MCMC->getX() << endl;	//<--
	}
	out.close();										//<--
	return sum/L;
};		// key function

/* FUNCTIONS TO COMPUTE AND PRINT MEANS VIA BLOCKING METHODS */
// NB ... TO TEST!!!!
void PrintProg( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int wd = 16;
	int L = int(nthrows/nblock);
	double sum = 0., sum_prog = 0., su2_prog = 0., err_prog = 0.;
	WriteResults << "Structure of the OUTPUT:" << endl;
	WriteResults << setw(wd) << "Throws" << setw(wd) << "Block_mean" << setw(wd) << "Prog_mean" << setw(wd) << "Prog_err" << endl;
	for (int i=0; i<nblock; i++){
		sum = f->BlockMean( L );
		sum_prog += sum;
		su2_prog += sum*sum;
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);
		WriteResults << setw(wd) << (i+1)*L << setw(wd) << sum << setw(wd) << sum_prog/(i+1) << setw(wd) << err_prog << endl;
	}
	WriteResults.close();
};

void PrintProg( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename, double &mean, double &rms){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int wd = 16;
	int L = int(nthrows/nblock);
	double sum = 0., sum_prog = 0., su2_prog = 0., err_prog = 0.;
	WriteResults << "Structure of the OUTPUT:" << endl;
	WriteResults << setw(wd) << "Throws" << setw(wd) << "Block_mean" << setw(wd) << "Prog_mean" << setw(wd) << "Prog_err" << endl;
	for (int i=0; i<nblock; i++){
		sum = f->BlockMean( L );
		sum_prog += sum;
		su2_prog += sum*sum;
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);
		WriteResults << setw(wd) << (i+1)*L << setw(wd) << sum << setw(wd) << sum_prog/(i+1) << setw(wd) << err_prog << endl;
	}
	WriteResults.close();
	mean = sum_prog/nblock; rms = err_prog;
};

void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int wd = 16;
	int size = f.size();
	int L = int(nthrows/nblock);
	vector<double> sum(size, 0.), sum_prog(size, 0.), su2_prog(size, 0.), err_prog(size, 0.);
	WriteResults << "Structure of the OUTPUT:" << endl;
	WriteResults << setw(wd) << "Throws";
	for(int j=0; j<size; j++)
		WriteResults << setw(wd) << "Block_mean" << j << setw(wd) << "Prog_mean" << j << setw(wd) << "Prog_err" << j << endl;	
	for (int i=0; i<nblock; i++){
		WriteResults << setw(wd) <<(i+1)*L;
		for(int j=0; j<size; j++){
			sum[j] = f[j]->BlockMean( L );
			sum_prog[j] += sum[j];
			su2_prog[j] += sum[j]*sum[j];
			err_prog[j] = error(sum_prog[j]/(i+1), su2_prog[j]/(i+1), i);
			WriteResults << setw(wd) << sum[j] << setw(wd) << sum_prog[j]/(i+1) << setw(wd) << err_prog[j];
		}
		WriteResults << endl;
	}
	WriteResults.close();
};

void Print( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int wd = 16;
	int L = int(nthrows/nblock);
	double sum = 0., sum_prog = 0., su2_prog = 0., err_prog = 0.;
	WriteResults << "Structure of the OUTPUT:" << endl;
	WriteResults << setw(wd) << "Throws" << setw(wd) << "Block_mean" << setw(wd) << "Prog_mean" << setw(wd) << "Prog_err" << endl;
	for (int i=0; i<nblock; i++){
		sum = f->BlockMean( L );
		sum_prog += sum;
		su2_prog += sum*sum;
	}
	err_prog = error(sum_prog/nblock, su2_prog/nblock, nblock-1);
	WriteResults << setw(wd) << sum << setw(wd) << sum_prog/nblock << setw(wd) << err_prog << endl;
	WriteResults.close();
}

void Print( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int wd = 16;
	int size = f.size();
	int L = int(nthrows/nblock);
	vector<double> sum(size, 0.), sum_prog(size, 0.), su2_prog(size, 0.), err_prog(size, 0.);
	WriteResults << "Structure of the OUTPUT:" << endl;
	WriteResults << setw(wd) << "Throws";
	for(int j=0; j<size; j++)
		WriteResults << setw(wd) << "Prog_mean" << j << setw(wd) << "Prog_err" << j << endl;	
	for (int i=0; i<nblock; i++){
		WriteResults << setw(wd) <<(i+1)*L;
		for(int j=0; j<size; j++){
			sum[j] = f[j]->BlockMean( L );
			sum_prog[j] += sum[j];
			su2_prog[j] += sum[j]*sum[j];
			WriteResults << setw(wd) << sum[j] << setw(wd) << sum_prog[j]/(i+1) << setw(wd) << err_prog[j];
		}
	}
	for(int j=0; j<size; j++){
		err_prog[j] = error(sum_prog[j]/nblock, su2_prog[j]/nblock, nblock-1);
		WriteResults << setw(wd) << sum_prog[j]/nblock << setw(wd) << err_prog[j];
	}
	WriteResults << endl;
	WriteResults.close();
};


/*
double Mean( int nthrows, shared_ptr<BasicBlockMean> &f){
	double sum = 0.;
	for (int i=0; i<nthrows; i++)
		sum += f->BlockMean( 1 );
	return sum/nthrows;
};
*/