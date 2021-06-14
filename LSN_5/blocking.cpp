#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include<memory>
#include "random.h"
#include "MCMC.h"
#include "blocking.h"

using namespace std;

// error function to compute error from mean and mean squared
double error(double ave, double av2, int n){
	if (n == 0)
		return 0;
	else
		return sqrt((av2 - pow(ave,2))/n);
}

// destructor
BasicBlockMean::~BasicBlockMean(){};

double UnifMCMC :: BlockMean( int L ){
	ofstream out;
	out.open(m_MCMC->m_name, ios::app); //to print out the positions of the MCMC
	double sum = 0;
	for (int i =0; i<L; i++){
		m_MCMC->UnifStep( *m_rnd );
		sum += m_MCMC->getR();
		if (m_MCMC->m_count % 250 == 0)	//every ~500 steps
			out << m_MCMC->getX() << " " << m_MCMC->getY() << " " << m_MCMC->getZ() << endl;
	}
	out.close();
	return sum/L;
};		// key function

double GaussMCMC :: BlockMean( int L ){
	ofstream out;
	out.open(m_MCMC->m_name, ios::app);
	double sum = 0;
	for (int i =0; i<L; i++){
		m_MCMC->GaussStep( *m_rnd );
		sum += m_MCMC->getR();
		if (m_MCMC->m_count % 250 == 0)
			out << m_MCMC->getX() << " " << m_MCMC->getY() << " " << m_MCMC->getZ() << endl;
	}
	out.close();
	return sum/L;
};		// key function



/* FUNCTIONS TO COMPUTE AND PRINT MEANS VIA BLOCKING METHODS */
// NB ... TO TEST!!!!
void Print( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults ){
	int size = f.size();
	int L = int(nthrows/nblock);
	vector<double> sum(size, 0.);
	for (int i=0; i<nblock; i++){
		for(int j=0; j<size; j++){
			sum[j] += f[j]->BlockMean( L );
		}
	}
	WriteResults << nthrows << " ";
	for(int j=0; j<size; j++)
		WriteResults << sum[j]/nblock << " " << error(sum[j]/nblock, sum[j]*sum[j]/nblock, nblock) << " ";
	WriteResults << endl;
}

void Print( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
	int L = int(nthrows/nblock);
	double sum = 0., sum_prog = 0., su2_prog = 0., err_prog = 0.;
	for (int i=0; i<nblock; i++){
		sum = f->BlockMean( L );
		sum_prog += sum;
		su2_prog += sum*sum;
		err_prog = error(sum_prog/(i+1), su2_prog/(i+1), i);
		WriteResults << sum_prog/(i+1) << " " << err_prog;
	}
	WriteResults.close();
}

void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, string filename){
	ofstream WriteResults;
	WriteResults.open(filename);
		if(!WriteResults){
		cerr << "PROBLEM: Unable to open " << filename << "!" << endl;
	}	
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
	WriteResults.close();
};

/*  USEFUL ??
void BlockPrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults ){
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

