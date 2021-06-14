#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <memory>
#include <numeric>
#include <iomanip>
#include "random.h"
#include "blocking.h"
#include "MCMC.h"

using namespace std;

int main (int argc, char *argv[]){
	
	int nthrows = 1000000;
	int nblocks = 100;

	// output filename
	string filename = "es08_2.txt";
	
	//taking mu and sigma in input to compute the energy for a grid of values using a script
	double mu = atof(argv[1]);
	double sigma = atof(argv[2]);
	double en, err_en;

	ofstream out("params.out", ios::app);
	int wd = 16;

	//initializing pdf to sample and integrand function
	shared_ptr<Basic_pdf> modq_psi = make_shared<wavefunct>(mu, sigma);
	shared_ptr<Basic_pdf> int_psi = make_shared<integrand>(mu, sigma);

	//initializing MCMC object
	vector<double> x = {0.};
	shared_ptr<MCMC> DistrSample = make_shared<MCMC>(modq_psi, x, atof(argv[3]), "modq_psi.txt" );

	//initializing object to compute progressive mean
	shared_ptr<BasicBlockMean> myBlockMean = make_shared<MCMCIntegral>( DistrSample, int_psi );
	
	PrintProg( nthrows, nblocks, myBlockMean, filename, en, err_en );		// computing and printing progressive means
	
	out << setw(wd) << mu << setw(wd) << sigma << setw(wd) << DistrSample->m_step << setw(wd) << en << setw(wd) << err_en << setw(wd) << DistrSample->m_count/(double)nthrows << endl;

	out.close();
	return 0;
}


