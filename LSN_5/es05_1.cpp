#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <memory>
#include <numeric>
#include "random.h"
#include "blocking.h"
#include "MCMC.h"

using namespace std;

int main (int argc, char *argv[]){
	
	int nthrows = 1000000;
	int nblock = 100;

	// output filename
	string filename = "es05_1.txt";
	
	shared_ptr<Basic_pdf> psi1 = make_shared<psi100>();		//initializing functions psi100
	shared_ptr<Basic_pdf> psi2 = make_shared<psi210>();		//initializing functions psi210

	shared_ptr<MCMC> UnifPsi100 = make_shared<MCMC>(psi1, 0., 0., 0., 1.25, "UnifPsi100.txt" );	//Markov Chain Monte Carlo for psi100 with Metropolis uniform step
	shared_ptr<MCMC> GaussPsi100 = make_shared<MCMC>(psi1, 0., 0., 0., .75, "GaussPsi100.txt");	//Markov Chain Monte Carlo for psi100 with Metropolis gaussian step
	shared_ptr<MCMC> UnifPsi210 = make_shared<MCMC>(psi2, 0., 0., 0., 3., "UnifPsi210.txt");	//Markov Chain Monte Carlo for psi210 with Metropolis uniform step
	shared_ptr<MCMC> GaussPsi210 = make_shared<MCMC>(psi2, 0., 0., 0., 1.9, "GaussPsi210.txt");	//Markov Chain Monte Carlo for psi210 with Metropolis gaussian step

	vector<shared_ptr<BasicBlockMean>>  myBlockMean(4);		//objects for progressive means computation
	myBlockMean[0] = make_shared<UnifMCMC>( UnifPsi100 );
	myBlockMean[1] = make_shared<GaussMCMC>( GaussPsi100 );
	myBlockMean[2] = make_shared<UnifMCMC>( UnifPsi210 );
	myBlockMean[3] = make_shared<GaussMCMC>( GaussPsi210 );
	PrintProg( nthrows, nblock, myBlockMean, filename );		// computing and printing progressive means 

	cout << "Acceptance rates:" << endl;
	cout << "	Uniform_psi100 = " << UnifPsi100->m_count/(double)nthrows << endl;
	cout << "	Gaussian_psi100 = " << GaussPsi100->m_count/(double)nthrows << endl;
	cout << "	Uniform_psi210 = " << UnifPsi210->m_count/(double)nthrows << endl;
	cout << "	Gaussian_psi210 = " << GaussPsi210->m_count/(double)nthrows << endl;
	
	return 0;
}
