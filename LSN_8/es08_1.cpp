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
	int nblocks = 100;

	// output filename
	string filename = "es08_1.txt";
	
	double mu = 1.;
	double sigma = 0.5;

	//initializing pdf to sample and integrand function
	shared_ptr<Basic_pdf> modq_psi = make_shared<wavefunct>(mu, sigma);
	shared_ptr<Basic_pdf> int_psi = make_shared<integrand>(mu, sigma);

	//initializing MCMC object
	vector<double> x = {0.};
	shared_ptr<MCMC> DistrSample = make_shared<MCMC>(modq_psi, x, 2.5, "modq_psi.txt" );
	//DistrSample->Stabilize(100000, 0.1, 0.05);

	//initializing object to compute progressive means
	vector<shared_ptr<BasicBlockMean>>  myBlockMean(1);
	myBlockMean[0] = make_shared<MCMCIntegral>( DistrSample, int_psi );

	PrintProg( nthrows, nblocks, myBlockMean, filename );		// computing and printing progressive means 
	cout << "Acceptance rate:\t" << DistrSample->m_count/(double)nthrows << endl;
	
	return 0;
}
