#ifndef __Blocking__
#define __Blocking__

#include "MCMC.h"
#include "random.h"
#include <cstdlib>
#include <vector>
#include <memory>

using namespace std;

// error function
double error(double ave, double av2, int n);

// Base class -> BasicBlockMean
// Abstract class from which one can create his own derived class to solve a specific problem
class BasicBlockMean {
	public:
		virtual ~BasicBlockMean() = 0;				// virtual destructor
		virtual double BlockMean( int L ) = 0;		// key function
		/*
			BlockMean:
			It's a function that calculate the mean value of a certain process on a single block
			as a function of the block length L. It is a virtual function that can be overloaded
			in a derived class to compute a specific simulation/solve a particular problem.		
		*/
};

class MCMCIntegral : public BasicBlockMean {

	public:
		MCMCIntegral( shared_ptr<MCMC> myMCMC, shared_ptr<Basic_pdf> myint){
			m_MCMC = myMCMC;
			m_integrand = myint;
		};
		~MCMCIntegral(){};					// virtual destructor
		double BlockMean( int L );		// key function
	private:
		shared_ptr<MCMC> m_MCMC;
		shared_ptr<Basic_pdf> m_integrand;
};

// computing and calculating progressive mean and error
void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, string filename);		// multiple progmean, new file
void PrintProg( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename);				// single progmean, new file
void PrintProg( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename, double &mean, double &rms);
void Print( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults);	// multiple finalmean, open file
void Print( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, ofstream &WriteResults);	// multiple finalmean, open file
//double Mean( int nthrows, shared_ptr<BasicBlockMean> &f);
#endif //__Blocking__