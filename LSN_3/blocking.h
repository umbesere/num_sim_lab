#ifndef __Blocking__
#define __Blocking__

#include "brownian_motion.h"
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

// computing and calculating progressive mean and error
// void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults ); USEFUL??
void PrintProg( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, string filename);		// multiple progmean, new file
void PrintProg( int nthrows, int nblock, shared_ptr<BasicBlockMean> &f, string filename);				// single progmean, new file
void Print( int nthrows, int nblock, vector<shared_ptr<BasicBlockMean>> &f, ofstream &WriteResults);	// multiple finalmean, open file

///** CLASSES TO SIMULATE THE PRICE OF CALL/PUT VANNILA EU OPTIONS **///

// Call with single step
class StepCall : public BasicBlockMean{

	public:
		// constructor
		StepCall( double in_asset_price, double int_rate, double volatility, double time, double strike_price){
			m_rnd = make_shared<Random>("seed.out");										// to gen random numbers
			m_BM = make_shared<BrownianMotion>(in_asset_price, int_rate, volatility);		// to simulate GBM
			m_time = time;
			m_strike_price = strike_price;
		};
		// destructor
		~StepCall(){};
		// overload of BlockMean
		double BlockMean( int L );

	private:
		// params
		shared_ptr<BrownianMotion> m_BM;	
		shared_ptr<Random> m_rnd;
		double m_time;
		double m_strike_price;
};

// descrete Call with multiple steps
class PathCall : public BasicBlockMean{

	public:
		// constructur
		PathCall( double in_asset_price, double int_rate, double volatility, double time, double strike_price, int nstep){
			m_rnd = make_shared<Random>("seed.out");										// to gen random numbers
			m_BM = make_shared<BrownianMotion>(in_asset_price, int_rate, volatility);		// to simulate GBM
			m_time = time;																	
			m_nstep = nstep;																// #steps
			m_strike_price = strike_price;													
		};
		// destructur
		~PathCall(){};
		// overload of BlockMean
		double BlockMean( int L );

	private:
		// params
		shared_ptr<BrownianMotion> m_BM;	
		shared_ptr<Random> m_rnd;
		int m_nstep;
		double m_time;
		double m_strike_price;
};

// Put with single step
class StepPut : public BasicBlockMean{

	public:
		// constructor
		StepPut( double in_asset_price, double int_rate, double volatility, double time, double strike_price){
			m_rnd = make_shared<Random>("seed.out");										// to gen random numbers
			m_BM = make_shared<BrownianMotion>(in_asset_price, int_rate, volatility);		// to simulate GBM
			m_time = time;
			m_strike_price = strike_price;
		};
		// destructur
		~StepPut(){};
		// overload of BlockMean
		double BlockMean( int L );

	private:
		shared_ptr<BrownianMotion> m_BM;	
		shared_ptr<Random> m_rnd;
		double m_time;
		double m_strike_price;
};

// descrete Put with multiple steps 
class PathPut : public BasicBlockMean{

	public:
		// constructor
		PathPut( double in_asset_price, double int_rate, double volatility, double time, double strike_price, int nstep){
			m_rnd = make_shared<Random>("seed.out");										// to gen random numbers
			m_BM = make_shared<BrownianMotion>(in_asset_price, int_rate, volatility);		// to simulate GBM
			m_time = time;
			m_nstep = nstep;																// #steps
			m_strike_price = strike_price;
		};
		// destructur
		~PathPut(){};
		// overload of BlockMean
		double BlockMean( int L );

	private:
		shared_ptr<BrownianMotion> m_BM;
		shared_ptr<Random> m_rnd;
		int m_nstep;
		double m_time;
		double m_strike_price;
};

#endif //__Blocking