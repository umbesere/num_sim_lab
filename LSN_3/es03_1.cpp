#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <memory>
#include <numeric>
#include "random.h"
#include "brownian_motion.h"
#include "blocking.h"

using namespace std;

int main (int argc, char *argv[]){
	
	//** SIMULATION OF VANILLA EU OPTION'S PRICE **//

	// output filename
	string filename = "es03_1.txt";
	
	// blocking vars
	int nthrows = 100000;
	int nblock = 100;
	int nstep = 100;
	
	// option params
	double in_asset_price = 100;
	double delivery_time = 1;
	double int_rate = 0.1;
	double volatility = 0.25;
	double strike_price = 100;


	vector<shared_ptr<BasicBlockMean>>  myBlockMean(4);		// vector of pointer to BasicBlockMean 
	myBlockMean[0] = make_shared<StepCall>(in_asset_price, int_rate, volatility, delivery_time, strike_price);				// Call
	myBlockMean[1] = make_shared<StepPut>(in_asset_price, int_rate, volatility, delivery_time, strike_price);				// Put
	myBlockMean[2] = make_shared<PathCall>(in_asset_price, int_rate, volatility, delivery_time, strike_price, nstep);		// Descrete Call
	myBlockMean[3] = make_shared<PathPut>(in_asset_price, int_rate, volatility, delivery_time, strike_price, nstep);		// Descrete Put
	
	PrintProg( nthrows, nblock, myBlockMean, filename );		// computing and printing progressive means 
	
	return 0;
}
