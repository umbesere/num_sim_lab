#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <iomanip>
#include "random.h"
#include "genetic.h"

using namespace std;

int main (int argc, char *argv[]){
	
	//SIMULATED ANNEALING//
	// output filename
	string filename = "es10_1.txt";
	//setting up the temperature schedule
	int ncities = 32;
	int niter = 100;
	int ntemp = 10000;
	double temp = 10.001;
	double deltat = 0.001;
	Random rnd("seed.in","seed.out");

	//initializing path
	double length, newlength;
	vector<int> path = StartingInd( ncities, rnd);
	vector<int> newpath;
	//Check_Individual( path );
	
	//initializing map
	vector<vector<double>> map = Square( ncities, rnd );
	ofstream outcities;
	outcities.open("map.txt");
	for ( auto city : map )
		outcities << city[0] << " " << city[1] << endl;
	outcities.close();

	//setting up the output
	ofstream out;
	out.open(filename);
	int wd = 18;
	out << setw(wd) << "ntotiter" << endl;
	out << setw(wd) << niter*ntemp << endl;
	out << setw(wd) << "Niter" << setw(wd) << "Route_length" << endl;

	//temperature routine
	for(int j=0; j<ntemp; j++){
		temp -= deltat; //decreases temperature of deltat
		for (int i=0; i<niter; i++){	//iterates at fixed temp
			// Swap
			newpath = Mutation_swap( path, rnd );
			length = IndFitting(path, map);
			newlength = IndFitting(newpath, map);
			if ( newlength > length ){
				if ( rnd.Rannyu() < exp(-(newlength-length)/temp) )
					path = newpath;
			} else {
				path = newpath;
			}
			//Check_Individual( path );
			// Shift
			newpath = Mutation_shift( path, rnd );
			length = IndFitting(path, map);
			newlength = IndFitting(newpath, map);
			if ( newlength > length ){
				if ( rnd.Rannyu() < exp(-(newlength-length)/temp) )
					path = newpath;
			} else {
				path = newpath;
			}
			//Check_Individual( path );
			// Swapcont
			newpath = Mutation_swapcont( path, rnd );
			length = IndFitting(path, map);
			newlength = IndFitting(newpath, map);
			if ( newlength > length ){
				if ( rnd.Rannyu() < exp(-(newlength-length)/temp) )
					path = newpath;
			} else {
				path = newpath;
			}
			//Check_Individual( path );
			// Inversion
			newpath = Mutation_inversion( path, rnd );
			length = IndFitting(path, map);
			newlength = IndFitting(newpath, map);
			if ( newlength > length ){
				if ( rnd.Rannyu() < exp(-(newlength-length)/temp) )
					path = newpath;
			} else {
				path = newpath;
			}
			//Check_Individual( path );
			cout << "temp: " << temp << ", " << "iter " << niter*j+i << "/" << niter*ntemp << endl;
			out << setw(wd) << niter*j+i << setw(wd) << IndFitting( path, map ) << endl;
		}
	}
	
	//printing results
	for( auto x : path )
		out << setw(wd/3) << x;
	out << endl;
	
	out.close();	
	return 0;
}
