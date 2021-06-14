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
	
	// output filename
	string filename = "es09.txt";
	//setting up initial population and map
	int ncities = 32;
	int nroutes = 900;
	int niter = 1000;
	Random rnd("seed.out");
	vector<vector<int>> pop = StartingPop( nroutes, ncities, rnd );
	//Check_Population( pop );
	vector<vector<double>> map = Square( ncities, rnd );
	//printing the map
	ofstream outcities;
	outcities.open("map.txt");
	for ( auto city : map )
		outcities << city[0] << " " << city[1] << endl;
	outcities.close();
	//setting up vectors, vars and output
	vector<vector<int>> newpop;
	vector<double> distance;
	vector<double> avdistance;
	vector<int> x,y;
	vector<vector<int>> cross;
	double sum = 0.;
	ofstream out;
	out.open(filename);
	int wd = 18;
	out << setw(wd) << "ntotiter" << endl;
	out << setw(wd) << niter << endl;
	out << setw(wd) << "Niter" << setw(wd) << "Best_route" << setw(wd) << "Elite_average" << endl;

	double scale = 1./2.;

	double p_swap = 0.05;
	double p_shift = 0.05;
	double p_swapcont = 0.05;
	double p_inversion = 0.05;
	double p_crossover = 0.8;
	
	// sorting the population
	FittingSort( pop, map );
	for (int i=0; i<niter; i++){
		sum = 0;
		cout << "iter " << i << "/" << niter << endl;
		pop.assign(pop.begin(),pop.begin()+nroutes); //eliminating individuals in excess
		for( int i=0; i<pop.size()/2; i++)
			sum += IndFitting(pop[i], map);
		out << setw(wd) << i <<  setw(wd) << IndFitting(pop[0], map) << setw(wd) << 2 * sum / pop.size() << endl;	//printing results
		//creating new population
		while (newpop.size() < pop.size()){
			// Swap
			if ( rnd.Rannyu() < p_swap){
				newpop.push_back(  Mutation_swap( Selection( pop, scale, rnd ), rnd ) );
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Shift
			if ( rnd.Rannyu() < p_shift){
				newpop.push_back(  Mutation_shift( Selection( pop, scale, rnd ), rnd ) );
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Swapcont
			if ( rnd.Rannyu() < p_swapcont){
				newpop.push_back(  Mutation_swapcont( Selection( pop, scale, rnd ), rnd ) );
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Inversion
			if ( rnd.Rannyu() < p_inversion){
				newpop.push_back(  Mutation_inversion( Selection( pop, scale, rnd ), rnd ) );
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Crossover
			if ( rnd.Rannyu() < p_crossover){
				for( int i=0; i<4; i++){
					x = Selection( pop, scale, rnd );
					y = Selection( pop, scale, rnd );
					cross = Crossover( x, y, rnd );
					newpop.push_back(cross[0]);
					//Check_Individual( newpop[newpop.size()-1] );
					newpop.push_back(cross[1]);
					//Check_Individual( newpop[newpop.size()-1] );
				}
			}
		}
		//Check_Population( newpop );
		pop = newpop;
		newpop.clear();
		FittingSort( pop, map );	//sortin the new population
		newpop.clear();
	}
	
	//printing results
	for( auto x : pop[0] )
		out << setw(wd/3) << x;
	out << endl;
	
	out.close();	
	return 0;
}
