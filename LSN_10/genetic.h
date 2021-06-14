#ifndef __GeneticAlg__
#define __GeneticAlg__

#include <cstdlib>
#include <vector>
#include <memory>
#include <string>
#include "random.h"

using namespace std;

//to create the maps
vector<vector<double>> Square( int ncities, Random &rnd );
vector<vector<double>> Circle( int ncities, Random &rnd );

//sorting functions
double IndFitting( vector<int> x, vector<vector<double>> map );
void FittingSort( vector<vector<int>> &pop, vector<vector<double>> map );

//to generate the starting population/individual
vector<vector<int>> StartingPop( int pop_size, int ncities, Random &rnd );
vector<int> StartingInd( int ncities, Random &rnd );

//checking functions
void Check_Individual( vector<int> x );
void Check_Population( vector<vector<int>> x );

//selection, mutation and crossover
vector<int> Selection( vector<vector<int>> pop, double scale, Random &rnd );
vector<int> Mutation_swap( vector<int> x, Random &rnd);
vector<int> Mutation_shift( vector<int> x, Random &rnd);
vector<int> Mutation_swapcont( vector<int> x, Random &rnd );
vector<int> Mutation_inversion( vector<int> x, Random &rnd );
vector<vector<int>> Crossover( vector<int> x, vector<int> y, Random &rnd );


#endif //__GeneticAlg__