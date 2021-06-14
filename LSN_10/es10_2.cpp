#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <iomanip>
#include <mpi.h>
#include "random.h"
#include "genetic.h"

using namespace std;

//transforms matrix into array 1D
template <class T>
T *resize(vector<vector<T>> mat)
{
	T *arr = new T[mat.size() * mat[0].size()];
	for (int i = 0; i < mat.size(); i++)
		for (int j = 0; j < mat[i].size(); j++)
			arr[i * mat[0].size() + j] = mat[i][j];
	return arr;
};

//retransforms array 1D into matrix
template <class T>
vector<vector<T>> resize(T *v, int size, int rowsize)
{
	vector<vector<T>> mat(int(size / rowsize), vector<T>(rowsize, 0));
	for (int i = 0; i < mat.size(); i++)
		for (int j = 0; j < mat[i].size(); j++)
			mat[i][j] = v[i * rowsize + j];
	return mat;
};

int main(int argc, char *argv[])
{

	//PARALLEL GENETIC ALGORITHM//
	//setting up parallel vars
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;
	MPI_Request req;

	//setting up the GA
	vector<vector<double>> map;
	int ncities = 32;
	int nroutes = 500;
	int niter = 1000;
	int wd = 18;
	Random rnd("seed" + to_string(rank) + ".in", "seed" + to_string(rank) + ".out");
	string filename = "es10_2rank" + to_string(rank) + ".txt";
	//creating the starting population
	vector<vector<int>> pop = StartingPop(nroutes, ncities, rnd);
	double *map1D;	//1D version of the map
	if (rank == 0)	//uncomment to start from a new map 
	{
		//map = Square(ncities, rnd);
		ifstream in("map.txt", ios::in);
		for(int i=0;!in.eof();i++){
			map.push_back(vector<double>(2,0));
			in >> map[i][0] >> map[i][1];
		}
		in.close();
		map1D = resize<double>(map);
		/*
		ofstream outcities;
		outcities.open("map.txt");
		for (int i = 0; i < ncities; i++)
			outcities << setw(wd) << map[i][0] << setw(wd) << map[i][1] << endl;
		outcities.close();
		*/
	}
	else
	{
		map1D = new double[2*ncities];
		for (int i = 0; i < 2 * ncities; i++)
			map1D[i] = 0;
	}
	MPI_Bcast(map1D, 2 * ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);	//cast the map from root to other processes
	if (rank != 0)
		map = resize<double>(map1D, 2 * ncities, 2); //resizing the map
	//Check_Population( pop );
	
	//setting up vectors for GA, output and vars
	vector<vector<int>> newpop;
	vector<int> x, y;
	vector<vector<int>> cross;
	double sum = 0.;
	ofstream out;
	out.open(filename);
	out << setw(wd) << "ntotiter" << endl;
	out << setw(wd) << niter << endl;
	out << setw(wd) << "Niter" << setw(wd) << "Best_route" << setw(wd) << "Elite_average" << endl;
	//initializing migration vars
	int ip, shift, nmigr;
	double fractomigr = 1./25.;
	vector<vector<int>> migrants;
	int *migrants1D_send;
	int *migrants1D_rec;

	//probabilities
	double scale = 1. / 2.;
	double p_swap = 0.05;
	double p_shift = 0.05;
	double p_swapcont = 0.05;
	double p_inversion = 0.05;
	double p_crossover = 0.8;

	for (int i = 0; i < niter; i++)
	{
		FittingSort(pop, map);		//sorting the populations
		pop.assign(pop.begin(),pop.begin()+nroutes);	//eliminating individuals in excess
		sum = 0;
		if ( rank == 0) cout << "iter" << i << "/" << niter << endl;
		for (int j = 0; j < pop.size() / 2; j++)
			sum += IndFitting(pop[j], map);
		out << setw(wd) << i << setw(wd) << IndFitting(pop[0], map) << setw(wd) << 2 * sum / pop.size() << endl;	//output
		sum = 0;
		//creating the new population
		while (newpop.size() < pop.size())
		{
			// Swap
			if (rnd.Rannyu() < p_swap)
			{
				newpop.push_back(Mutation_swap(Selection(pop, scale, rnd), rnd));
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Shift
			if (rnd.Rannyu() < p_shift)
			{
				newpop.push_back(Mutation_shift(Selection(pop, scale, rnd), rnd));
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Swapcont
			if (rnd.Rannyu() < p_swapcont)
			{
				newpop.push_back(Mutation_swapcont(Selection(pop, scale, rnd), rnd));
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Inversion
			if (rnd.Rannyu() < p_inversion)
			{
				newpop.push_back(Mutation_inversion(Selection(pop, scale, rnd), rnd));
				//Check_Individual( newpop[newpop.size()-1] );
			}
			// Crossover
			if (rnd.Rannyu() < p_crossover)
			{
				for (int i = 0; i < 4; i++)
				{
					x = Selection(pop, scale, rnd);
					y = Selection(pop, scale, rnd);
					cross = Crossover(x, y, rnd);
					newpop.push_back(cross[0]);
					//Check_Individual( newpop[newpop.size()-1] );
					newpop.push_back(cross[1]);
					//Check_Individual( newpop[newpop.size()-1] );
				}
				x.clear();
				y.clear();
				cross.clear();
			}
		}
		//Check_Population( newpop );
		FittingSort(newpop, map);	//sorting the new populations
		if ( i % 10 == 0 )		//migration
		{
			//uncomment to have random exchange
			//if (rank == 0) shift = rnd.RannyuDiscr(1,4);
			//MPI_Bcast(&shift, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
			shift = 1;
			nmigr = nroutes*fractomigr;
			migrants1D_rec = new int[nmigr*ncities];	//reciving var
			for (int k = 0; k < nmigr; k++)	//selecting the migrants
			{
				migrants.push_back(newpop[k]);
				newpop.erase(newpop.begin()+k);
			}
			migrants1D_send = resize<int>(migrants);
			for (int k = 0; k < size; k++)		//migrating
			{
				if (rank == k)
					MPI_Send(migrants1D_send, nmigr * ncities, MPI_INTEGER, (k + shift) % size, k, MPI_COMM_WORLD);
				if (rank == (k + shift) % size)
					MPI_Recv(migrants1D_rec, nmigr * ncities, MPI_INTEGER, k, k, MPI_COMM_WORLD, &stat);
			}
			migrants.clear();
			migrants = resize<int>(migrants1D_rec, nmigr * ncities, ncities);
			newpop.insert(newpop.begin(), migrants.begin(), migrants.end());
			delete[] migrants1D_rec;
			delete[] migrants1D_send;
			migrants.clear();
		}
		pop = newpop;
		newpop.clear();
	}

	//printing the output
	for (auto x : pop[0])
		out << setw(wd / 3) << x;
	out << endl;

	out.close();
	MPI_Finalize();
	return 0;
}
