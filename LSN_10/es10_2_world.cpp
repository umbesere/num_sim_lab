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

template <class T>
T *resize(vector<vector<T>> mat)
{
	T *arr = new T[mat.size() * mat[0].size()];
	for (int i = 0; i < mat.size(); i++)
		for (int j = 0; j < mat[i].size(); j++)
			arr[i * mat[0].size() + j] = mat[i][j];
	return arr;
};

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

	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;
	MPI_Request req;

	vector<vector<double>> map;
	int ncities = 32;
	int nroutes = 125;
	int niter = 300;
	int wd = 18;
	Random rnd("seed" + to_string(rank) + ".in", "seed" + to_string(rank) + ".out");
	double *map1D;
	vector<vector<int>> world;
	int *world1D;
	if (rank == 0)
	{
		world1D = new int[ncities*nroutes*size];
		map = Square(ncities, rnd);
		map1D = resize<double>(map);
		ofstream outcities;
		outcities.open("map.txt");
		for (int i = 0; i < ncities; i++)
			outcities << setw(wd) << map[i][0] << setw(wd) << map[i][1] << endl;
		outcities.close();
	}
	else
	{
		map1D = new double[2*ncities];
		for (int i = 0; i < 2 * ncities; i++)
			map1D[i] = 0;
	}
	MPI_Bcast(map1D, 2 * ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (rank != 0)
		map = resize<double>(map1D, 2 * ncities, 2);

	string filename = "es10_2rank" + to_string(rank) + ".txt";
	string worldname = "world.txt";
	vector<vector<int>> pop = StartingPop(nroutes, ncities, rnd);
	int *pop1D = resize<int>(pop);
	MPI_Gather(pop1D, nroutes * ncities, MPI_INTEGER, world1D, nroutes * ncities, MPI_INTEGER, 0, MPI_COMM_WORLD);
	delete[] pop1D;
	
	//if (rank == 1 || rank == 2 || rank == 3){double a = rnd.Rannyu();while (a > 0.0000000001) a = rnd.Rannyu();}
	
	
	if (rank == 0)
	{
		world = resize<int>(world1D, nroutes * ncities * size, ncities);
		FittingSort(world, map);
		world.assign(world.begin(), world.begin() + nroutes);
	}
	//Check_Population( pop );
	vector<vector<int>> newpop;
	vector<int> x, y;
	vector<vector<int>> cross;
	double sum = 0.;
	ofstream out;
	out.open(filename);
	out << setw(wd) << "ntotiter" << endl;
	out << setw(wd) << niter << endl;
	out << setw(wd) << "Niter" << setw(wd) << "Best_route" << setw(wd) << "Elite_average" << endl;
	ofstream outworld;
	if (rank == 0)
	{
		outworld.open(worldname);
		outworld << setw(wd) << "ntotiter" << endl;
		outworld << setw(wd) << niter << endl;
		outworld << setw(wd) << "Niter" << setw(wd) << "Best_route" << setw(wd) << "Elite_average" << endl;
		outworld.close();
	}

	int ip, shift, nmigr;
	double fractomigr = 1./25.;
	vector<vector<int>> migrants;
	int *migrants1D_send;
	int *migrants1D_rec;

	double scale = 1. / 2.;

	double p_swap = 0.05;
	double p_shift = 0.05;
	double p_swapcont = 0.05;
	double p_inversion = 0.05;
	double p_crossover = 0.8;

	//FittingSort(pop, map);
	//cout << "check" << rank << endl;
	for (int i = 0; i < niter; i++)
	{
		FittingSort(pop, map);
		pop.assign(pop.begin(),pop.begin()+nroutes);
		//cout << "startfor" << rank << " " << i << endl;
		sum = 0;
		for (int j = 0; j < pop.size() / 2; j++)
			sum += IndFitting(pop[j], map);
		out << setw(wd) << i << setw(wd) << IndFitting(pop[0], map) << setw(wd) << 2 * sum / pop.size() << endl;
		sum = 0;
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
		FittingSort(newpop, map);
		//pop.assign(newpop.begin(), newpop.begin() + nroutes);
		if ( i % 10 == 0 )
		{
			//shift = 2*rnd.RannyuDiscr()-1;
			//if (rank == 0) shift = rnd.RannyuDiscr(1,4);
			//MPI_Bcast(&shift, 1, MPI_INTEGER, 0, MPI_COMM_WORLD);
			shift = 1;
			nmigr = nroutes*fractomigr;
			migrants1D_rec = new int[nmigr*ncities];
			for (int k = 0; k < nmigr; k++)
			{
				migrants.push_back(newpop[k]);
				newpop.erase(newpop.begin()+k);
			}
			migrants1D_send = resize<int>(migrants);
			//cout << "sendcheck" << rank << endl;
			for (int k = 0; k < size; k++)
			{
				if (rank == k)
					MPI_Isend(migrants1D_send, nmigr * ncities, MPI_INTEGER, (k + shift) % size, k, MPI_COMM_WORLD, &req);
				if (rank == (k + shift) % size)
					MPI_Recv(migrants1D_rec, nmigr * ncities, MPI_INTEGER, k, k, MPI_COMM_WORLD, &stat);
			}
			//cout << "postcheck" << rank << endl;			
			//if (rank == 1 ){double a = rnd.Rannyu();while (a > 0.0000000001) a = rnd.Rannyu();}
			migrants.clear();
			migrants = resize<int>(migrants1D_rec, nmigr * ncities, ncities);
			newpop.insert(newpop.begin(), migrants.begin(), migrants.end());
			//FittingSort(newpop, map);
			//newpop.assign(newpop.begin(), newpop.begin() + nroutes);
			delete[] migrants1D_rec;
			//cout << "check" << rank << endl;
			//cout << "doublecheck" << rank << endl;
			delete[] migrants1D_send;
			migrants.clear();
			//cout << "migrcheck" << rank << endl;
		}
		pop = newpop;
		newpop.clear();
		pop1D = resize<int>(pop);
		//cout << "gathercheck" << rank << endl;
		MPI_Gather(pop1D, nroutes * ncities, MPI_INTEGER, world1D, nroutes * ncities, MPI_INTEGER, 0, MPI_COMM_WORLD);
		delete[] pop1D;
		if (rank == 0)
		{
			world = resize<int>(world1D, nroutes * ncities * size, ncities);
//			delete[] world1D;
			FittingSort(world, map);
//			world.assign(world.begin(), world.begin() + nroutes);
			cout << "iter " << i << "/" << niter << endl;
			for (int j = 0; j < pop.size() / 2; j++)
				sum += IndFitting(world[j], map);
			outworld.open(worldname, ios::app);
			outworld << setw(wd) << i << setw(wd) << IndFitting(world[0], map) << setw(wd) << 2 * sum / world.size() << endl;
			outworld.close();
		}
		//cout << "endfor" << rank << endl;
	}

	for (auto x : pop[0])
		out << setw(wd / 3) << x;
	out << endl;

	if (rank == 0)
	{
		outworld.open(worldname, ios::app);
		for (auto x : world[0])
			outworld << setw(wd / 3) << x;
		outworld << endl;
		outworld.close();
	}

	outworld.clear();
	out.close();
	MPI_Finalize();
	return 0;
}
