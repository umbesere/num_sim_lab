#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "random.h"
#include "genetic.h"

using namespace std;

//Generates square map
vector<vector<double>> Square( int ncities, Random &rnd ){
    vector<vector<double>> v(ncities, vector<double>(2));
    for( int i=0; i<ncities; i++)
        for( int j=0; j<2; j++)
            v[i][j] = rnd.Rannyu(-1.,1.);
    return v;
};

//Generates circle map
vector<vector<double>> Circle( int ncities, Random &rnd ){
    vector<vector<double>> v(ncities, vector<double>(2));
    double angle;
    for( int i=0; i<ncities; i++){
        angle = rnd.Rannyu(0.,2*M_PI);
        v[i][0] = cos(angle);
        v[i][1] = sin(angle);
    }
    return v;
};

//generates the starting population
vector<vector<int>> StartingPop( int pop_size, int ncities, Random &rnd){
    vector<vector<int>> v(pop_size, vector<int>(ncities));
    for(int j=0; j<ncities; j++)
        v[0][j] = j;
    for(int i=1; i<pop_size; i++)
        v[i] = Mutation_swap(v[i-1],rnd);
    return v;
};

//checks an individual
void Check_Individual( vector<int> x ){
    if (x[0] != 0)
        cerr << "you are not starting from 0!!!" << endl;
    for(int k=0; k<x.size()-1; k++){
        for(int j=k+1; j<x.size(); j++){
            if (x[j] == x[k]){
                cerr << "Error in mutation/crossover!!!" << endl;
                cerr << j << " " << k << endl;
                exit(1);
            }
        }
    }
};

//checks the entire poppulation
void Check_Population( vector<vector<int>> x ){
    for(int i=0; i<x.size(); i++){
        if (x[i][0] != 0)
            cerr << "you are not starting from 0!!!" << endl;
        for(int k=0; k<x[i].size()-1; k++){
            for(int j=k+1; j<x[i].size(); j++){
                if (x[i][j] == x[i][k]){
                    cerr << "Error in mutation/crossover!!!" << endl;
                    cerr << i << " " << j << " " << k << endl;
                    exit(1);
                }
            }
        }
    }
}

//computes path length of an individual
double IndFitting(vector<int> route, vector<vector<double>> map){
    double sum = 0;
    for(int i=0; i<route.size()-1; i++)
        for(int j=0; j<2; j++)
            sum += pow(map[route[i]][j]-map[route[i+1]][j],2);
    sum += pow(map[route[0]][0]-map[route[route.size()-1]][0],2) + pow(map[route[0]][1]-map[route[route.size()-1]][1],2);
    return sum;
};

//rearranges the population
void FittingSort(vector<vector<int>> &pop, vector<vector<double>> map){
    vector<int> a;
    for(int i=0; i<pop.size()-1; i++){
        for(int j=i+1; j<pop.size(); j++){
            if ( IndFitting(pop[i], map) > IndFitting(pop[j], map) ){
                a = pop[i];
                pop[i] = pop[j];
                pop[j] = a;
            }
        }
    }
};

//selection operator
vector<int> Selection( vector<vector<int>> pop, double scale, Random &rnd ){
    int ip;
    while (ip >= pop.size())
        ip = int( rnd.Exponential( 1./(scale*pop.size()) ) );
    return pop[ip];
};

//Mutations
vector<int> Mutation_swap( vector<int> x, Random &rnd){
    int ip1 = rnd.RannyuDiscr(1, x.size());
    int ip2 = rnd.RannyuDiscr(1, x.size());
    double a = x[ip1];
    x[ip1] = x[ip2];
    x[ip2] = a;
    return x;
};

vector<int> Mutation_shift( vector<int> x, Random &rnd ){
    vector<int> y(x.begin()+1, x.end());
    int start = rnd.RannyuDiscr(0, y.size());
    int shift = rnd.RannyuDiscr(1, y.size());
    int numtoshift = rnd.RannyuDiscr(1, y.size());
    vector<int> a;
    for(int i=0; i<numtoshift; i++)
        a.push_back(y[(start+i)%y.size()]);
    for(int i=0; i<shift; i++)
        y[(start+i)%y.size()] = y[(start+i+numtoshift)%y.size()];
    for(int i=0; i<numtoshift; i++)
        y[(start+shift+i)%y.size()] = a[i];
    y.insert(y.begin(), 0);
    return y;
};

vector<int> Mutation_swapcont( vector<int> x, Random &rnd ){
    vector<int> y(x.begin()+1, x.end());
    int start1 = rnd.RannyuDiscr(0, y.size()-1);
    int num = rnd.RannyuDiscr(0, y.size()/2.);
    int start2;
    if (start1 < (start1+num)%y.size()){
        start2 = rnd.RannyuDiscr(start1+num+1, y.size()+start1)%y.size();
    } else {
        start2 = rnd.RannyuDiscr((start1+num+1)%y.size(),start1);
    }
    vector<int> a(num+1,0);
    for(int i=0; i<num+1; i++){
        a[i] = y[(start1+i)%y.size()];
        y[(start1+i)%y.size()] = y[(start2+i)%y.size()];
        y[(start2+i)%y.size()] = a[i];
    }
    y.insert(y.begin(), 0);
    return y;    
}

vector<int> Mutation_inversion( vector<int> x, Random &rnd ){
    vector<int> y(x.begin()+1, x.end());
    int start = rnd.RannyuDiscr(0, y.size());
    int num = rnd.RannyuDiscr(0, y.size());
    int a;
    for(int i=0; i<(num+1)/2; i++){
        a = y[(start+i)%y.size()];
        y[(start+i)%y.size()] = y[(start+num-i)%y.size()];
        y[(start+num-i)%y.size()] = a;
    }
    y.insert(y.begin(), 0);
    return y;
}

//Crossover operator
vector<vector<int>> Crossover( vector<int> x, vector<int> y, Random &rnd ){
    int cut = rnd.RannyuDiscr(1, x.size());
    vector<int> x1(x.begin(), x.end()-cut);
    vector<int> x2(x.end()-cut, x.end());
    vector<int> y1(y.begin(), y.end()-cut);
    vector<int> y2(y.end()-cut, y.end());
    for (int i=0; i<y.size(); i++)
        for(int j=0; j<cut; j++)
            if ( x2[j] == y[i] )
                x1.push_back(y[i]);
    for (int i=0; i<x.size(); i++)
        for(int j=0; j<cut; j++)
            if ( y2[j] == x[i] )
                y1.push_back(x[i]);
    vector<vector<int>> res = {x1,y1};
    return res;
};
