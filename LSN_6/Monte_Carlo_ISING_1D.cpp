#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{
  Input();                                 //Inizialization
  for (int iblk = 1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk); //Reset block averages
    for (int istep = 1; istep <= nstep; ++istep)
    {
      Move(metro);
      Measure();
      Accumulate(); //Update block averages
    }
    Averages(iblk); //Print results for current block
  }
  ConfFinal(); //Write final configuration
  PrintFinal();
  return 0;
}

void Input(void)
{
  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl
       << endl;
  cout << "Nearest neighbour interaction      " << endl
       << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl
       << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

  //Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed, p1, p2);
  input.close();

  //Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0 / temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl
       << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> restart; //if=1 restart from config.0 else create a new one
  cout << restart << endl;
  ReadInput >> nblk;

  ReadInput >> nstep;

  if (metro == 1)
    cout << "The program perform Metropolis moves" << endl;
  else
    cout << "The program perform Gibbs moves" << endl;
  if (restart == 1)
    cout << "The program restarts from config.0" << endl;
  else
    cout << "The program starts from a new random configuration" << endl;

  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl
       << endl;
  ReadInput.close();

  //Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

  //initial configuration
  ifstream config0;
  config0.open("config.0");
  if (restart == 0)
  {
    //to start from a new configuration
    for (int i = 0; i < nspin; ++i)
    {
      //s[i]=1;
      if (rnd.Rannyu() >= 0.5)
        s[i] = 1;
      else
        s[i] = -1;
      //cout << s[i];
    }
    //cout << endl;
  }
  else
  {
    //to restart from an old configuration
    for (int i = 0; i < nspin; ++i)
    {
      config0 >> s[i];
      //cout << s[i];
    }
    //cout << endl;
  }
  config0.close();

  //Evaluate energy etc. of the initial configuration
  Measure();

  //Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu] / (double)nspin << endl;
}

void Move(int metro)
{
  int o;
  double p, energy_old, energy_new;
  double energy_up, energy_down;

  for (int i = 0; i < nspin; ++i)
  {
    //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu() * nspin);

    if (metro == 1) //Metropolis
    {
      attempted++;
      energy_old = Boltzmann(s[o], o);
      energy_new = Boltzmann(-s[o], o);
      p = fmin(1, exp(-beta * (energy_new - energy_old)));
      if (rnd.Rannyu() < p)
      {
        s[o] = -s[o];
        accepted++;
      }
    }
    else //Gibbs sampling
    {
      attempted++;
      accepted++;
      energy_up = Boltzmann(1., o);
      energy_down = Boltzmann(-1.,o);
      p = 1. / (1 + exp(-beta * (energy_down-energy_up)));      
      if (rnd.Rannyu() < p)
        s[o] = 1;
      else
        s[o] = -1;
    }
  }
}

double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * (s[Pbc(ip - 1)] + s[Pbc(ip + 1)]) - h * sm;
  return ene;
}

void Measure()
{
  double u = 0.0, m = 0.0;
  double en;

  //cycle over spins
  for (int i = 0; i < nspin; ++i)
  {
    en = -J * s[i] * s[Pbc(i + 1)] - 0.5 * h * (s[i] + s[Pbc(i + 1)]);
    u += en;
    m += s[i];
  }
  walker[iu] = u;
  walker[ic] = u*u;
  walker[im] = m;
  walker[ix] = m*m;
}

void Reset(int iblk) //Reset block averages
{

  if (iblk == 1)
  {
    for (int i = 0; i < n_props; ++i)
    {
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}

void Accumulate(void) //Update block averages
{

  for (int i = 0; i < n_props; ++i)
  {
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}

void Averages(int iblk) //Print results for current block
{

  ofstream Ene, Heat, Mag, Chi;
  const int wd = 18;

  cout << "Block number " << iblk << endl;
  cout << "Acceptance rate " << accepted / attempted << endl
       << endl;

  //calculating and printing Energy, Heat capacity, Magnetization, Magnetic Susceptibility (progressive mean)
  Ene.open("output.ene.0", ios::app);
  stima_u = blk_av[iu] / blk_norm / (double)nspin; //Energy
  glob_av[iu] += stima_u;
  glob_av2[iu] += stima_u * stima_u;
  err_u = Error(glob_av[iu], glob_av2[iu], iblk);
  Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu] / (double)iblk << setw(wd) << err_u << endl;
  Ene.close();

  Heat.open("output.heat.0", ios::app);
  stima_c = beta * beta * (blk_av[ic] / blk_norm - pow(blk_av[iu] / blk_norm,2)) / (double)nspin; //Heat capacity
  glob_av[ic] += stima_c;
  glob_av2[ic] += stima_c * stima_c;
  err_c = Error(glob_av[ic], glob_av2[ic], iblk);
  Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << glob_av[ic] / (double)iblk << setw(wd) << err_c << endl;
  Heat.close();

  Mag.open("output.mag.0", ios::app);
  stima_m = blk_av[im] / blk_norm / (double)nspin; //Magnetization
  glob_av[im] += stima_m;
  glob_av2[im] += stima_m * stima_m;
  err_m = Error(glob_av[im], glob_av2[im], iblk);
  Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im] / (double)iblk << setw(wd) << err_m << endl;
  Mag.close();

  Chi.open("output.chi.0", ios::app);
  stima_x = beta * (blk_av[ix] / blk_norm / (double)nspin); //Magnetic susceptibility
  glob_av[ix] += stima_x;
  glob_av2[ix] += stima_x * stima_x;
  err_x = Error(glob_av[ix], glob_av2[ix], iblk);
  Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix] / (double)iblk << setw(wd) << err_x << endl;
  Chi.close();

  cout << "----------------------------" << endl
       << endl;
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl
       << endl;
  WriteConf.open("config.final");
  for (int i = 0; i < nspin; ++i)
  {
    WriteConf << s[i] << endl;
    cout << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i) //Algorithm for periodic boundary conditions
{
  if (i >= nspin)
    i = i - nspin;
  else if (i < 0)
    i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk)
{
  if (iblk == 1)
    return 0.0;
  else
    return sqrt((sum2 / (double)iblk - pow(sum / (double)iblk, 2)) / (double)(iblk - 1));
}

void PrintFinal(){
  //printing final results
  const int wd = 18;
  ofstream out;
  string outfile;
  if (metro==1)
    outfile = "metro";
  else
    outfile = "gibbs";
  out.open("tempseries_"+outfile+".txt",ios::app);
  out << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << setw(wd)
   << glob_av[im]/(double)nblk << setw(wd) << err_m << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << setw(wd) << endl;
  out.close();
}


