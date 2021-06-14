#include<string>

//parameters, observables
const int m_props=1000;
int n_props;
int iv,ik,it,ie,igofr;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_gofr;
int nbins; double bin_size, deltaV;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props], block_av[m_props];
double err_pot, err_kin, err_etot, err_temp, err_gofr;

//configuration
const int m_part=108;
bool restart, eq; // implementing the possibility to restart and rescale velocites
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed, nblock;
double delta;

//pigreco
const double pi=3.1415927;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(int, int);
double Pbc(double);
double Error(double, double, int);  // error function
void Print(int, int, std::string);  // pirnting progressive means

