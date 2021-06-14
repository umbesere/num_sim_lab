#include <stdlib.h> // srand, rand: to generate random number
#include <iostream> // cin, cout: Standard Input/Output Streams Library
#include <fstream>	// Stream class to both read and write from/to files.
#include <cmath>	// rint, pow
#include <string>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main()
{
	Input(); //Inizialization
	int nconf = 1;
	for (int istep = 1; istep <= nstep; ++istep)
	{
		Move(); //Move particles with Verlet algorithm
		if (istep % iprint == 0)
			cout << "Number of time-steps: " << istep << endl;
		if (istep % 10 == 0)
		{
			Measure(); //Properties measurement
			//ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"!
			nconf += 1;
		}
	}
	ConfFinal(); //Write final configuration to restart
	int nblock = 50;
	Print(nconf - 1, nblock, "etot");
	Print(nconf - 1, nblock, "epot");
	Print(nconf - 1, nblock, "ekin");
	Print(nconf - 1, nblock, "temp");
	
	delete[] stima_gofr;
	return 0;
}

void Input(void)
{ //Prepare all stuff for the simulation
	ifstream ReadInput, ReadConf, ReadOldConf;
	double ep, ek, pr, et, vir;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl
		 << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl
		 << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	seed = 1;	 //Set seed for random numbers
	srand(seed); //Initialize random number generator

	ReadInput.open("input.dat"); //Read input

	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart / rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol, 1.0 / 3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> iprint;
	ReadInput >> restart;		// reading the restart var
	ReadInput >> eq;			// reading the equilibration var

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl
		 << endl;
	ReadInput.close();

	//Prepare array for measurements
	iv = 0;		 //Potential energy
	ik = 1;		 //Kinetic energy
	ie = 2;		 //Total energy
	it = 3;		 //Temperature
	n_props = 4; //Number of observables

	//measurement of g(r)
  	igofr = 2;
  	nbins = 100;
	stima_gofr = new double[nbins];
  	n_props = n_props + nbins;
  	bin_size = (box/2.0)/(double)nbins;
	glob_av = new double[n_props];
	glob_av2 = new double[n_props];

	//Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl
		 << endl;
	ReadConf.open("config.0");
	for (int i = 0; i < npart; ++i)
	{
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();

	double sumv2 = 0.0, fs;
	if (!restart)
	{
		//Prepare initial velocities
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl
			 << endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i = 0; i < npart; ++i)
		{
			vx[i] = rand() / double(RAND_MAX) - 0.5;
			vy[i] = rand() / double(RAND_MAX) - 0.5;
			vz[i] = rand() / double(RAND_MAX) - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim = 0; idim < 3; ++idim)
			sumv[idim] /= (double)npart;
		for (int i = 0; i < npart; ++i)
		{
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
		}
		sumv2 /= (double)npart;
		fs = sqrt(3 * temp / sumv2); // fs = velocity scale factor
		for (int i = 0; i < npart; ++i)
		{
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		}
	}
	else
	{	// restart using an old config

		//Read old configuration
		cout << "Read old configuration from file old.0 " << endl;
		ReadOldConf.open("old.0");
		for (int i = 0; i < npart; ++i)
		{
			ReadOldConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}
		ReadOldConf.close();

		if (eq)	// if one wants to equilibrate it is useful to rescale velocities
		{
			cout << "Rescaling velocities" << endl;
			Move();

			double sumv2 = 0.;
			for (int i = 0; i < npart; ++i)
			{
				vx[i] = Pbc(x[i] - xold[i]) / (delta);
				vy[i] = Pbc(y[i] - yold[i]) / (delta);
				vz[i] = Pbc(z[i] - zold[i]) / (delta);

				sumv2 += vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];
			}
			sumv2 /= (double)npart;

			double t = 0.;
			t = 0.5 * sumv2;
			t = (2.0 / 3.0) * t;
			fs = sqrt(temp / t); // fs = velocity scale factor

			for (int i = 0; i < npart; ++i)
			{
				vx[i] *= fs;
				vy[i] *= fs;
				vz[i] *= fs;

				xold[i] = Pbc(x[i] - vx[i] * delta);
				yold[i] = Pbc(y[i] - vy[i] * delta);
				zold[i] = Pbc(z[i] - vz[i] * delta);
			}
		}
		return;
	}
}

void Move(void)
{ //Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for (int i = 0; i < npart; ++i)
	{ //Force acting on particle i
		fx[i] = Force(i, 0);
		fy[i] = Force(i, 1);
		fz[i] = Force(i, 2);
	}

	for (int i = 0; i < npart; ++i)
	{ //Verlet integration scheme

		xnew = Pbc(2.0 * x[i] - xold[i] + fx[i] * pow(delta, 2));
		ynew = Pbc(2.0 * y[i] - yold[i] + fy[i] * pow(delta, 2));
		znew = Pbc(2.0 * z[i] - zold[i] + fz[i] * pow(delta, 2));

		vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
		vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
		vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}

double Force(int ip, int idir)
{ //Compute forces as -Grad_ip V(r)
	double f = 0.0;
	double dvec[3], dr;

	for (int i = 0; i < npart; ++i)
	{
		if (i != ip)
		{
			dvec[0] = Pbc(x[ip] - x[i]); // distance ip-i in pbc
			dvec[1] = Pbc(y[ip] - y[i]);
			dvec[2] = Pbc(z[ip] - z[i]);

			dr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
			dr = sqrt(dr);

			if (dr < rcut)
			{
				f += dvec[idir] * (48.0 / pow(dr, 14) - 24.0 / pow(dr, 8)); // -Grad_ip V(r)
			}
		}
	}

	return f;
}

void Measure()
{ //Properties measurement
	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	ofstream Epot, Ekin, Etot, Temp, Gofr;

	Epot.open("output_epot.dat", ios::app);
	Ekin.open("output_ekin.dat", ios::app);
	Temp.open("output_temp.dat", ios::app);
	Etot.open("output_etot.dat", ios::app);
	Gofr.open("output_gofr.dat", ios::app);

	int wd = 18;
	for (int i = 0; i < nbins; ++i){
		Gofr << setw(wd) << (i+1./2.)*bin_size;
	}
	Gofr << endl;

	v = 0.0; //reset observables
	t = 0.0;

	//cycle over pairs of particles
	for (int i = 0; i < npart - 1; ++i)
	{
		for (int j = i + 1; j < npart; ++j)
		{

			dx = Pbc(xold[i] - xold[j]); // here I use old configurations [old = r(t)]
			dy = Pbc(yold[i] - yold[j]); // to be compatible with EKin which uses v(t)
			dz = Pbc(zold[i] - zold[j]); // => EPot should be computed with r(t)

			dr = dx * dx + dy * dy + dz * dz;
			dr = sqrt(dr);

			//update of the histogram of g(r)
    		if (dr <= box/2)
    		{
      			for(int k=0; k<nbins; k++){
        			if ( k == int(dr/bin_size) )
          				stima_gofr[k] += 2;
				}
    		}

			if (dr < rcut)
			{
				vij = 4.0 / pow(dr, 12) - 4.0 / pow(dr, 6);

				//Potential energy
				v += vij;
			}
		}
	}

	//Kinetic energy
	for (int i = 0; i < npart; ++i)
		t += 0.5 * (vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i]);

	for (int i = 0; i < nbins; ++i){
		deltaV = 4./3.* pi * ( pow((i+1)*bin_size,3) - pow(i*bin_size,3) );
    	stima_gofr[i] /= rho * m_part * deltaV;
		Gofr << setw(wd) << (i+1./2.)*bin_size << setw(wd) << stima_gofr[i];
	}
	Gofr << endl;

	stima_pot = v / (double)npart;				  //Potential energy per particle
	stima_kin = t / (double)npart;				  //Kinetic energy per particle
	stima_temp = (2.0 / 3.0) * t / (double)npart; //Temperature
	stima_etot = (t + v) / (double)npart;		  //Total energy per particle

	Epot << stima_pot << endl;
	Ekin << stima_kin << endl;
	Temp << stima_temp << endl;
	Etot << stima_etot << endl;

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();
	Gofr.close();

	return;
}

void ConfFinal(void)
{ //Write final configuration
	ofstream WriteConf;
	ofstream WriteOldConf;

	cout << "Print final configuration to file config.final " << endl
		 << endl;
	WriteConf.open("config.final");
	cout << "Print pre-final configuration to file old.final " << endl
		 << endl;
	WriteOldConf.open("old.final");

	for (int i = 0; i < npart; ++i)
	{
		WriteOldConf << xold[i] / box << "   " << yold[i] / box << "   " << zold[i] / box << endl;
		WriteConf << x[i] / box << "   " << y[i] / box << "   " << z[i] / box << endl;
	}
	WriteConf.close();
	WriteOldConf.close();

	return;
}

void ConfXYZ(int nconf)
{ //Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i = 0; i < npart; ++i)
	{
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " << Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double Pbc(double r)
{ //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r / box);
}

// printing progressive means
void Print(int nthrows, int nblock, string input)
{
	ifstream in;
	ofstream out;
	int L = int(nthrows / nblock);
	double mea;
	double sum = 0., Sum = 0., Su2 = 0., err = 0.;
	in.open("output_" + input + ".dat");
	out.open("average_" + input + ".out");
	for (int i = 0; i < nblock; i++)
	{
		sum = 0;
		for (int k = 0; k < L; k++)
		{
			in >> mea;
			sum += mea;
		}
		sum /= L;
		Sum += sum;
		Su2 += sum * sum;
		err = error(Sum / (i + 1), Su2 / (i + 1), i);
		out << sum << " " << Sum / (i + 1) << " " << err << endl;
	}
	in.close();
	out.close();
}

// error function
double error(double ave, double av2, int n)
{
	if (n == 0)
		return 0;
	else
		return sqrt((av2 - pow(ave, 2)) / n);
}