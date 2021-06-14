#include<iostream>
#include<cmath>
#include<cstdlib>
#include"MCMC.h"

// Position

// constructors
Position :: Position() {
	m_pos.push_back(0.);
	m_pos.push_back(0.);
	m_pos.push_back(0.);
};

Position :: Position(double x, double y, double z) {
	m_pos.push_back(x);
	m_pos.push_back(y);
	m_pos.push_back(z);
};

// destructor
Position :: ~Position() {

};

// setting cartesian coordinates
void Position :: SetCoord(){
	m_pos[0] = 0.;
	m_pos[1] = 0.;
	m_pos[2] = 0.;
};

void Position :: SetCoord( double x, double y, double z ){
	m_pos[0] = x;
	m_pos[1] = y;
	m_pos[2] = z;
};

// getting spherical coordinates
double Position :: getR() const {
	return sqrt(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2]);
};

double Position :: getTheta() const {
	return acos(m_pos[3]/(sqrt(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2])));
};

double Position :: getPhi() const {
	return asin(m_pos[2]/sqrt(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1]));
};

// getting radius of cylindrical coordinates
double Position :: getRho() const {
	return sqrt( m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] );
};

Basic_pdf :: ~Basic_pdf(){};

//Metropolis uniform step for MCMC
void MCMC :: UnifStep( Random &rnd ){
	Position p(rnd.Rannyu(m_pos[0]-m_step,m_pos[0]+m_step), rnd.Rannyu(m_pos[1]-m_step,m_pos[1]+m_step), rnd.Rannyu(m_pos[2]-m_step,m_pos[2]+m_step));
	if ( rnd.Rannyu() < fmin( 1, m_f->Eval(p)/m_f->Eval(*this) ) ){
		m_pos = p.getPos();
		m_count++;
	}
};

//Metropolis gaussian step for MCMC
void MCMC :: GaussStep( Random &rnd ){
	Position p(rnd.Gauss(m_pos[0], m_step), rnd.Gauss(m_pos[1], m_step), rnd.Gauss(m_pos[2], m_step));
	if ( rnd.Rannyu() < fmin( 1, m_f->Eval(p)/m_f->Eval(*this) ) ){
		m_pos = p.getPos();
		m_count++;
	}
};

