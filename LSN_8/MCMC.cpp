#include<iostream>
#include<cmath>
#include<cstdlib>
#include<vector>
#include"MCMC.h"

// Position

// constructors
Position :: Position() {
	m_pos.push_back(0.);
	m_pos.push_back(0.);
	m_pos.push_back(0.);
	m_inpos = m_pos;
};

Position :: Position(std::vector<double> x) {
	for( int i=0; i<x.size(); i++)
		m_pos.push_back(x[i]);
	m_inpos = m_pos;
};

Position :: Position(double x, double y, double z) {
	m_pos.push_back(x);
	m_pos.push_back(y);
	m_pos.push_back(z);
	m_inpos = m_pos;
};

// destructor
Position :: ~Position() {
};

// setting cartesian coordinates
void Position :: SetCoord(){
	for( int i=0; i<m_pos.size(); i++)
		m_pos[i] = 0.;
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

//constructors
MCMC :: MCMC(std::shared_ptr<Basic_pdf> f, std::string name) : Position(){ 
    m_f = f;
    m_count = 0;
    m_name = name;
};
MCMC :: MCMC(std::shared_ptr<Basic_pdf> f, double x, double y, double z, double step, std::string name) : Position(x,y,z){
    m_f = f;
    m_step = step;
    m_count = 0;
    m_name = name;

};
MCMC :: MCMC(std::shared_ptr<Basic_pdf> f, std::vector<double> x, double step, std::string name) : Position(x){
	m_f = f;
    m_step = step;
    m_count = 0;
    m_name = name;

};
MCMC :: MCMC(std::shared_ptr<Basic_pdf> f, std::vector<double> x, double step) : Position(x){
	m_f = f;
    m_step = step;
    m_count = 0;

};

//Metropolis uniform step
void MCMC :: UnifStep(){
	std::vector<double> x(m_pos.size());
	for (int i=0; i<x.size(); i++)
		x[i] = m_pos[i]+m_rnd.Rannyu(-m_step/2,m_step/2.);
	Position pos(x);
	if ( m_rnd.Rannyu() < fmin( 1, m_f->Eval(pos)/m_f->Eval(*this) ) ){
		m_pos = pos.getPos();
		m_count++;
	}
};