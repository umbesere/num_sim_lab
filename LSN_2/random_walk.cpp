#include<iostream>
#include<cmath>
#include<cstdlib>
#include"random.h"
#include"random_walk.h"

// RandomWalk

// constructors
RandomWalk :: RandomWalk() {
	m_pos[0] = 0;
	m_pos[0] = 0;
	m_pos[0] = 0;
};

RandomWalk :: RandomWalk(double x, double y, double z) {
	m_pos[0] = x;
	m_pos[1] = y;
	m_pos[2] = z;
};

// destructor
RandomWalk :: ~RandomWalk() {

};

// setting cartesian coordinates
void RandomWalk :: SetCoord(){
	m_pos[0] = 0.;
	m_pos[1] = 0.;
	m_pos[2] = 0.;
};

void RandomWalk :: SetCoord( double x, double y, double z ){
	m_pos[0] = x;
	m_pos[1] = y;
	m_pos[2] = z;
};

// getting spherical coordinates
double RandomWalk :: getR() const {
	return sqrt(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2]);
};

double RandomWalk :: getTheta() const {
	return acos(m_pos[3]/(sqrt(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] + m_pos[2]*m_pos[2])));
};

double RandomWalk :: getPhi() const {
	return asin(m_pos[2]/sqrt(m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1]));
};

// getting radius of cylindrical coordinates
double RandomWalk :: getRho() const {
	return sqrt( m_pos[0]*m_pos[0] + m_pos[1]*m_pos[1] );
};

// discrete step
void RandomWalk :: DiscrStep( Random &rnd ){
	int step = rnd.RannyuDiscr()*2-1;		// generating random step
	int dir = rnd.RannyuDiscr(0,3);			// generating random direction
	m_pos[dir] += step;									// making the step
};

// continuous step
void RandomWalk :: ContStep( Random &rnd ){
	double phi, theta;
	rnd.SolidAngle(theta, phi);					// generating random solid angles
	m_pos[0] += sin(theta)*cos(phi);		// making the step
	m_pos[1] += sin(theta)*sin(phi);
	m_pos[2] += cos(theta);
};


