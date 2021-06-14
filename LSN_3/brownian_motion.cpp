#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <memory>
#include "brownian_motion.h"
#include "random.h"

using namespace std;

// standard BM
BrownianMotion :: BrownianMotion(){
	m_pos = 0.;
	m_drift = 0.;
	m_vol = 1.;
}

// BM constructor
BrownianMotion :: BrownianMotion(double pos, double drift, double vol){
	m_inpos = pos;
	m_pos = pos;
	m_drift = drift;
	m_vol = vol;
}

// destructor
BrownianMotion :: ~BrownianMotion(){}

// GBM step
void BrownianMotion :: Gstep( shared_ptr<Random> rnd, double time_step ){
	m_pos *= exp( (m_drift - m_vol*m_vol/2)*time_step + m_vol*rnd->Gauss()*sqrt(time_step) );
}

// nstep of GBM
void BrownianMotion :: Gpath( shared_ptr<Random> rnd, double tot_time, int nstep ){
	for (int i = 0; i < nstep; i++)
		Gstep(rnd, tot_time/nstep);
};



