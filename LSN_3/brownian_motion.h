#ifndef __BrownianMotion__
#define __BrownianMotion__

#include"random.h"
#include<memory>

class BrownianMotion {		//BrownianMotion class

	public:
	
		// constructors
		BrownianMotion();											// stadard BM
		BrownianMotion( double pos, double drift, double vol);		// BM
		~BrownianMotion();											// destructor
		
		void SetPos(){ m_pos = m_inpos; };		// reset position
		
		// getting params
		double GetPos() const { return m_pos; };
		double GetDrift() const { return m_drift; };
		double GetVol() const { return m_vol; };		
		
		// simulations
		void Gstep( std::shared_ptr<Random> rnd, double time_step );				// GBM step
		void Gpath( std::shared_ptr<Random> rnd, double tot_time, int nstep );		// nstep of GBM
				
	private:

		double m_inpos;		// initial position
		double m_pos;		// actual position
		double m_drift;		// drift
		double m_vol;		// volatility
};

#endif // __BrownianMotion__
