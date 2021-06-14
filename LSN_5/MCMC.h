#ifndef __MCMC_h__
#define __MCMC_h__

#include <memory>
#include <vector>
#include <cmath>
#include "random.h"

class Position {

  public:
  // constructor
  Position();
  Position(double x, double y, double z); 
  // destroyer
  ~Position();
  
  // methods
  std::vector<double> getPos() const { return m_pos; };
  double getX() const { return m_pos[0];};       // cartesian coordinates
  double getY() const { return m_pos[1];};
  double getZ() const { return m_pos[2];};
  
	void SetCoord();																// setting cartesian coordinates					
  void SetCoord( double x, double y, double z );
  void setX( double x) { m_pos[0] = x;};
  void setY( double y) { m_pos[1] = y;};  
  void setZ( double z) { m_pos[2] = z;};
  
  double getR() const;       // spherical coordinates
  double getPhi() const;
  double getTheta() const;
  double getRho() const;     // radius of cylindrical coordinates

  protected:
  std::vector<double> m_pos;	// actual position

};

//virtual class for wavefunctions
class Basic_pdf{
  public:
    virtual ~Basic_pdf() = 0;
    virtual double Eval( Position x ) = 0;
};

//ground state wavefunction
class psi100 : public Basic_pdf{
  public:
    ~psi100(){};
    double Eval( Position x ){return exp(-2*x.getR()) / M_PI;};
};

//210 wavefunction
class psi210 : public Basic_pdf{
  public:
    ~psi210(){};
    double Eval( Position x ){return pow(x.getZ(),2)*exp(-x.getR())/(32.*M_PI);};
};


class MCMC : public Position {
  public:
    MCMC(std::shared_ptr<Basic_pdf> f, std::string name) : Position(){ 
      m_f = f;
      m_count = 0;
      m_name = name;
    };
    MCMC(std::shared_ptr<Basic_pdf> f, double x, double y, double z, double step, std::string name) : Position(x,y,z){
      m_f = f;
      m_step = step;
      m_count = 0;
      m_name = name;
    };
    void UnifStep( Random &rnd );   //MCMC uniform step
    void GaussStep( Random &rnd );  //MCMC gaussian step
    int m_count;
    std::string m_name;

  protected:
    std::shared_ptr<Basic_pdf> m_f;
    double m_step;
};

#endif // __MCMC_h__
