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
  Position(std::vector<double> x);
  Position(double x, double y, double z); 
  // destroyer
  ~Position();
  
  // methods
  std::vector<double> getPos() const { return m_pos; };
  std::vector<double> getInPos() const { return m_inpos; };
  double getX() const { return m_pos[0];};       // cartesian coordinates
  double getY() const { return m_pos[1];};
  double getZ() const { return m_pos[2];};
  double getCoord(int ip) const{ return m_pos[ip]; };


	void SetCoord();																// setting cartesian coordinates					
  void SetCoord( double x, double y, double z );
  void SetCoord( std::vector<double> x );
  void setCoord( int ip, double p) { m_pos[ip] = p;};
  void SetInCoord() { m_pos = m_inpos; }
  void setX( double x) { m_pos[0] = x;};
  void setY( double y) { m_pos[1] = y;};  
  void setZ( double z) { m_pos[2] = z;};
  
  
  double getR() const;       // spherical coordinates
  double getPhi() const;
  double getTheta() const;
  double getRho() const;     // radius of cylindrical coordinates

  protected:
  std::vector<double> m_pos;	// actual position
  std::vector<double> m_inpos;	// actual position

};

class Basic_pdf{
  public:
    virtual ~Basic_pdf() = 0;
    virtual double Eval( Position x ) = 0;
    virtual void SetPar( double par1, double par2 ) = 0;
  protected:
    std::vector<double> m_parms;

};

class wavefunct : public Basic_pdf{
  public:
    wavefunct(double mu, double sigma) : Basic_pdf() {m_parms.push_back(mu); m_parms.push_back(sigma);};
    ~wavefunct(){};
    double Eval( Position x ){
      double exp1 = exp(-pow(x.getCoord(0)-m_parms[0],2)/(2*pow(m_parms[1],2)));
      double exp2 = exp(-pow(x.getCoord(0)+m_parms[0],2)/(2*pow(m_parms[1],2)));
      return pow(exp1+exp2,2);
    };
    void SetPar(double mu, double sigma){ m_parms[0]=mu; m_parms[1]=sigma; };
};

class integrand : public Basic_pdf{
  public:
    integrand(double mu, double sigma) : Basic_pdf() {m_parms.push_back(mu); m_parms.push_back(sigma);};
    ~integrand(){};
    double Eval( Position x ){
      double exp1 = exp(-pow(x.getCoord(0)-m_parms[0],2)/(2*pow(m_parms[1],2)));
      double exp2 = exp(-pow(x.getCoord(0)+m_parms[0],2)/(2*pow(m_parms[1],2)));
      double num1 = (pow(x.getCoord(0)-m_parms[0],2)/pow(m_parms[1],2)-1)*exp1;
      double num2 = (pow(x.getCoord(0)+m_parms[0],2)/pow(m_parms[1],2)-1)*exp2;
      double den = exp1+exp2;
      return -(num1+num2)/den/(2*pow(m_parms[1],2)) + pow(x.getCoord(0),4) - 5./2.*pow(x.getCoord(0),2);
    };
    void SetPar(double mu, double sigma){ m_parms[0]=mu; m_parms[1]=sigma; };
};


class MCMC : public Position {
  public:
    MCMC(std::shared_ptr<Basic_pdf> f, std::string name);
    MCMC(std::shared_ptr<Basic_pdf> f, double x, double y, double z, double step, std::string name);
    MCMC(std::shared_ptr<Basic_pdf> f, std::vector<double> x, double step, std::string name);
    MCMC(std::shared_ptr<Basic_pdf> f, std::vector<double> x, double step);
    void SetStep( double step );
    void UnifStep();
    //void GaussStep();
    //void Stabilize(int num, double stepvar, double prec);
    int m_count;
    //double m_rate;
    double m_step;
    std::string m_name;

  protected:
    std::shared_ptr<Basic_pdf> m_f;
    Random m_rnd = Random("seed.out"); //NBNB RANDOM SPOSTATO DA BLOCKING A QUI (HA PIU' SENSO)
};

#endif // __MCMC_h__
