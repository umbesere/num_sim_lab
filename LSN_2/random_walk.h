#ifndef __RandomWalk_h__
#define __RandomWalk_h__

class RandomWalk {

public:

  // constructur
  RandomWalk();
  RandomWalk(double x, double y, double z); 
  // destroyer
  ~RandomWalk();
  
  // methods
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
  
  void DiscrStep( Random &rnd ); // simulating steps
  void ContStep( Random &rnd );
  
private:

  double m_pos[3];	// actual position

};

#endif // __RandomWalk_h__
