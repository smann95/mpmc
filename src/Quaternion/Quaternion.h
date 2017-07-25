// ©2013 Adam Hogan
// Space Research Group
// Department of Chemistry
// University of South Florida

// see http://www.cprogramming.com/tutorial/3d/quaternions.html
// and http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation

#pragma once
class Quaternion
{
public:
	Quaternion( double x, double y, double z, double w, int mode );
	~Quaternion();

	static const int XYZW = 1;
	static const int AXIS_ANGLE_RADIAN = 2;
	static const int AXIS_ANGLE_DEGREE = 3;

	void normalize();
	Quaternion  operator*( const Quaternion &rhs );
	Quaternion& operator=( const Quaternion &rhs );
	Quaternion  conjugate();

	inline double X() { return x; }
	inline double Y() { return y; }
	inline double Z() { return z; }

private:
	double x, y, z, w;
  
};

