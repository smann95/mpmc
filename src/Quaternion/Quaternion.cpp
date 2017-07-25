// ©2013 Adam Hogan
// Space Research Group
// Department of Chemistry
// University of South Florida

// see http://www.cprogramming.com/tutorial/3d/quaternions.html
// and http://content.gpwiki.org/index.php/OpenGL:Tutorials:Using_Quaternions_to_represent_rotation

#include <cmath>

#include "Quaternion.h"
#include "constants.h"


Quaternion::Quaternion( double x, double y, double z, double w, int construction_mode )
{
	double angle     = w;
	double sinAngle  = 0,
	       magnitude = 0;

	switch( construction_mode ) {

		case XYZW:
			// Construct quaternion from components
			this->x = x;
			this->y = y;
			this->z = z;
			this->w = w;
			break;

		case AXIS_ANGLE_DEGREE:
			// Convert angle from degrees to radians, prior to construction of quaternion
			angle /= 57.2957795; // fall-through here is intentional...

		case AXIS_ANGLE_RADIAN:

			// Construct quaternion from an axis and an angle (in radians)
			// Normalizes the axis vector

			magnitude = sqrt( x*x + y*y + z*z );
			if( magnitude == 0.0 ) // edge case, if the axis to rotate around doesn't exist just return a quaternion with no rotation
			{
				this->x = 0;
				this->y = 0;
				this->z = 0;
				this->w = 1;
				return;
			}

			x = x/magnitude;
			y = y/magnitude;
			z = z/magnitude;

			sinAngle = sin(angle/2.0);
		
			this->x = x*sinAngle;
			this->y = y*sinAngle;
			this->z = z*sinAngle;
			this->w = cos(angle/2.0);
			break;
	
		default:
			throw invalid_quaternion_mode;
			break;
	}
}


Quaternion::~Quaternion()
{
}


// Normalize quaternion
void Quaternion::normalize()
{
	double magnitude = sqrt( x*x + y*y + z*z + w*w );
	x = x/magnitude;
	y = y/magnitude;
	z = z/magnitude;
	w = w/magnitude;
}

Quaternion& Quaternion::operator=(const Quaternion &rhs)
{
	if( &rhs != this ) {  // Avoids self assignment
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		w = rhs.w;
	}
	return *this; // allows chaining, i.e. x = y = z;
}


// QuaternionStore = Q1 * Q2
// Order matters!
Quaternion Quaternion::operator*( const Quaternion &rhs )
{
  double result_w = w*rhs.w - x*rhs.x - y*rhs.y - z*rhs.z;
  double result_x = w*rhs.x + x*rhs.w + y*rhs.z - z*rhs.y;
  double result_y = w*rhs.y - x*rhs.z + y*rhs.w + z*rhs.x;
  double result_z = w*rhs.z + x*rhs.y - y*rhs.x + z*rhs.w;
  
  return Quaternion( result_x, result_y, result_z, result_w, XYZW );
}

// A conjugate quaternion performs the opposite rotation
Quaternion Quaternion::conjugate()
{
	return Quaternion( -x, -y, -z, w, XYZW );
}
