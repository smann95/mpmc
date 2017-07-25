// Copyright 2015 Adam Hogan

#pragma once
class Vec
{
public:
	double components[3];

	
	~Vec();

	Vec();
	Vec ( double,double,double );
	void set ( double,double,double );
	double norm();
	double dot ( const Vec );
	Vec operator+ ( const Vec & );
	Vec operator* ( const double );
	
	Vec& operator= ( const Vec & );

};

Vec operator* ( const double x, const Vec right );