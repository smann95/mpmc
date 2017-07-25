// Copyright 2015 Adam Hogan

#include <math.h>

#include "Vec.h"


Vec::~Vec(){}

Vec::Vec () {
	for ( int i=0; i<3; i++ ) {
		components[i] = 0;
	}
}

Vec::Vec ( double x, double y, double z ) {
	components[0] = x;
	components[1] = y;
	components[2] = z;
}

void Vec::set ( double x, double y, double z ) {
	components[0] = x;
	components[1] = y;
	components[2] = z;
}

double Vec::dot ( const Vec other ) {
	double sum = 0.0;
	for ( int i=0; i<3; i++ ) {
		sum += components[i]*other.components[i];
	}
	return sum;
}

double Vec::norm() {
	return sqrt ( dot ( *this ) );
}

Vec Vec::operator+ ( const Vec &right ) {
	Vec result ( 0,0,0 );
	for ( int i=0; i<3; i++ ) {
		result.components[i] = components[i]+right.components[i];
	}
	return result;
}

Vec Vec::operator* ( const double x ) {
	Vec result ( 0,0,0 );
	for ( int i=0; i<3; i++ ) {
		result.components[i] = components[i]*x;
	}
	return result;
}

Vec operator* ( const double x, const Vec right ) {
	Vec result ( 0,0,0 );
	for ( int i=0; i<3; i++ ) {
		result.components[i] = right.components[i]*x;
	}
	return result;
}

Vec& Vec::operator= ( const Vec &right ) {
	if( &right != this ) {  // Avoids self assignment
		components[0] = right.components[0];
		components[1] = right.components[1];
		components[2] = right.components[2];
	}
	return *this; // allows chaining, i.e. x = y = z;
}
