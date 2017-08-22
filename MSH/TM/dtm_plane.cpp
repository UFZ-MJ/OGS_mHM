/*

			dtm_plane.cpp
			>Plane class

			last update : 2003.11.24

*/

#include "stdafx.h" /* MFC */


#include"dtm_plane.h"

namespace dtm{


// Plane
///03.11.24

Plane::Plane(Point *p0,Point *p1,Point *p2)
{
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;

	x0 = p0->getX();
	y0 = p0->getY();
    z0 = p0->getZ();
	x1 = p1->getX();
	y1 = p1->getY();
	z1 = p1->getZ();
	x2 = p2->getX();
	y2 = p2->getY();
	z2 = p2->getZ();
	plane_vector[0] = y0*z1 + y1*z2 + y2*z0 - y0*z2 - y1*z0 - y2*z1;
	plane_vector[1] = z0*x1 + z1*x2 + z2*x0 - z0*x2 - z1*x0 - z2*x1;
	plane_vector[2] = x0*y1 + x1*y2 + x2*y0 - x0*y2 - x1*y0 - x2*y1;

	/*//change basic vector
	double dist;
	dist = vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2];
	dist = sqrt(dist);
	vect[0] = vect[0]/dist;
	vect[1] = vect[1]/dist;
	vect[2] = vect[2]/dist;
	*//////////////////////
	plane_vector[3] = -plane_vector[0]*x0 - plane_vector[1]*y0 - plane_vector[2]*z0;

}


// setVector
///03.11.24

void Plane::setVector(Point *p0,Point *p1,Point *p2)
{
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;
	double a,b,c,d;

	x0 = p0->getX();
	y0 = p0->getY();
    z0 = p0->getZ();
	x1 = p1->getX();
	y1 = p1->getY();
	z1 = p1->getZ();
	x2 = p2->getX();
	y2 = p2->getY();
	z2 = p2->getZ();
	a = y0*z1 + y1*z2 + y2*z0 - y0*z2 - y1*z0 - y2*z1;
	b = z0*x1 + z1*x2 + z2*x0 - z0*x2 - z1*x0 - z2*x1;
	c = x0*y1 + x1*y2 + x2*y0 - x0*y2 - x1*y0 - x2*y1;
	/*//change basic vector
	double dist;
	dist = a*a+b*b+c*c;
	dist = sqrt(dist);
	a = a/dist;
	b = b/dist;
	c = c/dist;
	*//////////////////////
	d = -1.*(a*x0 + b*y0 + c*z0);
	plane_vector[0] = a;
	plane_vector[1] = b;
	plane_vector[2] = c;
	plane_vector[3] = d;


	return;
}


}

//////////////////////////////////////////////////////////////////EOF