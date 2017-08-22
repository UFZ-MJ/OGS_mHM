/*

			dtm_calclate.cpp

			last update : 2003.11.01		by t.manabe

*/

#include "stdafx.h" /* MFC */

#include"dtm_calculate.h"

namespace dtm{

// distance

double getDistance(Point *p1,Point *p2)
{
	double dist;
	double x1,y1,z1;
	double x2,y2,z2;

	if(!p1 || !p2) {
		Errors err(0,
			"double getDistance(Point *p1,Point *p2)",
			"point is null");
		throw err;
	}

	x1 = p1->getX();
	y1 = p1->getY();
	z1 = p1->getZ();
	x2 = p2->getX();
	y2 = p2->getY();
	z2 = p2->getZ();

	dist = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
	dist = sqrt(dist);

	return dist;
}

// getDistanceSquare
//not use sprt

double getDistanceSquare(Point *p1,Point *p2)
{
	double dist;
	double x1,y1,z1;
	double x2,y2,z2;

	if(!p1 || !p2) {
		Errors err(0,
			"double getDistanceSquare(Point *p1,Point *p2)",
			"point is null");
		throw err;
	}

	x1 = p1->getX();
	y1 = p1->getY();
	z1 = p1->getZ();
	x2 = p2->getX();
	y2 = p2->getY();
	z2 = p2->getZ();

	dist = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);

	return dist;
}


//getCenter

Point getCenter(Point *p1,Point *p2)
{
	Point cp;

	if(!p1 || !p2) {
		Errors err(0,
			"Point getCenter(Point *p1,Point *p2)",
			"point is null");
		throw err;
	}

	cp.setX((p1->getX() + p2->getX())*0.5);
	cp.setY((p1->getY() + p2->getY())*0.5);
	cp.setZ((p1->getZ() + p2->getZ())*0.5);

	return cp;
}

//getTetraVolume

double getTetraVolume(Point *p0,Point *p1,Point *p2,Point *p3)
{
//	double volume;
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;
	double x3,y3,z3;
	double va,vb,vc;
	double wa,wb,wc;

	
	if((!p0) || (!p1) || (!p2) || (!p3)) {
		Errors err(0,
			"double getTetraVolume(Point *p0,Point *p1,Point *p2,Point *p3)",
			"point is null");
		throw err;
	}
	
    x0 = p0->getX();
    y0 = p0->getY();
    z0 = p0->getZ();
    x1 = p1->getX();
    y1 = p1->getY();
    z1 = p1->getZ();
    x2 = p2->getX();
    y2 = p2->getY();
    z2 = p2->getZ();
    x3 = p3->getX();
    y3 = p3->getY();
    z3 = p3->getZ();
    va = (x1-x0)*(y2-y0)*(z3-z0);
    vb = (y1-y0)*(z2-z0)*(x3-x0);
    vc = (z1-z0)*(x2-x0)*(y3-y0);
    wa = (x1-x0)*(z2-z0)*(y3-y0);
    wb = (y1-y0)*(x2-x0)*(z3-z0);
    wc = (z1-z0)*(y2-y0)*(x3-x0);

    //volume = va+vb+vc-wa-wb-wc;
	//return volume;

	return (va+vb+vc-wa-wb-wc);
}



// getTetraSphereCenter

Point getTetraSphereCenter(Point *p0,Point *p1,Point *p2,Point *p3,
						   double determ)
{
	Point center;
	double x0,y0,z0;
	double x1,y1,z1;
	double x2,y2,z2;
	double x3,y3,z3;
	double xv,yv,zv;
	double p11,p12,p13;
	double p21,p22,p23;
	double p31,p32,p33;
	double xyza;
	double aa,bb,cc;


	if((!p0) || (!p1) || (!p2) || (!p3)) {
		Errors err(0,
			"Point getTetraSphereCenter(Point *p0,Point *p1,Point *p2,Point *p3,double determ)",
			"point is null");
		throw err;
	}

    x0 = p0->getX();
    y0 = p0->getY();
    z0 = p0->getZ();
    x1 = p1->getX();
    y1 = p1->getY();
    z1 = p1->getZ();
    x2 = p2->getX();
    y2 = p2->getY();
    z2 = p2->getZ();
    x3 = p3->getX();
    y3 = p3->getY();
    z3 = p3->getZ();

    p11 = (y2-y0)*(z3-z0) - (y3-y0)*(z2-z0);
    p12 = (x3-x0)*(z2-z0) - (x2-x0)*(z3-z0);
    p13 = (x2-x0)*(y3-y0) - (x3-x0)*(y2-y0);
    p21 = (y3-y0)*(z1-z0) - (y1-y0)*(z3-z0);
    p22 = (x1-x0)*(z3-z0) - (x3-x0)*(z1-z0);
    p23 = (x3-x0)*(y1-y0) - (x1-x0)*(y3-y0);
    p31 = (y1-y0)*(z2-z0) - (y2-y0)*(z1-z0);
    p32 = (x2-x0)*(z1-z0) - (x1-x0)*(z2-z0);
    p33 = (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0);

    xyza = x0*x0 + y0*y0 + z0*z0;
    aa = 0.5*(x1*x1 + y1*y1 + z1*z1 - xyza);
    bb = 0.5*(x2*x2 + y2*y2 + z2*z2 - xyza);
    cc = 0.5*(x3*x3 + y3*y3 + z3*z3 - xyza);

    xv = (p11*aa + p21*bb + p31*cc)/determ;
    yv = (p12*aa + p22*bb + p32*cc)/determ;
    zv = (p13*aa + p23*bb + p33*cc)/determ;

	center.setCoordinates(xv,yv,zv);

	return center;
}


// getTriangleArea

double getTriangleArea(Point *p1,Point *p2,Point *p3)
{
	double ar;
	double x1,y1,z1;
	double x2,y2,z2;
	double x3,y3,z3;
	double a,b,c,s;

	x1 = p1->getX();
	y1 = p1->getY();
	z1 = p1->getZ();
	x2 = p2->getX();
	y2 = p2->getY();
	z2 = p2->getZ();
	x3 = p3->getX();
	y3 = p3->getY();
	z3 = p3->getZ();

	a = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	b = sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3));
	c = sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+(z3-z1)*(z3-z1));

	s = (a+b+c)*0.5;
	ar = fabs( s*(s-a)*(s-b)*(s-c) );
	//ar =  s*(s-a)*(s-b)*(s-c);
	//if(ar < 0) return -1;

	ar = sqrt(ar);

	return ar;
}


// getTriangleCeter

Point getTriangleCenter(Triangle *tri)
{
	Point center;
	Point *p0,*p1,*p2;
	double cx,cy,cz;

	p0 = tri->getNode(0);
	p1 = tri->getNode(1);
	p2 = tri->getNode(2);
	cx = (p0->getX() + p1->getX() + p2->getX())/3.;
	cy = (p0->getY() + p1->getY() + p2->getY())/3.;
	cz = (p0->getZ() + p1->getZ() + p2->getZ())/3.;

	center.setCoordinates(cx,cy,cz);

	return center;
}


// getTriangleCenter

Point getTriangleCenter(Point *p0,Point *p1,Point *p2)
{
	Point center;
	double cx,cy,cz;

	cx = (p0->getX() + p1->getX() + p2->getX())/3.;
	cy = (p0->getY() + p1->getY() + p2->getY())/3.;
	cz = (p0->getZ() + p1->getZ() + p2->getZ())/3.;

	center.setCoordinates(cx,cy,cz);

	return center;
}


// getTetraCenter

Point getTetraCenter(Tetra *tet)
{
	Point center;
	Point *p0,*p1,*p2,*p3;
	double cx,cy,cz;

	p0 = tet->getNode(0);
	p1 = tet->getNode(1);
	p2 = tet->getNode(2);
	p3 = tet->getNode(3);
	cx = (p0->getX() + p1->getX() + p2->getX() + p3->getX())*0.25;
	cy = (p0->getY() + p1->getY() + p2->getY() + p3->getY())*0.25;
	cz = (p0->getZ() + p1->getZ() + p2->getZ() + p3->getZ())*0.25;

	center.setCoordinates(cx,cy,cz);

	return center;

}


// getTetraCenter

Point getTetraCenter(Point *p0,Point *p1,Point *p2,Point *p3)
{
	Point center;
	double cx,cy,cz;

	cx = (p0->getX() + p1->getX() + p2->getX() + p3->getX())*0.25;
	cy = (p0->getY() + p1->getY() + p2->getY() + p3->getY())*0.25;
	cz = (p0->getZ() + p1->getZ() + p2->getZ() + p3->getZ())*0.25;

	center.setCoordinates(cx,cy,cz);

	return center;
}


// quickSort

void quickSort(int f,int e,deque<int>&array,deque<int>&key)
{
	int i,k;
	int pv;
	int tmpa;
	int tmpk;
	int flg;


	if(f >= e) return;

	flg = 1;
	pv = key[f];
	for(i=f+1;i<=e;i++) {
		if(key[i] > pv) {
			pv = key[i];
			flg = 0;
			break;
		} else if(key[i] < pv) {
			flg = 0;
			break;
		}
	}
	if(flg)return;

	k = e;

	for(i=f;i<=k;i++) {
		if(key[i] >= pv) {
			tmpa = array[i];
			tmpk = key[i];
			array[i] = array[k];
			key[i] = key[k];
			array[k] = tmpa;
			key[k] = tmpk;
			k--;
			i--;
		}
	}

	dtm::quickSort(f,i-1,array,key);
	dtm::quickSort(i,e,array,key);

	return;
}


}

//////////////////////////////////////////////////////////////////EOF