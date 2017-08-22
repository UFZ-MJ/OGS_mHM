/*

			dtm_calculate.h
			>functions of calculation

			last update : 2003.11.18		by t.manabe

*/

#ifndef DTM_CALCULATE_H
#define DTM_CALCULATE_H	2

#include<iostream>
#include<cmath>

#include"dtm_error.h"
#include"dtm_point.h"
#include"dtm_triangle.h"
#include"dtm_tetra.h"


using namespace std;

namespace dtm{

	class Triangle;
	class Tetra;

	double getDistance(Point *p1,Point *p2);
	double getDistanceSquare(Point *p1,Point *p2);
	Point getCenter(Point *p1,Point *p2);

	void quickSort(int f,int e,deque<int>&array,deque<int>&key);

	double getTetraVolume(Point *p0,Point *p1,Point *p2,Point *p3);
	Point getTetraSphereCenter(Point *p0,Point *p1,Point *p2,Point *p3,
						   double determ);
	double getTriangleArea(Point *p1,Point *p2,Point *p3);

	Point getTriangleCenter(Triangle *tri);
	Point getTriangleCenter(Point *p0,Point *p1,Point *p2);
	Point getTetraCenter(Tetra *tet);
	Point getTetraCenter(Point *p0,Point *p1,Point *p2,Point *p3);

}
#endif


//////////////////////////////////////////////////////////////////EOF