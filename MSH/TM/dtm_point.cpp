/*

			dtm_point.cpp
			>Point class

			last update : 2003.11.23		by t.manabe

*/


#include "stdafx.h" /* MFC */


#include"dtm_point.h"

namespace dtm{


// distance
///03.11.23

double Point::distance(Point &pt)
{
	double dist;

	dist = (pt.cd_x - cd_x)*(pt.cd_x - cd_x)+(pt.cd_y - cd_y)*(pt.cd_y - cd_y)+(pt.cd_z - cd_z)*(pt.cd_z - cd_z);
	dist = sqrt(dist);

	return dist;
}


// view

void Point::view()
{
	cout << " x = " << cd_x << " y = " << cd_y << " z = " << cd_z << "\n";

	return;
}



}

//////////////////////////////////////////////////////////////////EOF