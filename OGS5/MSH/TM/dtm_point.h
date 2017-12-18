/*

			point.h
			>Point class
			>>base class of Node class

			last update : 2003.12.03		by t.manabe

*/


#ifndef DTM_POINT_H
#define DTM_POINT_H 0


#include<math.h>
#include<iostream>

#include"dtm_error.h"

using namespace std;

namespace dtm{


	class Point		//point
	{
	private:
		double cd_x;
		double cd_y;
		double cd_z;

	public:
		/* constructor and destructor */

		Point() { }
		Point(const Point &pt): cd_x(pt.cd_x) , cd_y (pt.cd_y) , cd_z(pt.cd_z) { }
		Point(double x,double y,double z): cd_x(x) , cd_y (y) , cd_z(z) { }
		~Point(){ }


		/* original methods */

		virtual void setCoordinates(double x,double y,double z)
		{ cd_x = x; cd_y = y; cd_z = z;}
		virtual void setX(double x){ cd_x = x; }
		virtual void setY(double y){ cd_y = y; }
		virtual void setZ(double z){ cd_z = z; }
		virtual double getX() const { return cd_x; }
		virtual double getY() const { return cd_y; }
		virtual double getZ() const { return cd_z; }
		virtual void moveX(double x){ cd_x += x; }
		virtual void moveY(double y){ cd_y += y; }
		virtual void moveZ(double z){ cd_z += z; }
		virtual void move(Point mp){ cd_x += mp.cd_x; cd_y += mp.cd_y; cd_z += mp.cd_z; }
		virtual void reduce(double rate){ cd_x /= rate; cd_y /= rate; cd_z /= rate;}	//point reduce
		virtual void expand(double rate){ cd_x *= rate; cd_y *= rate; cd_z *= rate;}	//point expand
		virtual void reverse(){	cd_x = -cd_x;	cd_y = -cd_y; cd_z = -cd_z;}//change plus to minus(minus to plus)
		virtual double distance(Point &pt);	//get distance
		Point& operator=(const Point &pt){
			cd_x = pt.cd_x; cd_y = pt.cd_y;	cd_z = pt.cd_z;
			return *this;
		}
	
		//for debug
		virtual void view();
	
	};

}














#endif


//////////////////////////////////////////////////////////////////EOF