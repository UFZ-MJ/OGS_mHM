/*

			dtm_plane.h
			>Plane class

			last update : 2003.11.24		by t.manabe

*/


#ifndef DTM_PLANE_H
#define DTM_PLANE_H 1

#include<math.h>

#include"dtm_error.h"
#include"dtm_point.h"

namespace dtm{


	class Plane		//Plane class
	{
	private:
		double plane_vector[4];	//vector of plane direction

	public:
		/* constructor and destructor */

		Plane(){ }
		Plane(Point *p0,Point *p1,Point *p2);


		/* original methods */

		void setVector(Point *p0,Point *p1,Point *p2);
		double getVector(int index){
			if(index > 3 || index < 0) {
				Errors err(0,
					"double getVector(int index)",
					"illegal index");
				throw err;
			}
			return plane_vector[index];
		}
		double getAnswer(Point *pt){
			double ans;
			ans = pt->getX()*plane_vector[0] + pt->getY()*plane_vector[1]
				+ pt->getZ()*plane_vector[2] + plane_vector[3];
			return ans;
		}
	};

}
#endif


//////////////////////////////////////////////////////////////////EOF