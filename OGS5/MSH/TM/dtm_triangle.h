/*

			dtm_triangleh
			>Triangle class
			>>based Figure Class

			last update : 2003.12.17		by t.manabe

*/


#ifndef DTM_TRIANGLE_H
#define DTM_TRIANGLE_H	2

#include<iostream>
#include<math.h>

#include"dtm_error.h"
#include"dtm_node.h"
#include"dtm_plane.h"
#include"dtm_figure.h"
#include"dtm_calculate.h"

#include"dtm_neighbor.h"

using namespace std;

namespace dtm{


	class Triangle : public Figure	//Triangle : Figure(node=3)
	{
	private:
		Neighbor neighbor[3];
		Plane plane;	//plane, determined this triangle
		int type;		//type of boundary 
						//outer surface or inner surface boundary = 1
						//domain boundary = 0
		bool visible;	//flg:use ray trace or not
		double maxl;	//maximun length
		double minl;	//minimun length
		char edge;		//widest edge number
		bool flg_pl;	//flg:plane is latest = 1, not = 0
		bool flg_lg;	//flg:length is latest = 1, not = 0

	public:
		Triangle();
		Triangle(int id);
		Triangle(int id,Node *p0,Node *p1,Node *p2);
		~Triangle();


		/* override of Element methods */

		void setId(int id);


		/* override of Figure methods */

		void setNode(int index,Node *pt);
		void renew(){
			flg_pl = false;
			flg_lg = false;
		}
		double getValue(){
			if(!flg_lg) setLength();
			return minl / maxl;
		}

		Point getCenterOfGravity();


		/* original methods */

		Triangle *resetTriangle(int id,Node *p0,Node *p1,Node *p2);
		void resetNode(Node *p0,Node *p1,Node *p2);

		void setNewNeighborElement(int index,Triangle *ptri);
		void resetNeighborElement(int index,Triangle *ptri);
		Triangle* getNeighborElement(int index,int element_number);
		void resetNeighbor(int index,Neighbor &nb);
		Neighbor& getNeighbor(int index);
		void eraseNeighborElement(int index,Triangle *ptri);
		void mergeNeighbor(int index,Neighbor &nb);
		int getNumberOfNeighborElement(int index);
		int getConection(Triangle *ptri);

		void setType(int n){ type = n; }
		int getType(){ return type;}

		void setLength();

		double getRate(){
			if(!flg_lg) setLength();
			return minl / (maxl * maxl);
		}
		double getMaxlength(){
			if(!flg_lg) setLength();
			return maxl;
		}
		double getMinlength(){
			if(!flg_lg) setLength();
			return minl;
		}
		int getEdge(){
			if(!flg_lg) setLength();
			return (int)edge;
		}		//get widethest edge number

		/*  */

		//double getEdgeAngle(int index);
	
		double getAngle(Triangle *ptri);	//get angle of each plane

		double getArea();		//get triangle area

		void reverseDirection();

		void setVisible(){ visible = true; }
		void setInvisible(){ visible = false; }
		bool isVisible(){ return visible; }

		void setPlane(){
			plane.setVector(getNode(0),getNode(1),getNode(2));
			flg_pl = true;
		}
		Plane getPlane(){ 
			if(!flg_pl) setPlane();
			return plane;
		}
		double getPlaneAnswer(Point *pt){
			if(!flg_pl) setPlane();
			return plane.getAnswer(pt);
		}
		double getPlaneVector(int index);


		/* friend function */

		friend bool operator==(Triangle &op1,Triangle &op2){
			return (op1.getId() == op2.getId());
		}

		///for debug///
		virtual void view();
		///////////////
	};


}

#endif

//////////////////////////////////////////////////////////////////EOF