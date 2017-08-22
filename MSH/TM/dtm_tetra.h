/*

			dtm_tetra.h
			>class Tetra
			>>based Figure Class

			last update : 2003.12.17			by t.manabe

*/

#pragma warning (disable : 4786) //for bug of Visual C++ 6.0

#ifndef DTM_TETRA_H
#define DTM_TETRA_H	2

#include<iostream>

#include"dtm_error.h"
#include"dtm_node.h"
#include"dtm_figure.h"
//#include"dtm_polyhedron.h"

#include"dtm_calculate.h"


using namespace std;

namespace dtm{

	class Tetra :	public Figure	//Tetra : Figure(node=4)
	{
	private:
		Tetra* neighbor[4];	//neighbors of tetra

		double volume;		//volume of tetra
		double radius;		//out boundary sphere radius
		double value;		//value of tetra

		Point center;		//center point of out boundary sphere

		bool flg_vol;		//flg:volume is latest = 1, not = 0
		bool flg_sph;		//flg:sphere is latest = 1, not = 0
		bool flg_lg;		//flg:length is latest = 1, not = 0


	public:
		/* constructor and destructor */

		Tetra();
		Tetra(int id);
		Tetra(int id,Node *p0,Node *p1,Node *p2,Node *p3);
		~Tetra();	


		/* override of Figure methos */

		void renew(){
			flg_vol = false;
			flg_sph = false;
			flg_lg = false;
		}

		double getValue(){
			if(!flg_lg) setLength();
			return value;
		}

		Point getCenterOfGravity();


		/* original methods */


		void resetTetra(int id,Node *p0,Node *p1,Node *p2,Node *p3);
		void resetNode(Node *p0,Node *p1,Node *p2,Node *p3);

		void setSphere();
		Point* getCenter(){
			if(!flg_sph) setSphere();
			return &center;
		}
		double getRadius(){
			if(!flg_sph) setSphere();
			return radius;
		}

		void setLength();

		/*  */

		double setVolume();
		double setVolume(double vol){
			volume = vol;
			flg_vol = true;
			flg_sph = false;
			return volume;
		}
		double getVolume(){
			if(!flg_vol) setVolume();
			return volume;
		}

		void resetNeighborElement(int index,Tetra *ptet);
		Tetra* getNeighborElement(int index);
		int getConection(Tetra *ptet);


		/* friend function */

		friend bool operator==(Tetra &op1,Tetra &op2){
			return (op1.getId() == op2.getId());
		}

		///for debug///
		void view();
		void viewVolume();
		void viewSphere();
		///////////////
	};


}



#endif


//////////////////////////////////////////////////////////////////EOF