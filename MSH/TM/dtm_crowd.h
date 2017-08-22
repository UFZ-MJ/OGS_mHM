/*

			dtm_crowd.h
			>Crowd class
			>>based GpNode

			last update : 2003.12.08		by t.manabe

*/


#ifndef DTM_CROWD_H
#define DTM_CROWD_H	3

#include<iostream>
#include<string>
#include<deque>

#include"dtm_point.h"
#include"dtm_node.h"
#include"dtm_laplas.h"

#include"dtm_gpnode.h"

using namespace std;

namespace dtm{



	class Crowd		:	public GpNode
	{
	private:
		deque<int>binlist;	//

		bool flg_std;		//flg:be normalized = 1, not = 0
		bool flg_wid;		//flg:every length are latest = 1, not = 0
		bool flg_den;		//

		Point stdpt;		//standard point for normalize(= -(minpt))
		double stdrate;		//standard rate for normalize(= maxwidth)

		Point maxpt;		//maximum point(maxX,maxY,maxZ)!=node
		Point minpt;		//minimum point(minX,minY,minZ)!=node
		double xwidth;		//maximum length of X direction
		double ywidth;		//maximum length of Y direction
		double zwidth;		//maximum length of Z direction
		double maxwidth;	//maximum length of each direction
		double minwidth;	//minimum length of each direction
		double maxdensity;	//maximum density
		double mindensity;	//minimum density
		double avedensity;	//average of density
		Node* maxaxispt[3];	//maximum node of each direction [0]:x,[1]:y,[2]:z
		Node* minaxispt[3];	//minimum node of each direction [0]:x,[1]:y,[2]:z


	public:
		/* constructor and destructor */

		Crowd();
		~Crowd();


		/* override of GpNode methods*/

		int erase(int id);
		void clear();
		Node* makeNewNode(double x,double y,double z,double density);


		/* original methods */

		virtual void initstd();	//must be called after input Node
		void initDensity();		//must be called after input Node`s density

		bool isStd(){ return flg_std; }

		double getStdRate();
		Point getStdPoint();
		Point getMaxPoint();
		Point getMinPoint();
		double getXWidth();
		double getYWidth();
		double getZWidth();
		double getMaxWidth();
		double getMinWidth();
		double getMaxDensity();
		double getMinDensity();
		double getAveDensity();
		Node* getMaxXAxisNode();
		Node* getMinXAxisNode();
		Node* getMaxYAxisNode();
		Node* getMinYAxisNode();
		Node* getMaxZAxisNode();
		Node* getMinZAxisNode();
		Point getCurrentMaxPoint();
		Point getCurrentMinPoint();
		double getCurrentXWidth();
		double getCurrentYWidth();
		double getCurrentZWidth();
		double getCurrentMaxWidth();
		double getCurrentMinWidth();
		double getCurrentMaxDensity();
		double getCurrentMinDensity();
		double getCurrentAveDensity();
		
		void normalize();	//normalization
		void restore();		//restore to original values
		void binsort();		//binsorting of nodes
		int getBinListNumber(int index);

		void laplasianTriangle(int first,int end);
		//virtual void laplasianTetra(int first,int end);

		void setDensity(double idensity,int flg = 1);

		void clearHanger(int fgType = -1);

		///for debug////
	//	void view();
		void viewHanger();
		////////////////
	};


}


#endif

//////////////////////////////////////////////////////////////////EOF