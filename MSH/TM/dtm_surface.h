/*

			dtm_surface.h
			>Surface class
			>>based Polygon class
			>>create new triangles on 2d(2.5d) delauny trianglation

			last update : 2003.12.08		by t.manabe

*/

#pragma warning (disable : 4786) //for bug of Visual C++ 6.0

#ifndef DTM_SURFACE_H
#define DTM_SURFACE_H	4

#include<iostream>
#include<string>
#include<deque>
#include<stack>
#include<vector>
#include<algorithm>
#include<math.h>

#include"dtm_fixed.h"

#include"dtm_node.h"
#include"dtm_plane.h"
#include"dtm_triangle.h"
#include"dtm_crowd.h"
#include"dtm_laplas.h"
#include"dtm_calculate.h"
#include"dtm_gptriangle.h"
#include"dtm_tetgen.h"

using namespace std;

namespace dtm{


	class Surface : public GpTriangle		//Surface class, derived Polygon class
	{
	private :
		int flg_lg;			//flg:length is latest = 1, not = 0

		double maxlength;	//maximum length of triangle side
		double minlength;	//minimum length of triangle side
		double maxaverage;	//average of maximum length
		double minaverage;	//average of minimum length

		Crowd *cr_crowd;

		int laplas_flg;			//do laplasian in fine = 1, not = 0
		double rate_density;	//permited rate of density
		double rate_minarea;	//permited rate of minimum area
		int stack_maxlong;		//max long of stack(using fine)

		vector<Point>center;	//buffer of center points of triangle
		vector<Triangle*>focus;	//buffer of triangles which used for judgeInOut
		deque<Triangle*>hittri;	//buffer of triangles which used for checkDirection

	protected:

		/* for selectDomain */

		int pickTriangle(Tetra *tet);
		int pickLine(Tetra *tet);
		int pickPoint(Tetra *tet);
		//judgement of tetra, if in domain or not
		int judgeInOut(Tetra *tet);


		/* for sortDirection*/

		void divideSection();
		void setDirection();
		int checkSamePlane(deque<Triangle*>&stktri,deque<double>&stkdist,
			Point &lmtpt);
		int correctTriangleDirection(deque<Triangle*>&que);


		/* for reMesh*/

		void qSortValue(int ,int e,deque<Triangle*>&array);
		void qSortRate(int ,int e,deque<Triangle*>&array);
		void qSortDist(int f,int e,deque<double>&dist,deque<Triangle*>&array);
		//remeshing of opposite side
		int remeshOppositeSide(Neighbor &neib,Triangle *ielm,Triangle *tt1,
			Node *iv1,Node *iv2,Node *iv3,
			deque<Triangle*>&istack,Node *ad);
		//swapping of triangles
		int swapTriangle(deque<Triangle*>&istack,Node *ad,double minarea);
		//create new node on widest side
		Node* createNode(Triangle* tri,int edge);
		int clearnStack(deque<Triangle*>&istack,Node *ad);
		int pushStack(deque<Triangle*>&istack,double minarea);


		/* for mergeNode */
		
		void divHanger(Node* po_nd,Node* ne_nd,deque<Triangle*>&pohang,
			deque<Triangle*>&nehang);
		int isPossibleMerge(Node *po_nd,Node *ne_nd,deque<Triangle*>&pohang,
			deque<Triangle*>&nehang);
		void swapLawson(deque<Triangle*>istack);



	public :
		/* constructor and destructor */

		Surface();
		Surface(Crowd *crowd);
		~Surface();


		/* original methods */

		virtual void setLength();
		virtual double getMaxlength();
		virtual double getMinlength();
		virtual double getMaxaverage();
		virtual double getMinaverage();

		void setCurrentCrowd(Crowd *crowd);

		virtual void setNeighbor();

		void setSurfaceParameter(int lplsflg = 1,double ratednst = 0.9,
			double ratearea = 0.05,int stklong = 1000);

		void pickDomain(vector<int>&area);
		void pickBoundary();

		void refineOnValue();
		void refineOnRate();
		void refineOnFast();
		void reMesh(Node* ad,Triangle* ielm,int edge);

		void sortDirection();
		void selectDomain(Tetgen *tetgen);

		void cleanBadTriangle(double min_value);
		Node* eraseNode(Node* nd);
		Node* mergeNode(Node* po_nd,Node* ne_nd);


	};

}

#endif

//////////////////////////////////////////////////////////////////EOF