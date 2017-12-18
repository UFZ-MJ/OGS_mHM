/*

			dtm_tetgen.h
			>Tetgen class
			>>based GpTetra
			>>create new tetrahedras on 3d delauny trianglation

			last update : 2003.12.08		by t.manabe

*/

#pragma warning (disable : 4786) //for bug of Visual C++ 6.0

#ifndef DTM_TETGEN_H
#define DTM_TETGEN_H 3

#include<iostream>
#include<deque>
#include<vector>
#include<stack>
#include<algorithm>
#include<functional>

#include"dtm_fixed.h"
#include"dtm_node.h"
#include"dtm_plane.h"
#include"dtm_tetra.h"
#include"dtm_crowd.h"
#include"dtm_calculate.h"

#include"dtm_gptetra.h"


using namespace std;

namespace dtm{


	class Tetgen	:	public GpTetra
	{
	private:
		double rate_minvolm;	//rate of minimum volume
		double rate_density;	//rate of density

		Node* sptet_nd[4];		//nodes of supertetrahedron

		//Adpt struct
		///used polygon(!= class Polygon) to tetra
		///for getting high speed of program
		struct Adpt{
			Node *face[3];
			Tetra *neib;
			int num;
			double vol;
			Adpt() {
				face[0] = NULL; face[1] = NULL;
				face[2] = NULL; neib = NULL; }
			Adpt(const Adpt& ad){
				face[0] = ad.face[0];
				face[1] = ad.face[1];
				face[2] = ad.face[2];
				neib = ad.neib; num = ad.num; vol = ad.vol; }
		};
		///hash table, used seach neighbor
		int tbl_long;				//table length
		int flg_tbl;				//flg:initialized = 1, not = 0
		vector<vector<Adpt> > tbl;	//hash table

		Crowd *cr_crowd;

	protected:

		/* for generator(delauny) and finTetra*/

		int getPolyTetra(vector<int>&tet_id_buff,Tetra *loc,
			Node *cr_nd);
		int correctPolyTetraAdd(vector<int>&tet_id_buff,
			vector<Adpt>&adpt_buff,Node *cr_nd);
		int correctPolyTetraCut(vector<int>&tet_id_buff,
			vector<Adpt>&adpt_buff,Node *cr_nd);
		int reMeshTetra(vector<int>&tet_id_buff,
			vector<Adpt>&adpt_buff,Node *cr_nd);
		int conectNewTetra(vector<int>&tet_id_buff,int nend);
		//add node in tetra(loc)
		int addNode(vector<int>&tetra_buff,Tetra *loc,Node *add);


	public:
		/* constructor and destructor */

		Tetgen();
		Tetgen(Crowd *crowd);
		~Tetgen();


		/* original methods */

		void setTetgenParameter(double minvol = 0.5,
			double rtdens = 0.8,int tbllg = 20);
		bool setCurrentCrowd(Crowd *crowd);
		void setDomain(int dm);
		void setTetHanger();

		Tetra* locate(Point *ppt);
		Tetra* optlocate(Point *ppt);

		void delaunay(Crowd *crowd);
		void delaunay();
		void generator(Node *cd_nd);

		void makeSuperTetra(double maxwidth,Point minpt);
		void makeSuperTetra();
		void removeSuperTetra();
		void removeTetrahedra();

		void fineTetra();

		void cleanBadTetra(double min_value);

		Node* eraseNode(Node *er_nd);
		Node* mergeNode(Node *po_nd,Node *ne_nd);

		void divHanger(Node *po_nd,Node *ne_nd,deque<Tetra*>&pohang,
						deque<Tetra*>&nehang);
		int isPossibleMerge(Node *po_nd,Node *ne_nd,deque<Tetra*>&pohang,
							deque<Tetra*>&nehang);
		///for debug///
		void view();
		///////////////
	};



}





#endif

//////////////////////////////////////////////////////////////////EOF