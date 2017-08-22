/*

			dtmesh.h
			>Dtmesh class
			>Delauny Trianglation 3d Mesh

			last update : 2003.11.01		by t.manabe

*/
#pragma warning (disable : 4786) //for bug of Visual C++ 6.0


#ifndef DTMESH_H
#define DTMESH_H	10

#include<conio.h>
#include<cstdio>
#include<cstring>
#include<vector>

#include"dtm_calculate.h"

#include"dtm_laplas.h"
#include"dtm_crowd.h"
#include"dtm_stdio.h"
#include"dtm_tetgen.h"
#include"dtm_surface.h"

#include"dtm_timer.h"

#include"msh_lib.h"
#include"msh_mesh.h"

namespace dtm{



	//Dtmesh class
	class Dtmesh
	{
	private:

		Crowd *crowd;			//
		Surface *surface;		//
		Tetgen *tetgen;			//
 
		//int state;				//
		/*non object = 0  crowd only = 1
		   surface    = 2  tetrahedra = 3 */

		int surface_refine_type;	//used fine type
		/* not = 0  fineValue(default)= 1 fineRate= 2 fineFast = 3 */

		int surface_laplas_turn;		//number of laplasian turn in triangle
		/* default = 5 */

		int flg_crowd;			//
		int flg_surface;		//
		int flg_tetgen;			//

		int domain_size;		//
		vector< vector<int> >divdomain;

	protected:

		void initialize();
		void terminate();

	public:
		/* constructor and destructor */

		Dtmesh();
		~Dtmesh();


		/* getter and setter */

		void setMaxNumberOfNodes(unsigned int n);
		int getMaxNumberOfNodes(){ return crowd->getMaxSize(); }
		void setMaxNumberOfTriangles(unsigned int n);
		int getMaxNumberOfTriangles(){ return surface->getMaxSize(); }
		void setMaxNumberOfTetrahedras(unsigned int n);
		int getMaxNumberOfTetrahedras(){ return tetgen->getMaxSize(); }

	//	void setState(int n){ state = n; }
	//	int getState(){ return state; }

		Crowd* getCrowd(){ return crowd; }
		Surface *getPolygon(){ return surface; }
		Tetgen *getTetgen(){ return tetgen; }

		int setSurfaceRefineType(int refine_type);
		int getSurfaceRefineType(){ return surface_refine_type; }

		int setSurfaceLaplasTurn(int laplas_turn);
		int getSurfaceLaplasTurn(){ return surface_laplas_turn; }

		int setDomain(deque<int>&buff);
		int setDomain(int n,int buff[]);
		int getDomainSize();
		int getDomainSurfaceSize(int dm);
		int getDomainNumber(int dm,int index);


		/* meshing */

		int Mesh(int refine_type = 1,int laplas_turn = 5);

		int refineSurface(int refine_type = 1,int laplas_turn = 5);
		int refineTetgen();
		int meshPointToTetra();
		int meshTriangleToTetra();


		/* input and output */

		int inputDli(char *fname);
		int inputMSH(CFEMesh* m_msh);
		int inputDtm(char *fname);

		int output3dtri(char *fname);
		int output3dtet(char *fname);
		int outputAvsTriangle(char *fname);
		int outputAvsTetra(char *fname);
		int outputDliTriangle(char *fname);
		int outputDliTetra(char *fname);
		int outputRfiTriangle(char *fname,char *version);
		int outputMSHTetra(char*fname);


		///for debug
		void view();
		///////////
	};

}

#endif

//////////////////////////////////////////////////////////////////EOF