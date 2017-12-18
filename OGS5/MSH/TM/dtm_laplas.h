/*

			dtm_laplas.h
			>laplasian of node

			last update : 2003.11.13			by t.manabe

*/


#ifndef DTM_LAPLAS_H
#define DTM_LAPLAS_H	2

#include"dtm_node.h"
#include"dtm_triangle.h"
#include"dtm_tetra.h"

#include<math.h>
#include<cstdio>

namespace dtm{

	void laplasInTri(Node *nd);
	void laplasInTriEdge(Node *nd);
	void laplasInTet(Node *nd);

}

#endif


//////////////////////////////////////////////////////////////////EOF