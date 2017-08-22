/*

			dtm_gptetra.h
			>GpTetra Class
			>>based Group Class
			>>manage all tetrahedra elements

			last update  : 2003.12.08			by t.manabe

  */

#ifndef DTM_GPTETRA_H
#define DTM_GPTETRA_H	2

#include<deque>

using namespace std;

#include"dtm_tetra.h"
#include"dtm_group.h"

namespace dtm{


	class GpTetra	:	public Group
	{
		deque<Tetra*>tetra;

	public:
		/* constructor and destructor */

		GpTetra();
		~GpTetra();


		/* override of Group methods */

		int getSize();
		void clear();
		int erase(int id);


		/* original methos */

		virtual Tetra* getTetra(int index);
		virtual Tetra* makeNewTetra(Node *nd1,Node *nd2,Node *nd3,Node *nd4);


	};




}

#endif

//////////////////////////////////////////////////////////////////EOF