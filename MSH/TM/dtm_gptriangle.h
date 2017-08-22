/*

			dtm_gptriangle.h
			>GpTriangle Class
			>>based Group Class
			>>manage all triangle elements

			last update : 2003.12.08			by t.manabe

  */

#ifndef DTM_GPTRIANGLE_H
#define DTM_GPTRIANGLE_H	2

#include<deque>

using namespace std;

#include"dtm_node.h"
#include"dtm_triangle.h"
#include"dtm_group.h"

namespace dtm{


	class GpTriangle	:	public Group
	{
	private:
		deque<Triangle*>triangle;

	public:
		/* constructor and destructor */

		GpTriangle();
		~GpTriangle();


		/* override of Group methods */

		virtual int getSize();
		virtual void clear();
		virtual int erase(int id);


		/* original meshods */

		virtual Triangle* getTriangle(int index);
		virtual Triangle* makeNewTriangle(Node *nd1,Node *nd2,Node *nd3);

		//
		virtual void view();
	};

}

#endif

//////////////////////////////////////////////////////////////////EOF