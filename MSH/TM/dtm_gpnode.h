/*

			dtm_gpnode.h
			>GpNode Class
			>>based Group Class
			>>manage all nodes

			last update  : 2003.12.08			by t.manabe

  */

#ifndef DTM_GPNODE_H
#define DTM_GPNODE_H	1

#include<deque>

using namespace std;

#include"dtm_error.h"
#include"dtm_group.h"
#include"dtm_node.h"

namespace dtm{


	class GpNode : public Group
	{
		deque<Node*>node;

	public:
		/* constructor and destructor */

		GpNode();
		~GpNode();


		/* override of Group methods */

		virtual int getSize();
		virtual void clear();
		virtual int erase(int id);


		/* original methods */

		virtual Node* getNode(int index);
		virtual Node* makeNewNode(double x,double y,double z,double density);


		//
		virtual void view();

	};



}

#endif

//////////////////////////////////////////////////////////////////EOF