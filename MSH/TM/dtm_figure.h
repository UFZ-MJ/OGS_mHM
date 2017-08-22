/*

			dtm_figure.h
			>FIgure class
			>>based Elememt class
			>>has nodes
			>>must be fixed the number of nodes when initialized

			last update : 2003.12.17			by t.manabe

  */

#pragma warning (disable : 4786) //for bug of Visual C++ 6.0

#ifndef DTM_FIGURE_H
#define DTM_FIGURE_H	1

#include"dtm_error.h"
#include"dtm_element.h"
#include"dtm_point.h"

namespace dtm{

	class Node;

	class Figure : public Element		//base class of all figure
	{
	private:
		const unsigned int number_of_node;
		Node **figure_node;				//nodes of figure
		int state;						//state of object
		int domain;						//domain of figure

	public:
		/* constructor and destructor */

		Figure(unsigned int nNode);
		Figure(unsigned int nNode,int id);
		~Figure();


		/* original methos */

		virtual void setNode(int index,Node *nd);
		virtual Node* getNode(int index);
		int getNumberOfNode(){ return number_of_node; }

		virtual int getNodeNumber(Node *nd);

		virtual void setState(int st){state = st;}
		virtual int getState(){return state;}

		virtual double getValue() = 0;

		virtual Point getCenterOfGravity() = 0;

		virtual void renew() = 0;


		virtual void setDomain(int dm) { domain = dm; }
		virtual int getDomain() { return domain; }


		/* friend function */

		friend bool operator==(Figure &op1,Figure &op2) {
			return op1.getId() == op2.getId();
		}

	};



}


#endif

//////////////////////////////////////////////////////////////////EOF