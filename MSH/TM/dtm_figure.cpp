/*

			dtm_figure.cpp
			>Figure class

			last update : 2003.12.17			by t.manabe

  */


#include "stdafx.h" /* MFC */


#include"dtm_figure.h"

namespace dtm{


//Figure
///03.11.30

Figure::Figure(unsigned int nNode)
: number_of_node(nNode), state(0), domain(0)
{
	figure_node =  new Node*[number_of_node];
	if(!figure_node) {
		Errors err(-1,
			"Figure::Figure(unsigned int nNode,unsigned int nNeighbor)",
			"failed making new object");
		throw err;
	}
	int n = number_of_node;
	for(int i=0;i<n;i++) figure_node[i] = NULL;
}

//Figure
///03.11.30

Figure::Figure(unsigned int nNode,int id)
: Element(id),
number_of_node(nNode), state(0), domain(0)
{
	figure_node =  new Node*[number_of_node];
	if(!figure_node) {
		Errors err(-1,
			"Figure::Figure(unsigned int nNode,unsigned int nNeighbor,int id,int st)",
			"failed making new object");
		throw err;
	}
	int n = number_of_node;
	for(int i=0;i<n;i++) figure_node[i] = NULL;
}


//~Figure
///03.11.23

Figure::~Figure()
{
	delete [] figure_node;
}


//setNode
///03.11.23

void Figure::setNode(int index,Node *nd)
{
	int n = number_of_node;
	if(index < 0 || index >= n) {
		Errors err(0,
			"void Figure::setNode(int index,Node *pt)",
			"illegal index");
		throw err;
	}

	figure_node[index] = nd;
}


//getNode
///03.11.23

Node* Figure::getNode(int index)
{
	int n = number_of_node;
	if(index < 0 || index >= n) {
		Errors err(0,
			"Node* Figure::getNode(int index)",
			"illegal index");
		throw err;
	}
	return figure_node[index];
}


//getNodeNumber
///03.12.15

int Figure::getNodeNumber(Node *nd)
{
	if(!nd) return -1;

	int n = number_of_node;
	for(int i=0;i<n;i++) {
		if(nd == (figure_node[i])) {
			return i;
		}
	}

	Errors err(0,
		"int Figure::getNodeNumber(Node *nd)",
		"nd is not used");
	throw err;

	//return -1;
}



}

//////////////////////////////////////////////////////////////////EOF