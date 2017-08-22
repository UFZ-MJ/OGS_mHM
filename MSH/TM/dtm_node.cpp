/*

			dtm_node.cpp
			>Node class

			last update : 2003.12.03		by t.manabe

*/


#include "stdafx.h" /* MFC */


#include"dtm_node.h"

namespace dtm{

// Node
///03.11.27

Node::Node()
: cr_density(0.) , type(0) //, flg_density(false)
{
	hanger.resize(NUM_FG_TYPE);
}


//Node
///03.11.27

Node::Node(int n)
: Element(n)
{
	cr_density = 0.;
	type = 0;
//	flg_density = false;
	hanger.resize(NUM_FG_TYPE);
}


//Node
///03.11.27

Node::Node(int n,double x,double y,double z)
: Element(n) , Point(x,y,z)
{
	cr_density = 0.;
	type = 0;
	//flg_density = false;
	hanger.resize(NUM_FG_TYPE);
}


//Node
///03.11.27

Node::Node(int n,double x,double y,double z,double density)
: Element(n) , Point(x,y,z)
{
	cr_density = density;
	type = 0;
//	flg_density = true;
	hanger.resize(NUM_FG_TYPE);
}


//Node
///03.11.27

Node::Node(int n,Point ppt,double density)
: Element(n) , Point(ppt)
{
	cr_density = density;
	type = 0;
//	flg_density = true;
	hanger.resize(NUM_FG_TYPE);
}


//Node
///03.11.27

Node::Node(int n,Point ppt)
: Element(n) , Point(ppt)
{
	type = 0;
//	flg_density = false;
	hanger.resize(NUM_FG_TYPE);
}


//getHanger
///03.11.23

Figure* Node::getHanger(int index,int fgType)
{
	int n = (int)hanger[fgType].size();
	if(index < 0 || index >= n) {
		Errors err(0,
			"int Node::getHanger(int index)",
			"illegal index");
		throw err;
	}

	return hanger[fgType][index];
}


//insertHanger
///03.11.23

void Node::insertHanger(Figure *hang,int fgType)
{
	int nid;
	vector<Figure*>::iterator p = hanger[fgType].begin();

	nid = hang->getId();
	int n = (int)hanger[fgType].size();
	for(int i=0;i<n;i++) {
		if(*hang == *(hanger[fgType][i]) ) {
			return;
		}else if(nid < hanger[fgType][i]->getId()) {
			p += i;
			hanger[fgType].insert(p,hang);
			return;
		}
	}
	hanger[fgType].push_back(hang);

	return;
}


//eraseHanger
///03.11.23

void Node::eraseHanger(int id,int fgType)
{
	if(fgType < 0 || fgType >= NUM_FG_TYPE) {
		Errors err(0,
			"void Node::eraseHanger(int id,int fgType)",
			"illegal figure type");
		throw err;
	}

	vector<Figure*>::iterator p = hanger[fgType].begin();

	int n = (int)hanger[fgType].size();
	for(int i=0;i<n;i++) {
		if(id == hanger[fgType][i]->getId()) {
			p += i;
			hanger[fgType].erase(p);
			return;
		}
	}

	return;
}


//clearHanger
///03.11.23

void Node::clearHanger(int fgType)
{
	int i;
	if(fgType < 0) {
		for(i=0;i<NUM_FG_TYPE;i++) hanger[i].clear(); 
	} else if(fgType >= NUM_FG_TYPE) {
		Errors err(0,
			"void Node::clearHanger(int fgType)",
			"illegal figure type");
		throw err;
	} else {
		hanger[fgType].clear(); 
	}
}


///for debug///
// viewHanger

void Node::viewHanger(int fgType)
{
	cout << "NODE ID = " << getId();
	cout << " ";
	int n = (int)hanger[fgType].size();
	for(int i=0;i<n;i++) {
		cout << " " << hanger[fgType][i]->getId();
	}
	cout << "\n";
/*	*/
	return;
}


// view

void Node::view()
{
	cout << "NODE ID " << getId();
	cout << " TP " << (int)type;
	cout << "  x " << getX() << "  y " << getY() << "  z " << getZ() << "  den "
		<< cr_density << "\n";
	return;
}


}

//////////////////////////////////////////////////////////////////EOF