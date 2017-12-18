/*

			dtm_triangle.cpp
			>Triangle class

			last update : 2003.12.17		by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtm_triangle.h"

namespace dtm{


// Triangle
///03.11.23

Triangle::Triangle()
:Figure(ND_NUM_TRI)
{
	///
	type = 1;
	setDomain(0);
	visible = true;

	maxl = 0.;
	minl = 0.;
	edge = 0;

	flg_pl = false;
	flg_lg = false;

}


// Triangle
///03.11.23

Triangle::Triangle(int id)
: Figure(ND_NUM_TRI,id)
{
	type = 1;//////////
	visible = true;

	maxl = 0.;
	minl = 0.;
	edge = 0;

	flg_pl = false;
	flg_lg = false;

	setDomain(0);
}


// Triangle
///03.11.23

Triangle::Triangle(int id,Node *p0,Node *p1,Node *p2)
: Figure(ND_NUM_TRI,id)
{
	Figure::setNode(0,p0);
	Figure::setNode(1,p1);
	Figure::setNode(2,p2);
	type = 1;/////////
	setDomain(0);
	visible = true;

	maxl = 0.;
	minl = 0.;
	edge = 0;

	flg_pl = false;
	flg_lg = false;

	///////
	p0->insertHanger(this,FG_TYPE_TRI);
	p1->insertHanger(this,FG_TYPE_TRI);
	p2->insertHanger(this,FG_TYPE_TRI);
	/////////
}


//~Triangle
///03.12.08

Triangle::~Triangle()
{
	Node *pnd;

	///////////
	pnd = getNode(0);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(1);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(2);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
}


//setId
///03.12.08

void Triangle::setId(int id)
{
	Element::setId(id);
	Node *pnd;

	///////////
	pnd = getNode(0);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(1);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(2);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);

	pnd = getNode(0);
	if(pnd) pnd->insertHanger(this,FG_TYPE_TRI);
	pnd = getNode(1);
	if(pnd) pnd->insertHanger(this,FG_TYPE_TRI);
	pnd = getNode(2);
	if(pnd) pnd->insertHanger(this,FG_TYPE_TRI);

	return;
}


// setNode
///03.11.23

void Triangle::setNode(int index,Node *pt)
{
	if(index < 0 || index > 2) {
		Errors err(0,
			"void Triangle::setNode(Node *pt,int index)",
			"illegal index");
		throw err;
	}
	Node *pnd;
	pnd = getNode(index);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);

	Figure::setNode(index,pt);

	flg_pl = false;
	flg_lg = false;

	//////////
	pt->insertHanger(this,FG_TYPE_TRI);
	////////

	return;
}


// resetTriangle
///03.11.23

Triangle* Triangle::resetTriangle(int id,Node *p0,Node *p1,Node *p2)
{
	Node *pnd;

	///////////
	pnd = getNode(0);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(1);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(2);
	if(pnd) pnd->eraseHanger(getId(),FG_TYPE_TRI);
	///////////

	setId(id);
	setState(0);
	Figure::setNode(0,p0);
	Figure::setNode(1,p1);
	Figure::setNode(2,p2);
	neighbor[0].clearElement();
	neighbor[1].clearElement();
	neighbor[2].clearElement();


	maxl = 0.;
	minl = 0.;
	edge = 0;

	flg_pl = false;
	flg_lg = false;

	///////
	p0->insertHanger(this,FG_TYPE_TRI);
	p1->insertHanger(this,FG_TYPE_TRI);
	p2->insertHanger(this,FG_TYPE_TRI);
	/////////

	return this;
}


// resetNode
///03.11.23

void Triangle::resetNode(Node *p0,Node *p1,Node *p2)
{
	Node *pnd;

	///////////
	pnd = getNode(0);
	pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(1);
	pnd->eraseHanger(getId(),FG_TYPE_TRI);
	pnd = getNode(2);
	pnd->eraseHanger(getId(),FG_TYPE_TRI);
	////////////////

	setState(0);
	Figure::setNode(0,p0);
	Figure::setNode(1,p1);
	Figure::setNode(2,p2);
	neighbor[0].clearElement();
	neighbor[1].clearElement();
	neighbor[2].clearElement();


	maxl = 0.;
	minl = 0.;
	edge = 0;

	flg_pl = false;
	flg_lg = false;

	///////
	p0->insertHanger(this,FG_TYPE_TRI);
	p1->insertHanger(this,FG_TYPE_TRI);
	p2->insertHanger(this,FG_TYPE_TRI);
	/////////

	return;
}


//setNeighborElement
///03.11.30

void Triangle::setNewNeighborElement(int index,Triangle *ptri)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"void Triangle::setNewNeighborElement(int index,Triangle *ptri)",
			"illegal index");
		throw err;
	}
	neighbor[index].pushElement(ptri);
}


//resetNeighborElement
///03.12.03

void Triangle::resetNeighborElement(int index,Triangle *ptri)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"void Triangle::resetNeighborElement(int index,Triangle *ptri)",
			"illegal index");
		throw err;
	}
	neighbor[index].clearElement();
	neighbor[index].pushElement(ptri);
}


//getNeighborElement
///03.12.03

Triangle* Triangle::getNeighborElement(int index,int element_number)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"Triangle* Triangle::getNeighborElement(int index,int element_number)",
			"illegal index");
		throw err;
	}
	return neighbor[index].getElement(element_number);
}


//resetNeitghor
///03.12.03

void Triangle::resetNeighbor(int index,Neighbor &nb)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"void Triangle::resetNeighbor(int index,Neighbor &nb)",
			"illegal index");
		throw err;
	}
	neighbor[index] = nb;
}


//getNeighbor
///03.12.03

Neighbor& Triangle::getNeighbor(int index)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"Neighbor& Triangle::getNeighbor(int index)",
			"illegal index");
		throw err;
	}
	return neighbor[index];
}


//eraseNeighborElement
///03.12.03

void Triangle::eraseNeighborElement(int index,Triangle *ptri)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"void Triangle::eraseNeighborElement(int index,Triangle *ptri)",
			"illegal index");
		throw err;
	}
	neighbor[index].eraseElement(ptri);
}


//mergeNeighbor
//03.12.10

void Triangle::mergeNeighbor(int index,Neighbor &nb)
{
	if(index > 3 || index < 0) {
		Errors err(0,
			"void Triangle::mergeNeighbor(int index,Neighbor &nb)",
			"illegal index");
		throw err;
	}
	int n = nb.getNumberOfElement();
	for(int i=0;i<n;i++) {
		setNewNeighborElement(index,nb.getElement(i));
	}

	return;
}


//getNumberOfNeighborElement
///03.12.03

int Triangle::getNumberOfNeighborElement(int index)
{
	if(index >= 3 || index < 0) {
		Errors err(0,
			"int Triangle::getNumberOfNeighborElement(int index)",
			"illegal index");
		throw err;
	}
	return neighbor[index].getNumberOfElement();
}


// getConection
///03.12.03

int Triangle::getConection(Triangle *ptri)
{
	Triangle *itri;
	int num;

	for(int i=0;i<3;i++) {
		num = neighbor[i].getNumberOfElement();
		for(int j=0;j<num;j++) {
			itri = neighbor[i].getElement(j);
			if(*ptri == *itri) return i;
		}
	}

	/*
	Errors err(0,
		"int Triangle::getConection(Triangle *ptri)",
		"cannot get conection");
	throw err;
	*/

	return -1;
}


// setLength
///03.11.23

void Triangle::setLength()
{
	double l0,l1,l2;

	l0 = getDistance(getNode(1),getNode(2));
	l1 = getDistance(getNode(2),getNode(0));
	l2 = getDistance(getNode(0),getNode(1));
	if(l0 > l1) {
		maxl = l0;
		minl = l1;
		edge = 0;
	} else {
		maxl = l1;
		minl = l0;
		edge = 1;
	}
	if(l2 > maxl) {
		maxl = l2;
		edge = 2;
	} else if(l2 < minl) {
		minl = l2;
	}

	flg_lg = true;

	return;
}


Point Triangle::getCenterOfGravity()
{
	return getTriangleCenter(this);
}


// getArea
///03.11.23

double Triangle::getArea()
{
	return getTriangleArea(getNode(0),getNode(1),getNode(2));
}



// getAngle
///03.11.24

double Triangle::getAngle(Triangle *ptri)
{
	double a,b,c,d;
	double dist;
	double cosa,cosb,cosc;
	double cosa1,cosb1,cosc1;
	double cosn;


	if(!flg_pl) setPlane();

	
	a = getPlaneVector(0);
	b = getPlaneVector(1);
	c = getPlaneVector(2);
	d = getPlaneVector(3);
	dist = a*a+b*b+c*c;
	dist = sqrt(dist);
	cosa = a/dist;
	cosb = b/dist;
	cosc = c/dist;
	a = ptri->getPlaneVector(0);
	b = ptri->getPlaneVector(1);
	c = ptri->getPlaneVector(2);
	d = ptri->getPlaneVector(3);
	dist = a*a+b*b+c*c;
	dist = sqrt(dist);
	cosa1 = a/dist;
	cosb1 = b/dist;
	cosc1 = c/dist;

	cosn = cosa*cosa1+cosb*cosb1+cosc*cosc1;

	return cosn;
}


// getPlaneVector
///03.11.23

double Triangle::getPlaneVector(int index)
{
	if(index < 0 || index > 3) {
		Errors err(0,
			"double Triangle::getPlaneVector(int index)",
			"illegal index");
		throw err;
	}
	if(!flg_pl) setPlane();

	return plane.getVector(index);
}


//reverseDirection
///03.12.03

void Triangle::reverseDirection()
{
	Node *p1,*p2;
	Neighbor t1,t2;

	p1 = getNode(1);
	p2 = getNode(2);
	t1 = neighbor[1];
	t2 = neighbor[2];

	Figure::setNode(1,p2);
	Figure::setNode(2,p1);
	neighbor[1] = t2;
	neighbor[2] = t1;

	flg_pl = false;
	flg_lg = false;

	return;
}



///for debug///
// view

void Triangle::view()
{
	Triangle* ptri;

	cout << "TRI ID " << getId();
	cout << " ST " << getState();
	cout << " TP " << type;
	cout << " DM " << getDomain();
	cout << "  MTJ";
	for(int i=0;i<3;i++) {
		if(!getNode(i)) {
			cout << " 0";
			continue;
		}
		cout << " " << getNode(i)->getId();
	}
	
	cout << "  JAC";
	for(int i=0;i<3;i++) {
		int n = getNumberOfNeighborElement(i);
		cout << " [";
		for(int j=0;j<n;j++) {
			ptri = getNeighborElement(i,j);
			if(!ptri) {
				cout << " 0";
				continue;
			}
			else cout << " " << ptri->getId();
		}
		cout << "]";
	}
	
	cout << "\n";

	return;
}


}
//////////////////////////////////////////////////////////////////EOF