/*

			dtm_tetra.cpp
			>Tetra class

			last update : 2003.12.03		by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtm_tetra.h"

namespace dtm{




// Tetra
///03.11.23

Tetra::Tetra()
:Figure(ND_NUM_TET)
{
	setDomain(0);

	neighbor[0] = NULL;
	neighbor[1] = NULL;
	neighbor[2] = NULL;
	neighbor[3] = NULL;
	volume = 0.;
	radius = 0.;
	value = 0.;

	flg_vol = false;
	flg_sph = false;
	flg_lg = false;
}


// Tetra
///03.11.23

Tetra::Tetra(int id)
:Figure(ND_NUM_TET,id)
{
	setDomain(0);

	neighbor[0] = NULL;
	neighbor[1] = NULL;
	neighbor[2] = NULL;
	neighbor[3] = NULL;
	volume = 0.;
	radius = 0.;
	value = 0.;

	flg_vol = false;
	flg_sph = false;
	flg_lg = false;
}


// Tetra
///03.11.23

Tetra::Tetra(int id,Node *p0,Node *p1,Node *p2,Node *p3)
: Figure(ND_NUM_TET,id)
{
	Figure::setNode(0,p0);
	Figure::setNode(1,p1);
	Figure::setNode(2,p2);
	Figure::setNode(3,p3);

	neighbor[0] = NULL;
	neighbor[1] = NULL;
	neighbor[2] = NULL;
	neighbor[3] = NULL;

	setDomain(0);

	volume = 0.;
	radius = 0.;
	value = 0.;

	flg_vol = false;
	flg_sph = false;
	flg_lg = false;

}


//~Tetra

Tetra::~Tetra()
{
	Node *pnd;
	for(int i=0;i<4;i++) {
		pnd = getNode(i);
		if(!pnd)continue;
		pnd->eraseHanger(getId(),FG_TYPE_TET);
	}
}


//getCenterOfGravity
///


Point Tetra::getCenterOfGravity()
{
	return getTetraCenter(this);
}



// setTetra
///03.11.23

void Tetra::resetTetra(int id,Node *p0,Node *p1,Node *p2,Node *p3)
{
	setId(id);
	setState(0);
	Figure::setNode(0,p0);
	Figure::setNode(1,p1);
	Figure::setNode(2,p2);
	Figure::setNode(3,p3);
	neighbor[0] = NULL;
	neighbor[1] = NULL;
	neighbor[2] = NULL;
	neighbor[3] = NULL;

	volume = 0.;
	radius = 0.;
	value = 0.;

	flg_vol = false;
	flg_sph = false;
	flg_lg = false;


	return;
}


// setNode	
///03.11.23

void Tetra::resetNode(Node *p0,Node *p1,Node *p2,Node* p3)
{
	setState(0);
	Figure::setNode(0,p0);
	Figure::setNode(1,p1);
	Figure::setNode(2,p2);
	Figure::setNode(3,p3);
	neighbor[0] = NULL;
	neighbor[1] = NULL;
	neighbor[2] = NULL;
	neighbor[3] = NULL;

	volume = 0.;
	radius = 0.;
	value = 0.;

	flg_vol = false;
	flg_sph = false;
	flg_lg = false;

}


//getNeighorElement
///03.12.08

Tetra* Tetra::getNeighborElement(int index)
{
	if(index >= 4 || index < 0) {
		Errors err(0,
			"Tetra* Tetra::getNeighborElement(int index)",
			"illegal index");
		throw err;
	}
	return neighbor[index];
}


//resetNeighborElement
//03.12.08

void Tetra::resetNeighborElement(int index,Tetra *ptet)
{
	if(index >= 4 || index < 0) {
		Errors err(0,
			"void Tetra::resetNeighborElement(int index,Tetra *ptet)",
			"illegal index");
		throw err;
	}
	neighbor[index] = ptet;
	return;
}


// getConection
///03.11.23

int Tetra::getConection(Tetra *ptet)
{

	if(neighbor[0])if(*neighbor[0] == *ptet) return 0;
	if(neighbor[1])if(*neighbor[1] == *ptet) return 1;
	if(neighbor[2])if(*neighbor[2] == *ptet) return 2;
	if(neighbor[3])if(*neighbor[3] == *ptet) return 3;
	/*	
	Errors err(0,
		"int Tetra::getConection(Tetra *ptet)",
		"cannot get conection");
	throw err;
	*/

	return -1;
}


// setVolume
///03.11.23

double Tetra::setVolume()
{
    volume = getTetraVolume(getNode(0),getNode(1),getNode(2),getNode(3));

	flg_vol = true;

	return volume;
}


// setSphere
///03.11.23

void Tetra::setSphere()
{
	if(!flg_vol) setVolume();
	center = getTetraSphereCenter(getNode(0),getNode(1),getNode(2),
		getNode(3),volume);
    radius = getDistance(&center,getNode(0))*(1. + LIMIT);

	flg_sph = true;

	return;
}


// setLength
///03.11.23

void Tetra::setLength()
{
	double length[6];
	double minl,maxl;
	const int tl[6][2] = {0,1,0,2,0,3,1,2,1,3,2,3};

	for(int i=0;i<6;i++) {
		length[i] = getDistance( getNode(tl[i][0]) ,getNode(tl[i][1]) );
	}
	maxl = minl = length[0];
	if(length[1] > maxl)maxl = length[1];
	else if(length[1] < minl) minl = length[1];
	if(length[2] > maxl)maxl = length[2];
	else if(length[2] < minl) minl = length[2];
	if(length[3] > maxl)maxl = length[3];
	else if(length[3] < minl) minl = length[3];
	if(length[4] > maxl)maxl = length[4];
	else if(length[4] < minl) minl = length[4];
	if(length[5] > maxl)maxl = length[5];
	else if(length[5] < minl) minl = length[5];


	value = minl / maxl;

	flg_lg = true;

	return;
}


///for debug///
// view

void Tetra::view()
{
	Tetra* ptet;

	cout << "TET ID " << getId();
	cout << " ST " << getState();
	cout << "  MTJ";
	for(int i=0;i<4;i++) {
		if(!getNode(i)) {
			cout << " 0";
			continue;
		}
		cout << " " << getNode(i)->getId();
	}
	cout << "  JAC";
	for(int i=0;i<4;i++) {
		ptet = getNeighborElement(i);
		if(!ptet) {
			cout << " 0";
			continue;
		}
		cout << " " << ptet->getId();
	}
	cout << "\n";

	return;
}


// viewVolume

void Tetra::viewVolume()
{
	if(!flg_vol) setVolume();

	cout << "VOLUME : TETRA ID = " << getId();
	cout << " VOL = " << volume << "\n";

	return;
}


// viewSphere

void Tetra::viewSphere()
{
	if(!flg_sph) setSphere();

	cout << "SPHERE : TETRA ID = " << getId();
	cout << " RAD = " << radius;
	center.view();

	return;
}


}


//////////////////////////////////////////////////////////////////EOF