/*

			dtm_crowd.cpp
			>Crowd class

			last update : 2003.12.08		by t.manabe

*/


#include "stdafx.h" /* MFC */


#include"dtm_crowd.h"

namespace dtm{


// Crowd
///03.11.23

Crowd::Crowd()
: flg_std(false) , flg_wid(false) ,flg_den(false)
{

	stdpt.setCoordinates(0.,0.,0.);
	stdrate = 1.;
	maxpt.setCoordinates(0.,0.,0.);
	minpt.setCoordinates(0.,0.,0.);
	xwidth = 0.;
	ywidth = 0.;
	zwidth = 0.;
	maxwidth = 0.;
	minwidth = 0.;

	maxdensity = 0.;
	mindensity = 0.;
	avedensity = 0.;


}


// ~Crowd
///03.11.23

Crowd::~Crowd()
{

}


//erase
///03.12.09

int Crowd::erase(int id)
{
	if(GpNode::erase(id) == -1)return -1;

	int n = (int)binlist.size();
	int er_index = id-1;
	int fn_index = n-1;
	int bin;
	for(int i=0;i<n;i++) {
		bin = binlist[i];
		if(bin == er_index) {
			deque<int>::iterator p = binlist.begin() + i;
			binlist.erase(p);
			n = (int)binlist.size();
			i--;
		} else if(bin == fn_index) {
			binlist[i] = er_index;
		} else if(bin > er_index) {
			binlist[i] -= 1;
		}
	}



	return getSize();
}


//clear
///03.12.09

void Crowd::clear()
{
	GpNode::clear();
	binlist.clear();
}


//makeNewNode
///03.12.09

Node* Crowd::makeNewNode(double x,double y,double z,double density)
{
	Node *pnd;

	pnd = GpNode::makeNewNode(x,y,z,density);

	binlist.push_back(pnd->getId()-1);

	return pnd;
}



// initstd
///03.11.24

void Crowd::initstd()
{
	double ix,iy,iz;
	Node *pnd;

	if(!getSize()) return;
	pnd = getNode(0);
	ix = pnd->getX();
	iy = pnd->getY();
	iz = pnd->getZ();
	maxpt.setCoordinates(ix,iy,iz);
	minpt.setCoordinates(ix,iy,iz);
	xwidth = ywidth = zwidth = 0.;
	maxwidth = minwidth = 0.;

	maxdensity = getNode(0)->getDensity();
	mindensity = avedensity = maxdensity;

	int n = getSize();
	for(int i=1;i<n;i++) {
		pnd = getNode(i);
		ix = pnd->getX();
		iy = pnd->getY();
		iz = pnd->getZ();
		if(maxpt.getX() < ix) {
			maxpt.setX(ix);
			maxaxispt[0] = getNode(i);
		} else if(minpt.getX() > ix) {
			minpt.setX(ix);
			minaxispt[0] = getNode(i);
		}
		if(maxpt.getY() < iy) {
			maxpt.setY(iy);
			maxaxispt[1] = getNode(i);
		} else if(minpt.getY() > iy) {
			minpt.setY(iy);
			minaxispt[1] = getNode(i);
		}
		if(maxpt.getZ() < iz) {
			maxpt.setZ(iz);
			maxaxispt[2] = getNode(i);
		} else if(minpt.getZ() > iz) {
			minpt.setZ(iz);
			minaxispt[2] = getNode(i);
		}
		if(maxdensity < getNode(i)->getDensity()) {
			maxdensity = getNode(i)->getDensity();
		} else if(mindensity > getNode(i)->getDensity()) {
			mindensity = getNode(i)->getDensity();
		}
		avedensity += getNode(i)->getDensity();
	}
	xwidth = maxpt.getX() - minpt.getX();
	ywidth = maxpt.getY() - minpt.getY();
	zwidth = maxpt.getZ() - minpt.getZ();
	maxwidth = minwidth = xwidth;
	if(maxwidth < ywidth) maxwidth = ywidth;
	else if(minwidth > ywidth) minwidth = ywidth;
	if(maxwidth < zwidth) maxwidth = zwidth;
	else if(minwidth > zwidth) minwidth = zwidth;

	avedensity = avedensity / getSize();

	flg_wid = true;

	return;
}

//
/*
bool Crowd::checkDensity()
{
	bool flg;
	int n = node.size();
	for(int i=0;i<n;i++) {
		flg = node[i]->hasDensity();
		if(!flg) return false;
	}
	return true;
}

*/
//

void Crowd::initDensity()
{
	Node *pnd;

	if(!getSize()) return;

	maxdensity = getNode(0)->getDensity();
	mindensity = avedensity = maxdensity;

	int n = getSize();
	for(int i=1;i<n;i++) {
		pnd = getNode(i);
		if(maxdensity < pnd->getDensity()) {
			maxdensity = pnd->getDensity();
		} else if(mindensity > pnd->getDensity()) {
			mindensity = pnd->getDensity();
		}
		avedensity += pnd->getDensity();
	}
	avedensity = avedensity / n;
	flg_den = true;

	return;
}


//normalize
///normalization of nodes
///03.11.25

void Crowd::normalize()
{
	if(flg_std) return;
	if(!flg_wid) initstd();

	stdrate = maxwidth;
	stdpt.setCoordinates(-minpt.getX(),-minpt.getY(),-minpt.getZ());

	Node *pnd;
	int n = getSize();
	for(int i=0;i<n;i++) {
		pnd = getNode(i);
		pnd->move(stdpt);
		pnd->reduce(stdrate);

	}

	flg_std = true;

	return;
}


//restore
///restorment of nodes
///03.11.23

void Crowd::restore()
{
	if(!flg_std) return;

	stdpt.setCoordinates(-stdpt.getX(),-stdpt.getY(),
		-stdpt.getZ());

	Node *pnd;
	int n = getSize();
	for(int i=0;i<n;i++) {
		pnd = getNode(i);
		pnd->expand(stdrate);
		pnd->move(stdpt);
	}
	flg_std = false;


	return;
}


// getMaxPoint
///03.11.23

Point Crowd::getMaxPoint()
{
	if(!flg_wid) initstd();
	return maxpt;
}


// getMinPoint
///03.11.23

Point Crowd::getMinPoint()
{
	if(!flg_wid) initstd();
	return minpt;
}


// getXwidth
///03.11.23

double Crowd::getXWidth()
{
	if(!flg_wid) initstd();
	return xwidth;
}


// getYWidth
///03.11.23

double Crowd::getYWidth()
{
	if(!flg_wid) initstd();
	return ywidth;
}


// getZWidth
///03.11.23

double Crowd::getZWidth()
{
	if(!flg_wid) initstd();
	return zwidth;
}


// getMaxwidth
///03.11.23

double Crowd::getMaxWidth()
{
	if(!flg_wid) initstd();
	return maxwidth;
}


// getMinWidth
///03.11.23

double Crowd::getMinWidth()
{
	if(!flg_wid) initstd();
	return minwidth;
}


//getMaxDensity
///03.11.23

double Crowd::getMaxDensity()
{
	if(!flg_wid) initstd();
	return maxdensity;
}


//getMinDensity
///03.11.23

double Crowd::getMinDensity()
{
	if(!flg_wid) initstd();
	return mindensity;
}


//getAveDensity
///03.11.23

double Crowd::getAveDensity()
{
	if(!flg_wid) initstd();
	return avedensity;
}


//getStdRate
///03.11.23

double Crowd::getStdRate()
{
	if(flg_std) return stdrate;
	else return 1.;
}


//getStdPoint

Point Crowd::getStdPoint()
{
	if(flg_std) return stdpt;
	else {
		Point tmp(0.,0.,0.);
		return tmp;
	}
}


//

Node *Crowd::getMaxXAxisNode()
{
	if(!flg_wid) initstd();

	return maxaxispt[0];
}

//

Node *Crowd::getMinXAxisNode()
{
	if(!flg_wid) initstd();

	return minaxispt[0];
}


//

Node *Crowd::getMaxYAxisNode()
{
	if(!flg_wid) initstd();

	return maxaxispt[1];
}


//

Node *Crowd::getMinYAxisNode()
{
	if(!flg_wid) initstd();

	return minaxispt[1];
}


//

Node *Crowd::getMaxZAxisNode()
{
	if(!flg_wid) initstd();

	return maxaxispt[2];
}


//

Node *Crowd::getMinZAxisNode()
{
	if(!flg_wid) initstd();

	return minaxispt[2];
}

//

Point Crowd::getCurrentMaxPoint()
{

	if(!flg_wid) initstd();
	if(flg_std) {
		Point tmp(maxpt);
		tmp.reduce(stdrate);
		return tmp;
	}
	return maxpt;
}

Point Crowd::getCurrentMinPoint()
{

	if(!flg_wid) initstd();
	if(flg_std) {
		Point tmp(minpt);
		tmp.reduce(stdrate);
		return tmp;
	}
	return minpt;
}

double Crowd::getCurrentMaxWidth()
{

	if(!flg_wid) initstd();
	if(flg_std) return maxwidth / stdrate;
	return maxwidth;
}

double Crowd::getCurrentMinWidth()
{

	if(!flg_wid) initstd();
	if(flg_std) return minwidth / stdrate;
	return minwidth;
}

double Crowd::getCurrentXWidth()
{
	if(!flg_wid) initstd();
	if(flg_std) return xwidth / stdrate;
	return xwidth;

}

double Crowd::getCurrentYWidth()
{

	if(!flg_wid) initstd();
	if(flg_std) return ywidth / stdrate;
	return ywidth;
}

double Crowd::getCurrentZWidth()
{

	if(!flg_wid) initstd();
	if(flg_std) return zwidth / stdrate;
	return zwidth;
}

double Crowd::getCurrentMaxDensity()
{

	if(!flg_wid) initstd();
	if(flg_std) return maxdensity / stdrate;
	return maxdensity;
}

double Crowd::getCurrentMinDensity()
{

	if(!flg_wid) initstd();
	if(flg_std) return mindensity / stdrate;
	return mindensity;
}

double Crowd::getCurrentAveDensity()
{

	if(!flg_wid) initstd();
	if(flg_std) return avedensity / stdrate;
	return avedensity;
}

//




// binsort
///03.11.25

void Crowd::binsort()
{
	int i,j,k,l;
	double ndiv,factx,facty,factz;
	deque<int> ibin;
	Node *nd;
	int num;

	/*////
	double stx,sty,stz;
	stx = stdpt.getX()/stdrate;
	sty = stdpt.getY()/stdrate;
	stz = stdpt.getZ()/stdrate;
	*//////

	///binlist.clear();
	num = getSize();
	///for(i=0;i<num;i++) binlist.push_back(i);

	ndiv = (int)pow((double)num,0.1);
	factx = ndiv/((getXWidth())*1.01/getMaxWidth());
	facty = ndiv/((getYWidth())*1.01/getMaxWidth());
	factz = ndiv/((getZWidth())*1.01/getMaxWidth());

	for(l=0;l<num;l++) {
		nd = getNode(l);
		i = (int)((nd->getX())*factx);
		j = (int)((nd->getY())*factz);
		k = (int)((nd->getZ())*facty);
		//////
		//i = (int)((nd->getX()-stx)*factx);
		//j = (int)((nd->getY()-sty)*factz);
		//k = (int)((nd->getZ()-stz)*facty);
		///////////
		if((k%2) == 0) {
			if((j%2) == 0) {
				ibin.push_back((int)(k*ndiv*ndiv + j*ndiv + i + 1));
			} else {
				ibin.push_back((int)(k*ndiv*ndiv + (j+1)*ndiv - i));
			}
		} else {
			if((j%2) == 0) {
				ibin.push_back((int)(k*ndiv*ndiv + (ndiv-j)*ndiv - i));
			} else {
				ibin.push_back((int)(k*ndiv*ndiv + (ndiv-j-1)*ndiv + i + 1));
			}
		}
	}


	dtm::quickSort(0,num-1,binlist,ibin);


	return;
}


int Crowd::getBinListNumber(int index)
{
	int n =  (int)binlist.size();
	if(index < 0 || index >= n) {
		Errors err(0,
			"int Crowd::getBinListNumber(int index)",
			"illegal index");
	}
	return binlist[index];
}


// laplasianTriangle
///03.11.23

void Crowd::laplasianTriangle(int first,int end)
{
	if(first < 0) first = 0;
	else if(first >= getSize()) return;

	if(end >= getSize()) end = getSize();
	else if(end < 0) return;

	for(int i=first;i<=end;i++) {
		if(getNode(i)->getType() != 2) continue;
		laplasInTri(getNode(i));
	}

	return;
}


// laplasianTetra
/*
void Crowd::laplasianTetra(int first,int end)
{
	int i;

	if(first > end) return;

	if(first < 0) first = 0;
	else if(first >= node.size()) return;
	
	if(end >= node.size()) end = node.size();
	else if(end < 0) return;

	for(i=first;i<=end;i++) {
		if(node[i]->getType() != 3) continue;
		node[i]->laplasInTet();
	}

	return;
}
*/

// setDensity
///03.11.23

void Crowd::setDensity(double idensity,int flg)
{

	if(flg) {
		if(flg_std) {
			idensity = idensity / stdrate;
		}
	}

	int n = getSize();
	for(int i=0;i<n;i++) {
		getNode(i)->setDensity(idensity);
	}

	return;
}


//clearHanger
///03.11.23

void Crowd::clearHanger(int fgType)
{
	int n = getSize();
	for(int i=0;i<n;i++) {
		getNode(i)->clearHanger(fgType);
	}

//	conected = 0;

	return;
}



///for debug///
// view
/*
void Crowd::view()
{
	cout << node.size() << "\n";

	int n = node.size();
	for(int i=0;i<n;i++) {
		node[i]->view();
	}

	return;
}
*/

// viewTriHanger
/**/
void Crowd::viewHanger()
{
	cout << getSize() << "\n";
	int n = getSize();
	for(int i=0;i<n;i++) {
		getNode(i)->viewHanger(0);
	}


	return;
}


}

//////////////////////////////////////////////////////////////////EOF