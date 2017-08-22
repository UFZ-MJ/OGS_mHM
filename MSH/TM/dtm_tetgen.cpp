/*

			dtm_tetgen.cpp
			>Tetgen class

			last update : 2003.12.08		by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtm_tetgen.h"

namespace dtm{



// Tetgen
///03.11.23

Tetgen::Tetgen()
{

	rate_minvolm = 0.05;
	rate_density = 0.8;

	for(int i=0;i<4;i++) sptet_nd[i] = NULL;

	tbl_long = 20;
	flg_tbl = 0;

	cr_crowd = NULL;
}


//Tetgen
///03.12.08

Tetgen::Tetgen(Crowd *crowd)
{
	rate_minvolm = 0.05;
	rate_density = 0.8;

	for(int i=0;i<4;i++) sptet_nd[i] = NULL;

	tbl_long = 20;
	flg_tbl = 0;

	cr_crowd = crowd;
}



// ~Tetgen
///03.11.23

Tetgen::~Tetgen()
{
	int i;
	for(i=0;i<4;i++) {
		if(sptet_nd[i]) delete sptet_nd[i];
	}
}


//setCurrentCrowd
///03.12.08

bool Tetgen::setCurrentCrowd(Crowd *crowd)
{
	cr_crowd = crowd;

	return true;
}



//initialize
///03.11.23

void Tetgen::setTetgenParameter(double minvol,double rtdens,int tbllg)
{
	rate_minvolm = minvol;
	rate_density = rtdens;
	tbl_long = tbllg;
//	flg_tbl = 0;
}





// setTetHanger
///03.11.23

void Tetgen::setTetHanger()
{
	Tetra *crtet;
	Node *nd;

	int n = getSize();
	for(int i=0;i<n;i++) {
		crtet = getTetra(i);
		for(int j=0;j<4;j++) {
			nd = crtet->getNode(j);
			nd->insertHanger(crtet,FG_TYPE_TET);
		}
	}


	return;
}



//delaunay
///03.12.08

void Tetgen::delaunay(Crowd *crowd)
{
	setCurrentCrowd(crowd);
	delaunay();
}



// delauny
///03.11.23

void Tetgen::delaunay()
{
	Node *nd;
	int index;

	if(!cr_crowd) return;

	clear();

	//bin sorting of each nodes
	cr_crowd->binsort();

	if(cr_crowd->isStd()) {
		//large tetrahedron containing all points
		makeSuperTetra();
	} else {
		makeSuperTetra(cr_crowd->getMaxWidth(),cr_crowd->getMinPoint());
	}

	//mesh useing surface node
	int n = cr_crowd->getSize();
	for(int i=0;i<n;i++) {
		index = cr_crowd->getBinListNumber(i);
		nd = cr_crowd->getNode(index);
		//mesh generator
		generator(nd);
	}
	//removing supertetrahedra
	removeSuperTetra();

	return;
}




// generator
///03.11.23

void Tetgen::generator(Node *pnd)
{
	Tetra *loc;
	vector<int>tet_id_buff;
	vector<Adpt>adpt_buff;	//adapter of polygon(!= class Polygon) to tetrahedora


	loc = locate(pnd);

	getPolyTetra(tet_id_buff,loc,pnd);
	correctPolyTetraAdd(tet_id_buff,adpt_buff,pnd);
	reMeshTetra(tet_id_buff,adpt_buff,pnd);

	conectNewTetra(tet_id_buff,(int)adpt_buff.size());
	int f = (int)adpt_buff.size();
	int e = (int)tet_id_buff.size();
	for(int i=f;i<e;i++) erase(tet_id_buff[i]);



	return;
}


//addNode
///03.12.08

int Tetgen::addNode(vector<int>&tet_id_buff,Tetra *loc,Node *add_node)
{
	vector<Adpt>adpt_buff;

	getPolyTetra(tet_id_buff,loc,add_node);
	correctPolyTetraCut(tet_id_buff,adpt_buff,add_node);
	reMeshTetra(tet_id_buff,adpt_buff,add_node);
	conectNewTetra(tet_id_buff,(int)adpt_buff.size());
	int f = (int)adpt_buff.size();
	int e = (int)tet_id_buff.size();
	for(int i=f;i<e;i++) erase(tet_id_buff[i]);

	return (int)tet_id_buff.size();
}



// locate
///03.12.09
//search tetra including point
//this function used only simple(non hollow) object
//faster than optlocate

Tetra* Tetgen::locate(Point *ppt)
{
	Tetra *crtet;
	Node *p0,*p1,*p2;
	Plane plane;
	double ans;
	const int ij0[4] = {1,2,3,0};
	const int ij1[4] = {3,3,1,1};
	const int ij2[4] = {2,0,0,2};
	int counter = 0;


	int i = 0;
	crtet = getTetra(getSize()-1);
	while (i < 4) {
		p0 = crtet->getNode(ij0[i]);
		p1 = crtet->getNode(ij1[i]);
		p2 = crtet->getNode(ij2[i]);
		plane.setVector(p0,p1,p2); 
		ans = plane.getAnswer(ppt);
		if (ans < -LIMIT_10) {
			counter++;
			if(counter > LIMIT_NUMBER_OF_ELEMENTS) {
				Errors err(0,
					"Tetra* Tetgen::locate(Point *ppt)",
					"cannot find tetra");
				throw err;
			}
			crtet = crtet->getNeighborElement(i);
			if(!crtet) {
				Errors err(0,
					"Tetra* Tetgen::locate(Point *ppt)",
					"cannot find tetra");
				throw err;
			}
			i = 0;
			continue;
		}
		i++;
	} //while i<4


	return crtet;
}


// optlocate
///03.11.23
//search tetra including point
//this function can be used any type object
//sometime very slow

Tetra* Tetgen::optlocate(Point *ppt)
{
	Tetra *crtet = NULL;
	Node *p0,*p1,*p2;
	Plane plane;
	double ans;
	int flg;
	const int ij0[4] = {1,2,3,0};
	const int ij1[4] = {3,3,1,1};
	const int ij2[4] = {2,0,0,2};

	int n = getSize();
	for(int i=0;i<n;i++) {
		crtet = getTetra(i);
		flg = 1;
		for(int j=0;j<4;j++) {
			p0 = crtet->getNode(ij0[j]);
			p1 = crtet->getNode(ij1[j]);
			p2 = crtet->getNode(ij2[j]);

			plane.setVector(p0,p1,p2);
			ans = plane.getAnswer(ppt);

			if (ans < -LIMIT) {
				flg = 0;
				break;
			}
		}
		if(flg) break;
	}

	return crtet;
}


// getPolyTetra
///03.11.23
//get tetra for remesh

int Tetgen::getPolyTetra(vector<int>&tet_id_buff,Tetra *loc,Node *cr_nd)
{
//	stack<Tetra*, vector<Tetra*> >istack;
	vector<Tetra*>istack;	//stack of tetra
	Tetra *crtet,*nbtet;
	double radius,dist;
	Point* center;
	int domain;


	loc->setState(1);
	tet_id_buff.push_back(loc->getId());
	istack.push_back(loc);
	//////////
	domain = loc->getDomain();
	///////////

	while (istack.size()) {
		crtet = istack[istack.size()-1];
		istack.pop_back();
		//crtet = istack[0];
		//istack.pop_front();
		for (int i=0;i<4;i++) {
			nbtet = crtet->getNeighborElement(i);
			if (!nbtet) continue;
			if (nbtet->getState() == 1) continue;
			///////////
			if(nbtet->getDomain() != domain) continue;
			///////////
			radius = nbtet->getRadius();
			center = nbtet->getCenter();
			dist = getDistance(center,cr_nd);
			if(dist > radius) continue;
			nbtet->setState(1);
			tet_id_buff.push_back(nbtet->getId());
			istack.push_back(nbtet);
		} // for i
	} //while istack.size()


	return (int)tet_id_buff.size();
}


//correctPolyTetraAdd
///03.12.08

int Tetgen::correctPolyTetraAdd(vector<int>&tet_id_buff,
								vector<Adpt>&adpt_buff,Node *cr_nd)
{
	int i,j;
	Tetra *crtet,*nbtet;
	Node *nd0,*nd1,*nd2;
	Adpt tmp_adpt;
	const int ij0[4] = {1,2,3,0};
	const int ij1[4] = {3,3,1,1};
	const int ij2[4] = {2,0,0,2};

	//pick up surface triangle(=adpt_buff) of poly tetrahedra(polygon)(=tet_id_buff)
	i = 0;
	int n = (int)tet_id_buff.size();
	while (i < n) {
		crtet = getTetra(tet_id_buff[i]-1);
		for (j=0;j<4;j++) {
			nbtet = crtet->getNeighborElement(j);
			nd0 = crtet->getNode(ij0[j]);
			nd1 = crtet->getNode(ij1[j]);
			nd2 = crtet->getNode(ij2[j]);
			if (!nbtet) {
				tmp_adpt.face[0] = nd0;
				tmp_adpt.face[1] = nd1;
				tmp_adpt.face[2] = nd2;
				tmp_adpt.neib = NULL;
				tmp_adpt.num = 0;
				tmp_adpt.vol = getTetraVolume(nd0,nd1,nd2,cr_nd);
				adpt_buff.push_back(tmp_adpt);
			} //if !nbtet
			else if (nbtet->getState() != 1) {
				tmp_adpt.face[0] = nd0;
				tmp_adpt.face[1] = nd1;
				tmp_adpt.face[2] = nd2;
				tmp_adpt.neib = nbtet;
				tmp_adpt.num = nbtet->getConection(crtet);
				tmp_adpt.vol = getTetraVolume(nd0,nd1,nd2,cr_nd);
				adpt_buff.push_back(tmp_adpt);
				if (tmp_adpt.vol < LIMIT) {
					tet_id_buff.push_back(nbtet->getId());
					nbtet->setState(1);
					adpt_buff.clear();
					i = -1;
					n = (int)tet_id_buff.size();
					break;
				}
			} //else if nbtet->getState() != 1
		} //for j
		i++;
	} //while i < tet_id_buff.size()


	return (int)tet_id_buff.size();
}


//correctPolyTetraCut
///03.12.08

int Tetgen::correctPolyTetraCut(vector<int>&tet_id_buff,
								vector<Adpt>&adpt_buff,Node *cr_nd)
{

	int i,j;
	Tetra *crtet,*nbtet;
	Adpt tmp_adpt;
	const int ij0[4] = {1,2,3,0};
	const int ij1[4] = {3,3,1,1};
	const int ij2[4] = {2,0,0,2};
	///adpt.reserve(tet_id_buff.size()*2);

	//pick up surface triangle(=adpt_buff) of poly tetrahedra(polygon)(=tet_id_buff)
	//this part is diffarent from function Poly (other parts are almost same)
	i = 0;
	int n = (int)tet_id_buff.size();
	while (i < n) {
		crtet = getTetra(tet_id_buff[i]-1);
		for (j=0;j<4;j++) {
			nbtet = crtet->getNeighborElement(j);
			if (!nbtet) {
				tmp_adpt.face[0] = crtet->getNode(ij0[j]);
				tmp_adpt.face[1] = crtet->getNode(ij1[j]);
				tmp_adpt.face[2] = crtet->getNode(ij2[j]);
				tmp_adpt.neib = NULL;
				tmp_adpt.num = 0;
				tmp_adpt.vol = getTetraVolume(tmp_adpt.face[0],tmp_adpt.face[1],tmp_adpt.face[2],cr_nd);
				adpt_buff.push_back(tmp_adpt);
				if(tmp_adpt.vol < LIMIT) {
					crtet->setState(-1);
					tet_id_buff[i] = tet_id_buff[tet_id_buff.size()-1];
					tet_id_buff.pop_back();
					adpt_buff.clear();
					i = -1;
					n = (int)tet_id_buff.size();
					break;
				}
			} //if !nbtet
			else if (nbtet->getState() != 1) {
				tmp_adpt.face[0] = crtet->getNode(ij0[j]);
				tmp_adpt.face[1] = crtet->getNode(ij1[j]);
				tmp_adpt.face[2] = crtet->getNode(ij2[j]);
				tmp_adpt.neib = nbtet;
				tmp_adpt.num = nbtet->getConection(crtet);
				tmp_adpt.vol = getTetraVolume(tmp_adpt.face[0],tmp_adpt.face[1],tmp_adpt.face[2],cr_nd);
				adpt_buff.push_back(tmp_adpt);
				if (tmp_adpt.vol < LIMIT) {
					crtet->setState(-1);
					tet_id_buff[i] = tet_id_buff[tet_id_buff.size()-1];
					tet_id_buff.pop_back();
					adpt_buff.clear();
					n = (int)tet_id_buff.size();
					i = -1;
					break;
				}
			} //else if nbtet->getState() != 1
		} //for j
		i++;
	} //while i < tet_id_buff.size()

	return (int)tet_id_buff.size();
}


//reMeshTetra
///03.12.08

int Tetgen::reMeshTetra(vector<int>&tet_id_buff,vector<Adpt>&adpt_buff,Node *cr_nd)
{
	int i;
	Tetra *crtet;
	Node *nd0,*nd1,*nd2;
	int domain = 0;

	//refresh tetrahedra
	if(tet_id_buff.size() <= adpt_buff.size()) {
		//recycle old tetrahedra
		int np = (int)tet_id_buff.size();
		for(i=0;i<np;i++) {
			crtet = getTetra(tet_id_buff[i]-1);
			domain = crtet->getDomain();
			nd0 = adpt_buff[i].face[0];
			nd1 = adpt_buff[i].face[1];
			nd2 = adpt_buff[i].face[2];
			crtet->resetNode(nd0,nd1,nd2,cr_nd);
			crtet->setVolume(adpt_buff[i].vol);
			crtet->resetNeighborElement(3,adpt_buff[i].neib);
			if(adpt_buff[i].neib) {
				adpt_buff[i].neib->resetNeighborElement(adpt_buff[i].num,crtet);
			}
			//crtet->setSphere();
			crtet->setState(0);
			crtet->setDomain(domain);
		} //for i
		//if run out old one, make new tetrahedra
		int na = (int)adpt_buff.size();
		for(i=np;i<na;i++) {
			int ip = getSize()+1;
			nd0 = adpt_buff[i].face[0];
			nd1 = adpt_buff[i].face[1];
			nd2 = adpt_buff[i].face[2];
			crtet = makeNewTetra(nd0,nd1,nd2,cr_nd);
			crtet->setVolume(adpt_buff[i].vol);
			crtet->resetNeighborElement(3,adpt_buff[i].neib);
			if(adpt_buff[i].neib) {
				adpt_buff[i].neib->resetNeighborElement(adpt_buff[i].num,crtet);
			}
			//crtet->setSphere();
			crtet->setState(0);
			crtet->setDomain(domain);
			tet_id_buff.push_back(ip);
		} //for i
	} else {
		//recycle old tetrahedra
		int na = (int)adpt_buff.size();
		for(i=0;i<na;i++) {
			crtet = getTetra(tet_id_buff[i]-1);
			domain = crtet->getDomain();
			nd0 = adpt_buff[i].face[0];
			nd1 = adpt_buff[i].face[1];
			nd2 = adpt_buff[i].face[2];
			crtet->resetNode(nd0,nd1,nd2,cr_nd);
			crtet->setVolume(adpt_buff[i].vol);
			crtet->resetNeighborElement(3,adpt_buff[i].neib);
			if(adpt_buff[i].neib) {
				adpt_buff[i].neib->resetNeighborElement(adpt_buff[i].num,crtet);
			}
			//crtet->setSphere();
			crtet->setState(0);
			crtet->setDomain(domain);
		}
		sort(tet_id_buff.begin()+i,tet_id_buff.end(),greater<int>());
	} //if tet_id_buff.size() <= adpt_buff.size() , else

	return (int)tet_id_buff.size();
}


//conectNewTetra
///03.12.08

int Tetgen::conectNewTetra(vector<int>&tet_id_buff,int nend)
{
	int i,j,k,l;
	int flg;
	Tetra *crtet;
	Node *cr_nd1,*cr_nd2,*nb_nd1,*nb_nd2;
	Adpt tmp_adpt;
	const int jj0[3] = {1,2,0};
	const int jj1[3] = {2,0,1};

	if(flg_tbl) {
		for(i=0;i<tbl_long;i++) tbl[i].clear();
	} else {
		tbl.resize(tbl_long);
		flg_tbl = 1;
	}

	//conect each new tetrahedra
	//int np = tet_id_buff.size();
	for(i=0;i<nend;i++) {
		crtet = getTetra(tet_id_buff[i]-1);
		for(j=0;j<3;j++) {
			flg = 0;
			cr_nd1 = crtet->getNode(jj0[j]);
			cr_nd2 = crtet->getNode(jj1[j]);
			////////+10:for surertetra`s nodes
			l = (cr_nd1->getId() + cr_nd2->getId() + 10)%tbl_long;
			///////////
			int n = (int)tbl[l].size();
			for(k=0;k<n;k++) {
				tmp_adpt = tbl[l][k];
				nb_nd1 = tmp_adpt.face[0];
				nb_nd2 = tmp_adpt.face[1];
				if(cr_nd1 == nb_nd1 && cr_nd2 == nb_nd2) {
					crtet->resetNeighborElement(j,tmp_adpt.neib);
					tmp_adpt.neib->resetNeighborElement(tmp_adpt.num,crtet);
					tbl[l][k] = tbl[l][n-1];
					tbl[l].pop_back();
					flg = 1;
					break;
				}
			} //for k
			if(flg) continue;
			tmp_adpt.face[0] = cr_nd2;
			tmp_adpt.face[1] = cr_nd1;
			tmp_adpt.neib = crtet;
			tmp_adpt.num = j;
			tbl[l].push_back(tmp_adpt);

		} //for j
	} //for i

	return (int)tet_id_buff.size();
}


// makeSuperTetra
///03.11.23

void Tetgen::makeSuperTetra(double maxwidth,Point minpt)
{
	double ix,iy,iz;

	sptet_nd[0] = new Node(-1);
	sptet_nd[1] = new Node(-2);
	sptet_nd[2] = new Node(-3);
	sptet_nd[3] = new Node(-4);

	ix =   86.60*maxwidth+minpt.getX();
    iy =    0.00+minpt.getY();
    iz =  -50.00*maxwidth+minpt.getZ();
	sptet_nd[0]->setCoordinates(ix,iy,iz);
    ix =  -43.30*maxwidth+minpt.getX();
    iy =   75.00*maxwidth+minpt.getY();
    iz =  -50.00*maxwidth+minpt.getZ();
	sptet_nd[1]->setCoordinates(ix,iy,iz);
   ix =  -43.30*maxwidth+minpt.getX();
    iy =  -75.00*maxwidth+minpt.getY();
    iz =  -50.00*maxwidth+minpt.getZ();
	sptet_nd[2]->setCoordinates(ix,iy,iz);
    ix =    0.00+minpt.getX();
    iy =    0.00+minpt.getY();
    iz =  100.00*maxwidth+minpt.getZ();
	sptet_nd[3]->setCoordinates(ix,iy,iz);

	makeNewTetra(sptet_nd[0],sptet_nd[1],sptet_nd[2],sptet_nd[3]);
	getTetra(0)->setVolume();


	return;
}


// makeSuperTetra
///03.11.23

void Tetgen::makeSuperTetra()
{
	double ix,iy,iz;

	sptet_nd[0] = new Node(-1);
	sptet_nd[1] = new Node(-2);
	sptet_nd[2] = new Node(-3);
	sptet_nd[3] = new Node(-4);

	ix =   86.60;
    iy =    0.00;
    iz =  -50.00;
	sptet_nd[0]->setCoordinates(ix,iy,iz);
    ix =  -43.30;
    iy =   75.00;
    iz =  -50.00;
	sptet_nd[1]->setCoordinates(ix,iy,iz);
    ix =  -43.30;
    iy =  -75.00;
    iz =  -50.00;
	sptet_nd[2]->setCoordinates(ix,iy,iz);
   ix =    0.00;
    iy =    0.00;
    iz =  100.00;
	sptet_nd[3]->setCoordinates(ix,iy,iz);

	makeNewTetra(sptet_nd[0],sptet_nd[1],sptet_nd[2],sptet_nd[3]);
	getTetra(0)->setVolume();



	return;
}


// removeSuperTetra
///03.11.23

void Tetgen::removeSuperTetra()
{
	Node *pnd;

	int n = getSize();
	for(int i=0;i<n;i++) {
		getTetra(i)->setState(0);
		for(int j=0;j<4;j++) {
			pnd = getTetra(i)->getNode(j);
			if(pnd->getId() < 0) {
				getTetra(i)->setState(-1);
				break;
			}
		}
	}
	removeTetrahedra();

	for(int i=0;i<4;i++) {
		delete sptet_nd[i];
		sptet_nd[i] = NULL;
	}

	return;
}



// fineTetra
///03.11.23
//subsequent refinement of tetrahedra

void Tetgen::fineTetra()
{
	int i;
	int n;
	Tetra *crtet;
	Node *nd0,*nd1,*nd2,*nd3;
	Node *new_nd;
	Point tmp_pt;
	double vol;
	double dist;
	double mindens,minvol;
	double nx,ny,nz,nden;
	int flg;
	vector<int> tet_id_buff;


	if(!cr_crowd) return;

	mindens = cr_crowd->getCurrentMinDensity();
	minvol = mindens*mindens*mindens*rate_minvolm;


	flg = 1;
	while(flg){
		flg = 0;
		n = getSize();
		for(i=0;i<n;i++) {
			crtet = getTetra(i);
			vol = crtet->getVolume();
			if(vol < minvol) continue;
			nd0 = crtet->getNode(0);
			nd1 = crtet->getNode(1);
			nd2 = crtet->getNode(2);
			nd3 = crtet->getNode(3);
			nx = (nd0->getX() + nd1->getX() + nd2->getX() + nd3->getX())*0.25;
			ny = (nd0->getY() + nd1->getY() + nd2->getY() + nd3->getY())*0.25;
			nz = (nd0->getZ() + nd1->getZ() + nd2->getZ() + nd3->getZ())*0.25;
			tmp_pt.setCoordinates(nx,ny,nz);

			dist = getDistance(&tmp_pt,nd0);
			if(nd0->getDensity() > dist*rate_density) continue;
			dist = getDistance(&tmp_pt,nd1);
			if(nd1->getDensity() > dist*rate_density) continue;
			dist = getDistance(&tmp_pt,nd2);
			if(nd2->getDensity() > dist*rate_density) continue;
			dist = getDistance(&tmp_pt,nd3);
			if(nd3->getDensity() > dist*rate_density) continue;
			
			flg = 1;
			nden = (nd0->getDensity() + nd1->getDensity() + 
				nd2->getDensity() + nd3->getDensity())*0.25;

			new_nd = cr_crowd->makeNewNode(nx,ny,nz,nden);
			new_nd->setType(3);

			//remeshing
			tet_id_buff.clear();
			addNode(tet_id_buff,crtet,new_nd);

			if(!tet_id_buff.size()) {
				cr_crowd->erase(new_nd->getId());
				continue;
			}
			n = getSize();
		} //for i
	}  //while flg

	/*
	view();
	getch();
	setTetHanger();
	

		Node *nnd = cr_crowd->getNode(202);
	
		nnd->view();
		getch();
		eraseNode(nnd);

		nnd = cr_crowd->getNode(203);
	
		nnd->view();
		getch();
		eraseNode(nnd);
		getch();

	view();
	getch();
	*/

	return;
}


//removeTetrahedra
///03.12.08

void Tetgen::removeTetrahedra()
{
	int n = getSize();
	for(int i=0;i<n;i++) {
		if(getTetra(i)->getState() == -1) {
			n = erase(i+1);
			i--;
			continue;
		}
		getTetra(i)->setState(0);
	}

	return;
}


//setDomain
///03.12.08

void Tetgen::setDomain(int dm)
{
	int n = getSize();
	for(int i=0;i<n;i++) {
		if(getTetra(i)->getState() == 1) {
			getTetra(i)->setDomain(dm);
		}
		getTetra(i)->setState(0);
	}

	return;
}


//cleanBadTetra
///03.12.12

void Tetgen::cleanBadTetra(double min_value)
{
	Node *nd;
	Tetra *itet;
	double mindens,minvolm;


	mindens = cr_crowd->getCurrentMinDensity();
	minvolm = mindens*mindens*mindens*rate_minvolm;
	
	int n = getSize();
	for(int i=0;i<n;i++) {
		itet = getTetra(i);

		if(itet->getVolume() > minvolm &&
			itet->getValue() > min_value) continue;
		
		for(int j=0;j<4;j++) {
			nd = itet->getNode(j);
			if(nd->getType() == 0 )continue;
			if(eraseNode(nd))break;
		}
		n = getSize();
	}



	return;
}


//eraseNode
///03.12.11

Node* Tetgen::eraseNode(Node *er_nd)
{

	int i,j,k;
	int n,nn;
	Tetra *itet;
	Node *pnd,*tmpnd;
	double lg;
	deque<double>lg_buff;
	deque<Node*>po_nd_buff;
	deque<double>::iterator it_d;
	deque<Node*>::iterator it_n;
	int flg;
	int type;
	int tp_flg = 0;
	int size = 0;

	////
	printf("erase  ");
	er_nd->view();

	type = er_nd->getType();
	if(type == 3) {
		tp_flg = 0;
	} else if(type == 2) {
		tp_flg = 2;
	} else if(type == 1) {
		tp_flg = 1;
	} else {
		return NULL;
	}

	printf("%5d\n",tp_flg);


	n = er_nd->getHangerSize(FG_TYPE_TET);
	for(i=0;i<n;i++) {
		itet = (Tetra*)er_nd->getHanger(i,FG_TYPE_TET);
		for(j=0;j<4;j++) {
			pnd = itet->getNode(j);
			if(*pnd == *er_nd)continue;
			if(tp_flg == 1) {
				if(pnd->getType() == 2)continue;
				if(pnd->getType() == 3)continue;
			} else if(tp_flg == 2) {
				if(pnd->getType() == 3)continue;
			}

			flg = 0;
			nn = (int)po_nd_buff.size();
			for(k=0;k<nn;k++) {
				tmpnd = po_nd_buff[k];
				if(*tmpnd == *pnd) {
					flg = 1;
					break;
				}
			}
			if(flg)continue;
			lg = getDistanceSquare(pnd,er_nd);
			for(k=0;k<nn;k++) {
				if(lg < lg_buff[k])break;
			}
			if(k == 0) {
				po_nd_buff.push_front(pnd);
				lg_buff.push_front(lg);
			} else if(k == nn){
				po_nd_buff.push_back(pnd);
				lg_buff.push_back(lg);
			} else {
				it_n = po_nd_buff.begin() + k;
				po_nd_buff.insert(it_n,pnd);
				it_d = lg_buff.begin() + k;
				lg_buff.insert(it_d,lg);
			}
		}
	}

	printf("buf\n");
    size = (int)po_nd_buff.size();
	for(i=0;i<size;i++) {
		po_nd_buff[i]->view();
		printf("%f\n",lg_buff[i]);
	}/**/
	getch();

	//try merge in order of length
	n = (int)po_nd_buff.size();
	for(i=0;i<n;i++) {
		pnd = po_nd_buff[i];

		pnd->view();
		getch();

		if(!mergeNode(pnd,er_nd)) continue;
		return pnd;
	}

	return NULL;
}


//mergeNode
///03.12.11

Node* Tetgen::mergeNode(Node *po_nd,Node *ne_nd)
{
	int i,j;
	Tetra *potet,*netet,*nbtet;
	deque<Tetra*>pohang;
	deque<Tetra*>nehang;
	int edge,ndedge;
	int flg;
	int size = 0;
	int size2 = 0;

	//set pohang ang nehang
	//tetra of pohang state = 1
	//tetra of nehang state = -1
	divHanger(po_nd,ne_nd,pohang,nehang);

	//////
	printf("po ");
	po_nd->view();
	printf("ne ");
	ne_nd->view();
	printf("pohang\n");
    size = (int)pohang.size();
	for(i=0;i<size;i++) {
		pohang[i]->view();
	}
	printf("nehang\n");
    size2 = (int)nehang.size();
	for(i=0;i< size2;i++) {
		nehang[i]->view();
	}
//	getch();
/*	*/

	//check vol, if can merge or not
	if(!isPossibleMerge(po_nd,ne_nd,pohang,nehang)) {
		int nn = (int)pohang.size();
		for(i=0;i<nn;i++) {
			pohang[i]->setState(0);
		}
		nn = (int)nehang.size();
		for(i=0;i<nn;i++) {
			nehang[i]->setState(0);
		}
		return NULL;
	}


	int n = (int)pohang.size();
	for(i=0;i<n;i++) {
		potet = pohang[i];
		////
		printf("i %d\n",i);
		potet->view();

		flg = 0;
		for(j=0;j<4;j++) {
			netet = potet->getNeighborElement(j);
			if(!netet) continue;
			if(netet->getState() == -1) {
				netet->view();

				//conect nw_nb to potet
				edge = netet->getNodeNumber(ne_nd);
				nbtet = netet->getNeighborElement(edge);
				potet->resetNeighborElement(j,nbtet);

				potet->view();
				
				if(nbtet) {
					nbtet->view();
	
					edge = nbtet->getConection(netet);
					nbtet->resetNeighborElement(edge,potet);

					nbtet->view();
				}
			}
	
			getch();
		}
		//change node to po_nd
		ndedge = potet->getNodeNumber(ne_nd);
		potet->setNode(ndedge,po_nd);
		///potet->setState(0);

		///
	//	potet->view();
	//	getch();
	}
	n = (int)nehang.size();
	for(i=0;i<n;i++) {
		netet = nehang[i];
		erase(netet->getId());
	}
	cr_crowd->erase(ne_nd->getId());

	n = (int)pohang.size();
	for(i=0;i<n;i++) {
		potet = pohang[i];
		potet->setState(0);
	}


	//printf("END MERGE\n");


	return po_nd;
}


//divHanger
///03.12.10
///divide hanger to remain tetrahedras and erase tetrahedras

void Tetgen::divHanger(Node *po_nd,Node *ne_nd,deque<Tetra*>&pohang,
						deque<Tetra*>&nehang)
{
	int flg;
	Node *pnd;
	Tetra *hgtet;

	int n = ne_nd->getHangerSize(FG_TYPE_TET);
	for(int i=0;i<n;i++) {
		hgtet = (Tetra*)ne_nd->getHanger(i,FG_TYPE_TET);
		flg = 0;
		for(int j=0;j<4;j++) {
			pnd = hgtet->getNode(j);
			if(*pnd == *po_nd) {
				flg = 1;
				break;
			}
		}
		if(flg) {
			hgtet->setState(-1);
			nehang.push_back(hgtet);	//erase tetrahedras
		} else {
			hgtet->setState(1);
			pohang.push_back(hgtet);	//remain tetrahedras
		}
	}

	return;
}


//isPossibleMerge
///03.12.10

int Tetgen::isPossibleMerge(Node *po_nd,Node *ne_nd,deque<Tetra*>&pohang,
							deque<Tetra*>&nehang)
{
	double vol,pre_vol,pst_vol;
	Tetra *ptet;
	Node *nd1,*nd2,*nd3;
	int edge;
	double mindens,minvolm;
	const int ij0[4] = {1,2,3,0};
	const int ij1[4] = {3,3,1,1};
	const int ij2[4] = {2,0,0,2};

	mindens = cr_crowd->getCurrentMinDensity();
	minvolm = mindens*mindens*mindens*rate_minvolm*0.5;

	pre_vol = pst_vol = 0.;
	int n = (int)pohang.size();
	for(int i=0;i<n;i++) {
		ptet = pohang[i];
		pre_vol += ptet->getVolume();
		edge = ptet->getNodeNumber(ne_nd);
		nd1 = ptet->getNode(ij0[edge]);
		nd2 = ptet->getNode(ij1[edge]);
		nd3 = ptet->getNode(ij2[edge]);
		vol = getTetraVolume(nd1,nd2,nd3,po_nd);

		printf("pr %20.15f\n",pre_vol);
		printf("%5d%5d%5d%5d\n",po_nd->getId(),nd1->getId(),nd2->getId(),
			nd3->getId());
		printf("%20.15f\n",vol);
		//check area
		//if new area is too small( area < minvolm ), do not merge
		//but it possible to merge, if you want to do
		
		////if(vol < minvolm) return 0;
		if(vol < LIMIT_10) return 0;
		pst_vol += vol;
	}
	n = (int)nehang.size();
	for(int i=0;i<n;i++) {
		ptet = nehang[i];
		pre_vol += ptet->getVolume();
	}

	printf("pre %20.15f\n",pre_vol);
	printf("pst %20.15f\n",pst_vol);
	getch();

	//if pst_vol > pre_vol(or pst_vol < pre_vol ), cannot merge
	if(fabs(pre_vol - pst_vol) > LIMIT) return 0;


	return 1;
}




///for debug///
// view

void Tetgen::view()
{

	cout << "TETGEN SIZE = " << getSize() << "\n";
	int n = getSize();
	for(int i=0;i<n;i++) {
		getTetra(i)->view();
	}

	return;
}




}


//////////////////////////////////////////////////////////////////EOF