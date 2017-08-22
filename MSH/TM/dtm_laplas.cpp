/*

			dtm_laplas.h
			>laplasian of node
			
			last update : 2003.11.13			by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtm_laplas.h"


namespace dtm{


// laplasInTri

void laplasInTri(Node *nd)
{
	int i;
	double xc,yc,zc;
	Triangle *itri;
	Node *j1,*j2,*j3;
	double s;
	double cgrax,cgray,cgraz;
	double px,py,pz;
	double par,ar;

	double dens;
	double minarea;

	int num;

	if(nd->getType() != 2) return;
	num = nd->getHangerSize(FG_TYPE_TRI);
	if(!num) return;


	dens = nd->getDensity();
	if(dens < LIMIT_8) {
		minarea = LIMIT;
	} else { 
		minarea = dens*dens*0.05;
	}

	px = nd->getX();
	py = nd->getY();
	pz = nd->getZ();
	par = 0.0;
	cgrax = 0.;
	cgray = 0.;
	cgraz = 0.;
	for(i=0;i<num;i++) {
		itri = (Triangle*)nd->getHanger(i,FG_TYPE_TRI);
		j1 = itri->getNode(0);
		j2 = itri->getNode(1);
		j3 = itri->getNode(2);
		ar = itri->getArea();
		par += ar;
		if(ar < minarea) minarea = ar;
		xc = (j1->getX() + j2->getX() + j3->getX())/3; 
		yc = (j1->getY() + j2->getY() + j3->getY())/3; 
		zc = (j1->getZ() + j2->getZ() + j3->getZ())/3; 
		cgrax += xc*ar;
		cgray += yc*ar;
		cgraz += zc*ar;
	}

	cgrax = cgrax/par;
	cgray = cgray/par;
	cgraz = cgraz/par;

	nd->setCoordinates(cgrax,cgray,cgraz);

	ar = 0.0;
	for(i=0;i<num;i++) {		
		itri = (Triangle*)nd->getHanger(i,FG_TYPE_TRI);	
		s = itri->getArea();
		ar += s;
		if(s < minarea) {
			nd->setCoordinates(px,py,pz);
			return;
		
		}
	}

	if( fabs(par - ar) > LIMIT) {
		nd->setCoordinates(px,py,pz);
		return;
	}

	nd->renew();

	return;
}



void laplasInTriEdge(Node *nd)
{
	int i,j;
	double xc,yc,zc;
	Triangle *itri;
	Node *cp;
//	Node *j1,*j2,*j3;
	double s;
	double cgrax,cgray,cgraz;
	double px,py,pz;
	double par,ar;

	double ds,pds;

	double dens;
	double minarea;

	int num,n;

	if(nd->getType() != 1) return;
	num = nd->getHangerSize(FG_TYPE_TRI);
	if(!num) return;


	dens = nd->getDensity();
	if(dens < LIMIT_8) {
		minarea = LIMIT;
	} else { 
		minarea = dens*dens*0.05;
	}

	px = nd->getX();
	py = nd->getY();
	pz = nd->getZ();
	xc = 0.0;
	yc = 0.0;
	zc = 0.0;
	pds = 0.0;
	par = 0.0;
	n = 0;
	//prm = LIMIT_8;
	for(i=0;i<num;i++) {
		itri = (Triangle*)nd->getHanger(i,FG_TYPE_TRI);
		par += itri->getArea();
		for(j=0;j<3;j++) {
			cp = itri->getNode(j);
			if(cp == nd) {
				cp = itri->getNode((j+1)%3);
				if(cp->getType() == 0 || cp->getType() == 1) {
					ds = cp->getDensity();
					n++;
					pds += ds;
					xc += cp->getX()*ds;
					yc += cp->getY()*ds;
					zc += cp->getZ()*ds;
				break;
				}
			}
		}
	}

	if(n > 1) {
		cgrax = xc/pds;
		cgray = yc/pds;
		cgraz = zc/pds;
	} else return;

	nd->setCoordinates(cgrax,cgray,cgraz);

	ar = 0.0;
	for(i=0;i<num;i++) {		
		itri = (Triangle*)nd->getHanger(i,FG_TYPE_TRI);	
		s = itri->getArea();
		ar += s;
		if(s < minarea) {
			nd->setCoordinates(px,py,pz);
			return;
		
		}
	}

	if( fabs(par - ar) > LIMIT) {
		nd->setCoordinates(px,py,pz);
		return;
	}

	nd->renew();

	return;
}



// laplasInTet

void laplasInTet(Node *nd)
{
	int i;
	double gx,gy,gz;
	double xc,yc,zc;
	Tetra *itet;
	Node *j1,*j2,*j3,*j4;
	double s,ar,nar;
	double cgrax,cgray,cgraz;
	double px,py,pz;
	double prm;

	int num;



	if(nd->getType() != 3) return;
	num = nd->getHangerSize(FG_TYPE_TRI);
	if(!num) return;


	px = nd->getX();
	py = nd->getY();
	pz = nd->getZ();
	gx = 0.;
	gy = 0.;
	gz = 0.;
	ar = 0.;
	prm = LIMIT_8;


	for(i=0;i<num;i++) {
		itet = (Tetra*)nd->getHanger(i,FG_TYPE_TRI);
		j1 = itet->getNode(0);
		j2 = itet->getNode(1);
		j3 = itet->getNode(2);
		j4 = itet->getNode(3);

		s = itet->getVolume();

		xc = (j1->getX() + j2->getX() + j3->getX() + j4->getX())*0.25;
		yc = (j1->getY() + j2->getY() + j3->getY() + j4->getX())*0.25;
		zc = (j1->getZ() + j2->getZ() + j3->getZ() + j4->getX())*0.25;
		ar += s;
		gx += s*xc;
		gy += s*yc;
		gz += s*zc;
	}
	cgrax = gx/ar;
	cgray = gy/ar;
	cgraz = gz/ar;


	nd->setCoordinates(cgrax,cgray,cgraz);


	nar = 0.;
	for(i=0;i<num;i++) {		
		itet = (Tetra*)nd->getHanger(i,FG_TYPE_TRI);

		s = itet->getVolume();

		nar += s;
		if(s < prm) {
			nd->setCoordinates(px,py,pz);
			return;
		}
	}
	if(fabs(ar - nar) > LIMIT) {
		nd->setCoordinates(px,py,pz);
		return;
	}

	nd->renew();

	return;
}



}


//////////////////////////////////////////////////////////////////EOF