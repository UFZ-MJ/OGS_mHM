/*

			dtm_stdio.h
			>input and output

			last update : 2003.11.01		by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtm_stdio.h"

namespace dtm{

// inpDtm3dPoint

bool inpDtm3dPoint(char *fname,Crowd *crowd)
{
	int i;
	int tmpn,tmpid;
	double tmpx,tmpy,tmpz;
	
	ifstream in(fname);
	if(!in) {
		return false;
	}

	in >> tmpn;
	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tmpx;
		in >> tmpy;
		in >> tmpz;
		crowd->makeNewNode(tmpx,tmpy,tmpz,0.);
	}

	crowd->initstd();

	in.close();

	return true;
}


// inpDtm3dTetra

bool inpDtm3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd)
{
	int i;
	int tmpn,tmpid;
	double tmpx,tmpy,tmpz;
	int tt0,tt1,tt2,tt3;
	Node *n0,*n1,*n2,*n3;
	Tetra *tet;

	ifstream in(fname);
	if(!in) {
		return false;
	}

	in >> tmpn;

	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tmpx;
		in >> tmpy;
		in >> tmpz;
		crowd->makeNewNode(tmpx,tmpy,tmpz,0.);
	}

	in >> tmpn;
	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tt0;
		in >> tt1;
		in >> tt2;
		in >> tt3;
		n0 = crowd->getNode(tt0-1);
		n1 = crowd->getNode(tt1-1);
		n2 = crowd->getNode(tt2-1);
		n3 = crowd->getNode(tt3-1);
		tet = tetgen->makeNewTetra(n0,n1,n2,n3);
		tet->setVolume();
	}

	crowd->initstd();

	in.close();

	return true;
}


// inpDtm3dTriangle

bool inpDtm3dTriangle(char *fname,Surface *surface,Crowd *crowd)
{
	int i;
	int tmpn,tmpid;
	double tmpx,tmpy,tmpz;
	int tt0,tt1,tt2;
	Node *n0,*n1,*n2;
	Triangle *tri;

	ifstream in(fname);
	if(!in) {
		return false;
	}

	in >> tmpn;

	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tmpx;
		in >> tmpy;
		in >> tmpz;
		crowd->makeNewNode(tmpx,tmpy,tmpz,0.);
	}

	in >> tmpn;
	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tt0;
		in >> tt1;
		in >> tt2;
		n0 = crowd->getNode(tt0-1);
		n1 = crowd->getNode(tt1-1);
		n2 = crowd->getNode(tt2-1);
		tri = surface->makeNewTriangle(n0,n1,n2);
	}

	crowd->initstd();

	in.close();

	return true;
}


// inpDtm3dTriangleWithDensity

bool inpDtm3dTriangleWithDensity(char *fname,Surface *surface,Crowd *crowd)
{
	int i;
	int tmpn,tmpid;
	double tmpx,tmpy,tmpz;
	double tmpd;
	int tt0,tt1,tt2;
	Node *n0,*n1,*n2;
	Triangle *tri;

	ifstream in(fname);
	if(!in) {
		return false;
	}

	in >> tmpn;

	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tmpx;
		in >> tmpy;
		in >> tmpz;
		in >> tmpd;
		crowd->makeNewNode(tmpx,tmpy,tmpz,tmpd);
	}

	in >> tmpn;
	for(i=0;i<tmpn;i++) {
		in >> tmpid;
		in >> tt0;
		in >> tt1;
		in >> tt2;
		n0 = crowd->getNode(tt0-1);
		n1 = crowd->getNode(tt1-1);
		n2 = crowd->getNode(tt2-1);
		tri = surface->makeNewTriangle(n0,n1,n2);
	}


	crowd->initstd();
	crowd->initDensity();

	in.close();

	return true;
}


// outRfi3dTriangle

bool outRfi3dTriangle(char *fname,Surface *surface,Crowd *crowd,char *version)
{
	int i,j;
	Node *tmpnd;
	Triangle *tri;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << "#0#0#0#0#0.000000#0#";
	out << version;
	out << "####################################################\n";

	out << "0 ";
	out << crowd->getSize() << " ";
	out << surface->getSize() << "\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId()-1 << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	n = surface->getSize();
	for(i=0;i<n;i++) {
		tri = surface->getTriangle(i);
		out << tri->getId()-1 << " 0 -1 tri";
		for(j=0;j<3;j++) {
			tmpnd = tri->getNode(j);
			out << " ";
			out << tmpnd->getId()-1;
		}
		out << "\n";
	}

	out.close();

	return true;
}

// outRfi3dTetra

bool outMSH3dTetra(char *fname,Tetgen* tetgen,Crowd* crowd)
{
	int i,j;
	Node *tmpnd;
	Tetra *ptet;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << "#FEM_MSH\n";
	out << " $PCS_TYPE\n";
	out << "  NO_PCS\n";
	out << " $NODES\n";
	out << crowd->getSize() << "\n";
   	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId()-1 << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	out << " $ELEMENTS\n";
    out << tetgen->getSize() << "\n";
	n = tetgen->getSize();
	for(i=0;i<n;i++) {
		ptet = tetgen->getTetra(i);
		out << ptet->getId()-1 <<  " "<< ptet->getDomain() << " -1 tet";
		for(j=0;j<4;j++) {
			tmpnd = ptet->getNode(j);
			out << " ";
			out << tmpnd->getId()-1;
		}
		out << "\n";
	}
	out << "#STOP\n";
	out.close();


	return true;
}


// outDtm3dTetra

bool outDtm3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Tetra *ptet;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << "\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	out << tetgen->getSize() << "\n";
	n = tetgen->getSize();
	for(i=0;i<n;i++) {
		ptet = tetgen->getTetra(i);
		out << ptet->getId() << " ";
		for(j=0;j<4;j++) {
			tmpnd = ptet->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}

	out.close();


	return true;
}


// outDli3dTetra

bool outDli3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Tetra *ptet;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << " ";
	out << tetgen->getSize() << " 0 1 0\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << " ";
		out << tmpnd->getDensity() << "\n";
	}
	n = tetgen->getSize();
	for(i=0;i<n;i++) {
		ptet = tetgen->getTetra(i);
		out << ptet->getId() << " ";
		out << ptet->getDomain() << " 0 TET ";
		for(j=0;j<4;j++) {
			tmpnd = ptet->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}

	out.close();

	return true;
}


bool outDli3dTriangle(char *fname,Surface *surface,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Triangle *ptri;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << " ";
	out << surface->getSize() << "0 1 0\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << " ";
		out << tmpnd->getDensity() << "\n";
	}
	n = surface->getSize();
	for(i=0;i<n;i++) {
		ptri = surface->getTriangle(i);
		out << ptri->getId() << " ";
		out << ptri->getDomain() << " ";
		out << ptri->getType() << " TRI ";
		for(j=0;j<3;j++) {
			tmpnd = ptri->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}

	out.close();


	return true;

	//return true;
}


// outAvs3dTetra

bool outAvs3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Tetra *ptet;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << " ";
	out << tetgen->getSize() << " 1 0 0\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	n = tetgen->getSize();
	for(i=0;i<n;i++) {
		ptet = tetgen->getTetra(i);
		out << ptet->getId() << " 1 tet ";
		for(j=0;j<4;j++) {
			tmpnd = ptet->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}
	out << " 1 1\n";
	out << " z_coordinates,\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getZ() << "\n";
	}

	out.close();

	return true;
}



bool outAvs3dTetraWithType(char *fname,Tetgen *tetgen,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Tetra *ptet;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << " ";
	out << tetgen->getSize() << " 1 1 0\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	n = tetgen->getSize();
	for(i=0;i<n;i++) {
		ptet = tetgen->getTetra(i);
		out << ptet->getId() << " 1 tet ";
		for(j=0;j<4;j++) {
			tmpnd = ptet->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}
	out << " 1 1\n";
	out << " z_coordinates,\n";
	for(i=0;i<crowd->getSize();i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getZ() << "\n";
	}
	out << " 1 1\n";
	out << " type,\n";
	n = tetgen->getSize();
	for(i=0;i<n;i++) {
		ptet = tetgen->getTetra(i);
		out << i << " ";
		out << ptet->getDomain();
		out << "\n";
	}

	out.close();

	return true;
}



// outDtm3dTriangle

bool outDtm3dTriangle(char *fname,Surface* surface,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Triangle *ptri;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << "\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	out << surface->getSize() << "\n";
	n = surface->getSize();
	for(i=0;i<n;i++) {
		ptri = surface->getTriangle(i);
		out << ptri->getId() << " ";
		for(j=0;j<3;j++) {
			tmpnd = ptri->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}

	out.close();


	return true;
}


// outAvs3dTriangle

bool outAvs3dTriangle(char *fname,Surface *surface,Crowd *crowd)
{
	int i,j;
	Node *tmpnd;
	Triangle *ptri;
	int n;
	
	ofstream out(fname);
	if(!out) {
		return false;
	}

	out << crowd->getSize() << " ";
	out << surface->getSize() << " 1 1 0\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getX() << " ";
		out << tmpnd->getY() << " ";
		out << tmpnd->getZ() << "\n";
	}
	n = surface->getSize();
	for(i=0;i<n;i++) {
		ptri = surface->getTriangle(i);
		out << ptri->getId() << " 1 tri ";
		for(j=0;j<3;j++) {
			tmpnd = ptri->getNode(j);
			out << tmpnd->getId() << " ";
		}
		out << "\n";
	}
	out << " 1 1\n";
	out << " z_coordinates,\n";
	n = crowd->getSize();
	for(i=0;i<n;i++) {
		tmpnd = crowd->getNode(i);
		out << tmpnd->getId() << " ";
		out << tmpnd->getZ() << "\n";
	}
	out << " 1 1\n";
	out << " type,\n";
	n = surface->getSize();
	for(i=0;i<n;i++) {
		ptri = surface->getTriangle(i);
		out << i+1 << " ";
		out << ptri->getDomain();
	//	out << ptri->getType();
		out << "\n";
	}

	out.close();

	return true;
}




}








//////////////////////////////////////////////////////////////////EOF