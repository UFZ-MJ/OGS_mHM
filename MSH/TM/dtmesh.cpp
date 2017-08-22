/*

			dtmesh.h
			>Dtmesh class
			>>interface of dtmesh library

			last update : 2003.12.08		by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtmesh.h"


namespace dtm{


// Dtmesh
///03.12.09

Dtmesh::Dtmesh()
: surface_refine_type(1),surface_laplas_turn(5)//,state(0)
{
	crowd = NULL;
	surface = NULL;
	tetgen = NULL;

	flg_crowd = 0;
	flg_surface = 0;
	flg_tetgen = 0;

	domain_size = 0;
}


// ~Dtmesh
///03.12.09

Dtmesh::~Dtmesh()
{
	//delete objects
	terminate();
}


// initialize
///03.12.09

void Dtmesh::initialize()
{

	//try{
		if(tetgen) delete tetgen;
		if(surface) delete surface;
		if(crowd) delete crowd;

		crowd = new Crowd;
		surface = new Surface(crowd);
		tetgen = new Tetgen(crowd);

		flg_crowd = 0;
		flg_surface = 0;
		flg_tetgen = 0;
		domain_size = 0;
		divdomain.clear();
	/*}
	catch(Errors err) {
		cout << "******ERROR******\n";
		cout << "cannot initialize\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";
		exit(0);
	}*/
}


// terminate
///03.12.09

void Dtmesh::terminate()
{
	flg_crowd = 0;
	flg_surface = 0;
	flg_tetgen = 0;
	domain_size = 0;
	divdomain.clear();
	if(tetgen) delete tetgen;
	if(surface) delete surface;
	if(crowd) delete crowd;
}



//setMaxNumberOfNodes
///03.12.09

void Dtmesh::setMaxNumberOfNodes(unsigned int n)
{
	crowd->setMaxSize(n);
}


//setMaxNumberOfTriangles
///03.12.09

void Dtmesh::setMaxNumberOfTriangles(unsigned int n)
{
	surface->setMaxSize(n);
}


//setMaxNumberOfTetrahedras
///03.12.09

void Dtmesh::setMaxNumberOfTetrahedras(unsigned int n)
{
	tetgen->setMaxSize(n);
}

//setSurfaceRefineType
///03.12.09

int Dtmesh::setSurfaceRefineType(int refine_type)
{
	if(refine_type > 3 || refine_type < 0)
		return surface_refine_type;

	surface_refine_type = refine_type;
	return surface_refine_type;
}


//setSurfaceLaplasTurn
///03.12.09

int Dtmesh::setSurfaceLaplasTurn(int laplas_turn)
{
	if(laplas_turn > MAX_LAPLAS_TURN){
		surface_laplas_turn = MAX_LAPLAS_TURN;
		return surface_laplas_turn;
	} else if( laplas_turn < 0){
		surface_laplas_turn = 0;
		return surface_laplas_turn;
	}

	surface_laplas_turn = laplas_turn;
	return surface_laplas_turn;
}


//setDomain
///03.12.08

int Dtmesh::setDomain(deque<int>&buff)
{
	domain_size++;
	divdomain.resize(domain_size);
	int n = (int)buff.size();
	for(int i=0;i<n;i++) {
		divdomain[domain_size-1].push_back(buff[i]);
	}

	return (int)divdomain.size();
}

//setDomain
///03.12.08

int Dtmesh::setDomain(int n,int buff[])
{
	domain_size++;
	divdomain.resize(domain_size);
	for(int i=0;i<n;i++) {
		divdomain[domain_size-1].push_back(buff[i]);
	}

	return (int)divdomain.size();
}

//getDomainSize
///03.12.08

int Dtmesh::getDomainSize()
{
	return domain_size;
}


//getDomainSurfaceSize
///03.12.18

int Dtmesh::getDomainSurfaceSize(int dm)
{
	int n = (int)divdomain.size();
	if(dm < 0 || dm >= n) {
		return -1;
	}
	return (int)divdomain[dm].size();
}


//getDomainNumber
///03.12.08

int Dtmesh::getDomainNumber(int dm,int index)
{
	int n = (int)divdomain.size();
	if(dm < 0 || dm >= n) {
		return -1;
	}
	int nn = (int)divdomain[dm].size();
	if(index < 0 || index >= nn) {
		return -1;
	}

	return divdomain[dm][index];
}


// Mesh
///03.12.09

int Dtmesh::Mesh(int refine_type,int laplas_turn)
{
	int i;
	Point cpt;
	int prevnode1;

//	try{

		//normalize coordinates
		crowd->normalize();

		prevnode1 = crowd->getSize();
	    //subsequent refinement of triangle pairs
		switch(setSurfaceRefineType(refine_type)) {
		case 1:
			surface->refineOnRate();
			break;
		case 2:
			surface->refineOnValue();
			break;
		case 3:
			surface->refineOnFast();
			break;
		case 0:
		default:
			break;
		}

		//laplasian triangle
		int n = setSurfaceLaplasTurn(laplas_turn);
		for(i=0;i<n;i++)
			crowd->laplasianTriangle(prevnode1,crowd->getSize()-1);

		//meshing to tetrahedra
		tetgen->delaunay(crowd);

		//removing non-convex parts
		surface->pickBoundary();
		surface->sortDirection();
		surface->selectDomain(tetgen);
		tetgen->removeTetrahedra();

		//set domain of tetrahedra
		n = (int)divdomain.size();
		for(i=0;i<n;i++) {
			surface->pickDomain(divdomain[i]);
			surface->sortDirection();
			surface->selectDomain(tetgen);
			tetgen->setDomain(i+1);
		}

		surface->clear();

		//subsequent refinement of tetrahedra
		tetgen->fineTetra();

		//restore nodes` coordinates to original value
		crowd->restore();

		flg_tetgen = 1;

/*	}
	catch(Errors err) {
		cout << "******ERROR******\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";

		terminate();

		exit(0);
	}*/
	
	return tetgen->getSize();
}


//refineSurface
///03.12.09

int Dtmesh::refineSurface(int refine_type,int laplas_turn)
{
	if(!(surface->getSize())) return 0;

	int prevnode = crowd->getSize();
	//try{
		//normalize coordinates
		crowd->normalize();

	    //subsequent refinement of triangle pairs
		switch(setSurfaceRefineType(refine_type)) {
		case 1:
			surface->refineOnRate();
			break;
		case 2:
			surface->refineOnValue();
			break;
		case 3:
			surface->refineOnFast();
			break;
		case 0:
		default:
			break;
		}

		//laplasian triangle
		int n = setSurfaceLaplasTurn(laplas_turn);
		for(int i=0;i<n;i++)
			crowd->laplasianTriangle(prevnode,crowd->getSize()-1);

	//	surface->cleanBadTriangle();

		//restore nodes` coordinates to original value
		crowd->restore();
	/*}
	catch(Errors err) {
		cout << "******ERROR******\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";

		terminate();

		exit(0);
	}*/


	return surface->getSize();
}


//meshPointToTetra
///03.12.09

int Dtmesh::meshPointToTetra()
{
	if(!(crowd->getSize())) return 0;

	//try{
		//normalize coordinates
		crowd->normalize();

		//meshing to tetrahedra
		tetgen->delaunay(crowd);

		//restore nodes` coordinates to original value
		crowd->restore();

		flg_tetgen = 1;

	/*}
	catch(Errors err) {
		cout << "******ERROR******\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";

		terminate();

		exit(0);
	}*/


	return tetgen->getSize();
}


//meshTriangleToTetra
///03.12.09

int Dtmesh::meshTriangleToTetra()
{
	if(!(surface->getSize())) return 0;
//	try{
		//normalize coordinates
		crowd->normalize();

		//meshing to tetrahedra
		tetgen->delaunay(crowd);

		//removing non-convex parts
		surface->pickBoundary();
		surface->sortDirection();
		surface->selectDomain(tetgen);
		tetgen->removeTetrahedra();
	
		//restore nodes` coordinates to original value
		crowd->restore();

		flg_tetgen = 1;
	/*}
	catch(Errors err) {
		cout << "******ERROR******\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";

		terminate();

		exit(0);
	}*/


	return tetgen->getSize();
}


//refineTetgen
///03.12.09

int Dtmesh::refineTetgen()
{
	if(!(tetgen->getSize())) return 0;

   // try{
		//normalize coordinates
		crowd->normalize();

		//subsequent refinement of tetrahedra
		tetgen->fineTetra();

		//restore nodes` coordinates to original value
		crowd->restore();
	/*}
	catch(Errors err) {
		cout << "******ERROR******\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";

		terminate();

		exit(0);
	}*/

	return tetgen->getSize();
}

//outputRfiTriangle
///03.12.09

int Dtmesh::outputRfiTriangle(char *fname,char *version)
{
	if(!flg_surface) return 0;

	if(!outRfi3dTriangle(fname,surface,crowd,version)) {
		return 0;
	}
	return 1;
}

//outputMSHTetra
///03.12.09

int Dtmesh::outputMSHTetra(char*fname)
{
	if(!flg_tetgen) return 0;

	if(!outMSH3dTetra(fname,tetgen,crowd)) {
		return 0;
	}
	return 1;
}


//outputAvsTetra
///03.12.09

int Dtmesh::outputAvsTetra(char *fname)
{
	if(!flg_tetgen) return 0;

	if(!outAvs3dTetraWithType(fname,tetgen,crowd)) {
		return 0;
	}

	return 1;
}


//outputAvsTriangle
///03.12.09

int Dtmesh::outputAvsTriangle(char *fname)
{
	if(!flg_surface) return 0;

	if(!outAvs3dTriangle(fname,surface,crowd)) {
		return 0;
	}
	
	return 1;
}


//outputDliTetra
//03.12.18

int Dtmesh::outputDliTetra(char *fname)
{
	if(!flg_tetgen) return 0;

	if(!outDli3dTetra(fname,tetgen,crowd))
		return 0;

	return 1;
}


//outputDliTetra
//03.12.18

int Dtmesh::outputDliTriangle(char *fname)
{
	if(!flg_surface) return 0;

	if(!outDli3dTriangle(fname,surface,crowd))
		return 0;

	return 1;
}


//output3dtri
///03.12.09

int Dtmesh::output3dtri(char *fname)
{
	if(!flg_surface) return 0;

	if(!outDtm3dTriangle(fname,surface,crowd)) {
		return 0;
	}

	return	1;
}


//output3dtet
///03.12.09

int Dtmesh::output3dtet(char *fname)
{
	if(!flg_tetgen) return 0;

	if(!outDtm3dTetra(fname,tetgen,crowd)) {
		return 0;
	}
	return 1;

}

// inpRfiTriangle
///03.12.09
/*
bool Dtmesh::inputRfiTriangle(char *fname)
{
	int i;
	int npt,nelm;
	int tt0,tt1,tt2;
	Node *n0,*n1,*n2;
	int id;
	double ix,iy,iz;
	double dens;
	int type1,type2;
	char etype[100];
	Triangle *itri;
	int tmp;
	char tmptitle[256];

	initialize();

	ifstream in(fname);
	if(!in) {
		return false;
	}

	in.getline(tmptitle,sizeof(tmptitle));

	in >> tmp;
	in >> npt;
	in >> nelm;

	dens = 0.;
	for(i=0;i<npt;i++) {
		in >> id;
		in >> ix;
		in >> iy;
		in >> iz;
		//in >> dens;

	//	crowd->makeNewNode(ix,iy,iz,dens);
		crowd->makeNewNode(ix,iy,iz,0.);
	}

	type1 = type2 = 0;
	for(i=0;i<nelm;i++) {
		in >> id >> type1 >> type2;
		in >> etype;
		if(!strcmpi(etype,"TRI")) {
			in >> tt0 >> tt1 >> tt2;
			n0 = crowd->getNode(tt0);
			n1 = crowd->getNode(tt1);
			n2 = crowd->getNode(tt2);
			itri = surface->makeNewTriangle(n0,n1,n2);
		} else {
			Errors err(0,
				"bool Dtmesh::inputRfiTriangle(char *fname)",
				"cannot read input file");
			throw err;
		}
	}

	crowd->initstd();
	//crowd->initDensity();
	surface->setNeighbor();

	flg_crowd = 1;
	flg_surface = 1;

	return true;
}
*/

int Dtmesh::inputMSH(CFEMesh* m_msh)
{
	int i;
	int npt,nelm;
	Node *n0,*n1,*n2,*n3;
	double ix,iy,iz;
	Triangle *itri;
	int flgtri,flgtet;
	int flg;

	initialize();

    npt = (int)m_msh->nod_vector.size();
	for(i=0;i<npt;i++) {      
        ix = m_msh->nod_vector[i]->X();
        iy = m_msh->nod_vector[i]->Y();
        iz = m_msh->nod_vector[i]->Z();
        crowd->makeNewNode(ix,iy,iz,0.);
	}

	flgtri = flgtet = 0;
    nelm = (int)m_msh->ele_vector.size();
	for(i=0;i<nelm;i++) {
        flg = m_msh->ele_vector[i]->GetElementType();
    	if(flg == 5) {
                n0 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(0));
				n1 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(1));
				n2 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(2));
				n3 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(3));
				tetgen->makeNewTetra(n0,n1,n2,n3);
				flgtet = 1;
		} 
        else if(flg == 4) {
                n0 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(0));
				n1 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(1));
				n2 = crowd->getNode(m_msh->ele_vector[i]->GetNodeIndex(2));
				itri = surface->makeNewTriangle(n0,n1,n2);
				itri->setType(1);
				itri->setDomain(0);
				flgtri = 1;
		} 
        else {
			return FILE_NON;
		}
	}


	crowd->initstd();
	if(flgtri) {
		surface->setNeighbor();
		flg_crowd = 1;
		flg_surface = 1;
		return FILE_RFI_TRI;
	}
	if(flgtet) {
		flg_crowd = 1;
		flg_tetgen = 1;
		return FILE_RFI_TET;
	}

	return FILE_NON;
}


// inputDti
///03.12.09
//input Dtmesh Library Interface file

int Dtmesh::inputDli(char *fname)
{
	int i;
	int npt,nelm;
	int tt0,tt1,tt2,tt3;
	Node *n0,*n1,*n2,*n3;
	int id;
	double ix,iy,iz;
	double dens;
	int type1,type2;
	int flg1,flg2;
	char etype[100];
	Triangle *itri;
	int ndm,nar,ar;
	int flgtri,flgtet;
	int flg;
	char buff[256];


	initialize();

	ifstream in(fname);
	if(!in) {
		return FILE_NON;
	}

	if(!in.getline(buff,sizeof(buff)))return FILE_NON;
	if(sscanf(buff,"%d%d%d%d%d",
		&npt,&nelm,&ndm,&flg1,&flg2) != 5)return FILE_NON;

	dens = 0.;
	type1 = type2 = 0;
	for(i=0;i<npt;i++) {
		if(!in.getline(buff,sizeof(buff)))return FILE_NON;
		if((flg1) && (flg2)) {
			if(sscanf(buff,"%d%lf%lf%lf%lf%d",
				&id,&ix,&iy,&iz,&dens,&type1) != 6)
				return FILE_NON;
			crowd->makeNewNode(ix,iy,iz,dens);
		} else if(flg1) {
			if(sscanf(buff,"%d%lf%lf%lf%lf%d",
				&id,&ix,&iy,&iz,&dens) != 5)
				return FILE_NON;
			crowd->makeNewNode(ix,iy,iz,dens);
		} else if(flg2) {
			if(sscanf(buff,"%d%lf%lf%lf%d",
				&id,&ix,&iy,&iz,&type1) != 5)
				return FILE_NON;
			crowd->makeNewNode(ix,iy,iz,0.);
		} else {
			if(sscanf(buff,"%d%lf%lf%lf",
				&id,&ix,&iy,&iz) != 4)
				return FILE_NON;
			crowd->makeNewNode(ix,iy,iz,0.);
		}
	}

	flgtri = flgtet = 0;
	type1 = type2 = 0;
	for(i=0;i<nelm;i++) {
		if(!in.getline(buff,sizeof(buff)))return FILE_NON;
		flg = sscanf(buff,"%d%d%d%s%d%d%d%d",
			&id,&type1,&type2,etype,&tt0,&tt1,&tt2,&tt3);
		if(flg == 8) {
			if(!strcmpi(etype,"TET")) {
				if(tt0 <= 0 || tt0 > npt) return FILE_NON;
				if(tt1 <= 0 || tt1 > npt) return FILE_NON;
				if(tt2 <= 0 || tt2 > npt) return FILE_NON;
				if(tt3 <= 0 || tt3 > npt) return FILE_NON;
				n0 = crowd->getNode(tt0-1);
				n1 = crowd->getNode(tt1-1);
				n2 = crowd->getNode(tt2-1);
				n3 = crowd->getNode(tt3-1);
				tetgen->makeNewTetra(n0,n1,n2,n3);
				flgtet = 1;
			} else {
				return FILE_NON;
			}
		} else if(flg == 7) {

			if(!strcmpi(etype,"TRI")) {
				if(tt0 <= 0 || tt0 > npt) return FILE_NON;
				if(tt1 <= 0 || tt1 > npt) return FILE_NON;
				if(tt2 <= 0 || tt2 > npt) return FILE_NON;
				n0 = crowd->getNode(tt0-1);
				n1 = crowd->getNode(tt1-1);
				n2 = crowd->getNode(tt2-1);
				itri = surface->makeNewTriangle(n0,n1,n2);
				itri->setType(type1);
				itri->setDomain(type2);
				flgtri = 1;
			} else {
				return FILE_NON;
			}
		} else {
			return FILE_NON;
		}
	}

	if(ndm > 0) {
		domain_size = ndm;
		divdomain.resize(ndm);
		for(i=0;i<ndm;i++) {
			in >> id >> nar >> type1 >> etype;
			if(strcmpi(etype,"DM"))return FILE_NON;
			for(int j=0;j<nar;j++) {
				in >> ar;
				divdomain[i].push_back(ar);
			}
		}
	}

	in.close();

	crowd->initstd();
	if(flg1) crowd->initDensity();
	if(flgtri) {
		surface->setNeighbor();
		flg_crowd = 1;
		flg_surface = 1;
	}
	if(flgtet) {
		flg_crowd = 1;
		flg_tetgen = 1;
	}

	if(!nelm) {
		if(flg1) return FILE_DLI_PTWD;
		else return FILE_DLI_PT;
	} else {
		if((flgtri) && (flgtet)) {
			if(flg1) return FILE_DLI_TRITETWD;
			else return FILE_DLI_TRITET;
		} else if(flgtri) {
			if(flg1) return FILE_DLI_TRIWD;
			else return FILE_DLI_TRI;
		} else if(flgtet) {
			if(flg1) return FILE_DLI_TETWD;
			else return FILE_DLI_TET;
		} else {
			return FILE_NON;
		}
	}

	//return FILE_NON;
}


//input3dtri
///03.12.09

int Dtmesh::inputDtm(char *fname)
{
	int i;
	int tmpnd,tmpne,tmpid;
	double tmpx,tmpy,tmpz;
	double tmpd;
	int tt0,tt1,tt2,tt3;
	Node *n0,*n1,*n2,*n3;
	Triangle *tri;
	int flg;
	int flgtet,flgtri;
	char buff[256];
	int flgdens;
	initialize();


	ifstream in(fname);
	if(!in) {
		return FILE_NON;
	}

	if(!in.getline(buff,sizeof(buff)))return FILE_NON;
	if(sscanf(buff,"%d",&tmpnd) != 1) return FILE_NON;

	flgdens = 0;
	for(i=0;i<tmpnd;i++) {
		if(!in.getline(buff,sizeof(buff)))return FILE_NON;
		flg = sscanf(buff,"%d%lf%lf%lf%lf",&tmpid,&tmpx,&tmpy,&tmpz,&tmpd);
		if(flg == 5) {
			crowd->makeNewNode(tmpx,tmpy,tmpz,tmpd);
			flgdens = 1;
		} else if(flg == 4) {
			crowd->makeNewNode(tmpx,tmpy,tmpz,0.);
		} else {
			return FILE_NON;
		}
	}

	flgtet = flgtri = 0;
	if(!in.getline(buff,sizeof(buff))){
		if(flgdens) {
			crowd->initDensity();
			return FILE_DTM_PTWD;
		} else {
			return FILE_DTM_PT;
		}
	}
	if(sscanf(buff,"%d",&tmpne) != 1) return FILE_NON;
	//in >> tmpn;
	for(i=0;i<tmpne;i++) {
		if(!in.getline(buff,sizeof(buff)))return FILE_NON;
		flg = sscanf(buff,"%d%d%d%d%d",&tmpid,&tt0,&tt1,&tt2,&tt3);
		if(flg == 5) {
			if(tt0 < 1 || tt0 > tmpnd) return FILE_NON;
			if(tt1 < 1 || tt1 > tmpnd) return FILE_NON;
			if(tt2 < 1 || tt2 > tmpnd) return FILE_NON;
			if(tt3 < 1 || tt3 > tmpnd) return FILE_NON;
			n0 = crowd->getNode(tt0-1);
			n1 = crowd->getNode(tt1-1);
			n2 = crowd->getNode(tt2-1);
			n3 = crowd->getNode(tt3-1);
			tetgen->makeNewTetra(n0,n1,n2,n3);
			flgtet = 1;
		} else if(flg == 4) {
			if(tt0 < 1 || tt0 > tmpnd) return FILE_NON;
			if(tt1 < 1 || tt1 > tmpnd) return FILE_NON;
			if(tt2 < 1 || tt2 > tmpnd) return FILE_NON;
			n0 = crowd->getNode(tt0-1);
			n1 = crowd->getNode(tt1-1);
			n2 = crowd->getNode(tt2-1);
			tri = surface->makeNewTriangle(n0,n1,n2);
			tri->setType(1);
			tri->setDomain(0);
			flgtri = 1;
		} else {
			return FILE_NON;
		}
	}

	in.close();


	crowd->initstd();
	if(flgtri) {
		surface->setNeighbor();
		flg_crowd = 1;
		flg_surface = 1;
		if(flgdens) {
			crowd->initDensity();
			return FILE_DTM_TRIWD;
		} else {
			return FILE_DTM_TRI;
		}
	} else if(flgtet) {
		flg_crowd = 1;
		flg_tetgen = 1;
		if(flgdens) {
			crowd->initDensity();
			return FILE_DTM_TETWD;
		} else {
			return FILE_DTM_TET;
		}
	} 

	return FILE_NON;
}


///for debug////
void Dtmesh::view()
{
	cout << "NODE : " << crowd->getSize();
	cout << " TRI : " << surface->getSize();
	cout << " TET : " << tetgen->getSize();
	cout << "\n";

	return;
}

}



//////////////////////////////////////////////////////////////////EOF