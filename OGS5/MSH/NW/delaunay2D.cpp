// delaunay2D.cpp : Okayama 2007 2D Triangulation.
//

#include "stdafx.h"
#include <stdio.h>
#include <string>
#include <algorithm>
#include <cctype>
#include "delaunay2D.h"
#include "fdelaun2D.h"
#include "geo_pnt.h"
#include "geo_ply.h"
#include "geo_sfc.h"
#include <fstream>

using namespace std;

#ifdef MFC
//#include <windows.h>
#include <process.h>	//thread
#endif

#pragma comment(lib, "fdelaun2d-msvc.lib")

volatile static int returnCodeForThread;	//thread


/**************************************************************************/
/* MSH - Function: getIndexInList
/* Task:
/* History:
   03/2007     NW        Implementation
**************************************************************************/
int getIndexInList(vector<CGLPoint*> &vct_points, int vertexId) 
{
	for (unsigned int i=0; i<vct_points.size(); i++) {
		if (vertexId == vct_points.at(i)->id) {
			return i;
		}
	}
	return -1;
}

/**************************************************************************/
/* MSH - Function: calcTriangleSignedArea
/* Task:
/* History:
   03/2007     NW        Implementation
**************************************************************************/
double calcTriangleSignedArea(CGLPoint* p1, CGLPoint* p2)
{
	double area = (p1->x * p2->y) - (p2->x * p1->y);
	area *= 0.5;
	return area;
}

/**************************************************************************/
/* MSH - Function: calcPolylineSignedArea
/* Task:
/* History:
   03/2007     NW        Implementation
**************************************************************************/
double calcPolylineSignedArea(vector<CGLPoint*> &polyline)
{
	double area = 0.0;
	int n = (int)polyline.size();
	for (int i=0; i<n; i++) {
		area += calcTriangleSignedArea(polyline[i], polyline[(i+1)%n]);
	}
	return area;
}

/**************************************************************************/
/* MSH - Function: isPolylineClockwise
/* Task:
/* History:
   03/2007     NW        Implementation
**************************************************************************/
bool isPolylineClockwise(vector<CGLPoint*> &polyline)
{
	double dSignedArea = calcPolylineSignedArea(polyline);
	return dSignedArea < 0;
}

/**************************************************************************/
/* MSH - Function: isPolylineClockwise
/* Task:
/* History:
   03/2007     NW        Implementation
**************************************************************************/
bool isClosedSurface(Surface* psurface)
{
	int beginVertexID = psurface->polyline_of_surface_vector.at(0)->point_vector.front()->id;
	int endVertexID = psurface->polyline_of_surface_vector.at(0)->point_vector.back()->id;

	int n = (int)psurface->polyline_of_surface_vector.size();
	for (int i=1;i<n; i++) {
		if (psurface->polyline_of_surface_vector.at(i)->point_vector.front()->id == endVertexID) {
			endVertexID = psurface->polyline_of_surface_vector.at(i)->point_vector.back()->id;
		} else if (psurface->polyline_of_surface_vector.at(i)->point_vector.back()->id == endVertexID) {
			endVertexID = psurface->polyline_of_surface_vector.at(i)->point_vector.front()->id;
		} else {
			return false;
		}
	}

	return (beginVertexID == endVertexID);
}

/**************************************************************************/
/* MSH - Function: isOpenBoundary
/* Task:
/* History:
   03/2007     NW        Implementation
**************************************************************************/
bool isOpenBoundary(CGLPolyline *polyline)
{
	string polyline_name = polyline->name;
	std::transform(polyline_name.begin(), polyline_name.end(), polyline_name.begin(), std::toupper);

	const unsigned int postfix_size = (unsigned int)strlen(OPEN_BOUNDARY_POSTFIX);
	if (polyline_name.length() < postfix_size
		|| polyline_name.substr(polyline_name.length() - postfix_size, postfix_size).compare(OPEN_BOUNDARY_POSTFIX)!=0) {
		return false;
	}
	return true;
}

/**************************************************************************/
/* MSH - Function: ReadInputData
                                                                          */
/* Task: Read Function for gli data

/* History:
   03/2007     NW        1st Coding
**************************************************************************/
int ReadInputData(int *nex, int *nbk, int *nob, int *nib, int *ibex, int *ibno, int *nidm, int *ibreak, int *nbreak, 
				  double *px, double *py, double *pd)
{
	vector<CGLPoint*> vct_points;

	Surface* psurface = NULL;
	*nex = 0;
	*nbk = 0;
	for (int i=0; i<(int)surface_vector.size(); i++) {
		psurface = surface_vector.at(i);
		if (psurface->meshing_allowed == 0) {
			continue;
		}

		int numberOfpolyline = (int)psurface->polyline_of_surface_vector.size();
		if (numberOfpolyline == 0) {
			continue;
		}

		if (isClosedSurface(psurface)) {
			//#CLOSED BOUNDARY
			if(*nex == KBD){
				return 0;
			}
			ibex[*nex] = (int)psurface->polygon_point_vector.size();
			nidm[*nex] = i+1;
			bool clockwise = isPolylineClockwise(psurface->polygon_point_vector);
			for (int j=0; j<ibex[*nex]; j++) {
				CGLPoint* pt;
				if (clockwise) {
					pt = psurface->polygon_point_vector.at(j);
				} else {
					pt = psurface->polygon_point_vector.at(ibex[*nex]-j-1);
				}
				if (getIndexInList(vct_points, pt->id) == -1) {
					vct_points.push_back(pt);
				}
				*(ibno+KBD*j+*nex) = pt->id;
			}

			*nex = *nex + 1;
		}
	}

	for (int i=0; i<(int)polyline_vector.size(); i++) {
		if (!isOpenBoundary(polyline_vector.at(i))) {
			continue;
		}
		//#OPEN BOUNDARY
		if(*nbk == KBD){
			return 0;
		}
		ibreak[*nbk] = (int)polyline_vector.at(i)->point_vector.size();
		for (int j=0; j<ibreak[*nbk]; j++) {
			CGLPoint* pt = polyline_vector.at(i)->point_vector.at(j);
			if (getIndexInList(vct_points, pt->id) == -1) {
				vct_points.push_back(pt);
			}
			*(nbreak+KBD*j+*nbk) = pt->id;
		}

		*nbk = *nbk + 1;
	}

	//Reset vertex id
	for (int i=0; i<*nex; i++) {
		for (int j=0; j<ibex[i]; j++) {
			int vertexId = *(ibno+KBD*j+i);
			*(ibno+KBD*j+i) = getIndexInList(vct_points, vertexId)+1;
		}
	}
	for (int i=0; i<*nbk; i++) {
		for (int j=0; j<ibreak[i]; j++) {
			int vertexId = *(nbreak+KBD*j+i);
			*(nbreak+KBD*j+i) = getIndexInList(vct_points, vertexId)+1;
		}
	}
	
	//#POINT
	*nob = (int)vct_points.size();
	*nib = 0;

	for (unsigned int i=0; i<vct_points.size(); i++) {
		px[i] = vct_points.at(i)->x;
		py[i] = vct_points.at(i)->y;
		pd[i] = vct_points.at(i)->mesh_density;
	}

	return -1;
}

/**************************************************************************/
/* MSH - Function: WriteOutputFile
                                                                          */
/* Task: Output function

	 
/* History:
   02/2007     NW        1st Coding
   05/2007     NW        Modified how to set the material no.
**************************************************************************/
int WriteOutputFile(char* filepath, int *node, double *x, double *y, int *nelm, int *mtj, int *idm)
{
	FILE *msh_file;
	int i = 0;
	msh_file = fopen(filepath, "w");
	if (msh_file == NULL) {
		printf("ERROR: OUTPUT FILE OPEN ERROR - %s\n",filepath);
        return 0;
	}
    //write MSH Head
	   	fprintf( msh_file, "%s\n", "#FEM_MSH");
    //write PCS Type
        fprintf( msh_file, "%s\n", " $PCS_TYPE");
	   	fprintf( msh_file, "%s\n", "  NO_PCS");
    //Write Nodes
	   	fprintf( msh_file, "%s\n", " $NODES");
        fprintf(msh_file,"% ld\n",*node);
    //write Geometry    
	for (i=0; i<*node; i++) {
        fprintf( msh_file, "%d ", i);
		fprintf(msh_file, "%15.7lf%15.7lf%15.7lf\n",x[i],y[i],0.0);
	}
    //Write Elements
	   	fprintf( msh_file, "%s\n", " $ELEMENTS");
		fprintf(msh_file,"% ld\n",*nelm);
    //write Topology
    for (i=0; i<*nelm; i++) {
        fprintf( msh_file, "%d ", i);
        fprintf( msh_file, " %d tri ", idm[i]-1);	//NW Marerial No. starts zero 
    	fprintf(msh_file, "%d %d %d\n",*(mtj+i)-1,*(mtj+KTE+i)-1,*(mtj+KTE*2+i)-1);	
		/*
		//Prepared for Neighbourhood write out but no read function for this yet
		// There is a neighbourhood data structrure filled when reading the MSH file
		fprintf(msh_file, "%7d%7d%7d%7d%7d%7d%7d%7d\n",
			*(mtj+i)-1,*(mtj+*maxelm+i)-1,*(mtj+(*maxelm)*2+i)-1,*(mtj+(*maxelm)*3+i)-1,
			*(jac+i)-1,*(jac+*maxelm+i)-1,*(jac+(*maxelm)*2+i)-1,*(jac+(*maxelm)*3+i)-1);
	    */
	}
    //Write File Terminator
	fprintf( msh_file, "%s\n", "#STOP");
	
	fclose(msh_file);
	return -1;
}

/**************************************************************************/
/* MSH - Function: WriteOutputFileBIN
                                                                          */
/* Task: Output binary mesh file

	 
/* History:
   04/2007     NW        1st Coding
**************************************************************************/
int WriteOutputFileBIN(char* filepath, int *node, double *x, double *y, int *nelm, int *mtj, int *idm)
{
	fstream *fem_msh_file = new fstream(filepath, ios::binary|ios::out);
	if( fem_msh_file->fail() ) {
		printf("ERROR: OUTPUT FILE OPEN ERROR - %s\n",filepath);
        return 0;
	}
	//KEYWORD
	char binary_char_9[9] = "#FEM_MSH";
	fem_msh_file->write((char*)(&binary_char_9),sizeof(char[9]));
	//--------------------------------------------------------------------
	// TYPE
	char binary_char_10[10] = "$PCS_TYPE";
	fem_msh_file->write((char*)(&binary_char_10),sizeof(char[10]));
    char dyn_char_G[17] = "GROUNDWATER_FLOW";
    fem_msh_file->write((char*)(&dyn_char_G),sizeof(char[17]));
	//--------------------------------------------------------------------
	// NODES
	char binary_char_7[7] = "$NODES";
	fem_msh_file->write((char*)(&binary_char_7),sizeof(char[7]));
	long binary_long = (long)(*node);
	fem_msh_file->write((char*)(&binary_long),sizeof(binary_long));
	double binary_double = 0.0;
	for(int i=0;i<(long)*node;i++){
		fem_msh_file->write((char*)(&i),sizeof(long));
		fem_msh_file->write((char*)(&x[i]),sizeof(double));
		fem_msh_file->write((char*)(&y[i]),sizeof(double));
		binary_double = 0.0;
		fem_msh_file->write((char*)(&binary_double),sizeof(double));
	}
	//--------------------------------------------------------------------
	// ELEMENTS
	char binary_char[10] = "$ELEMENTS";
	fem_msh_file->write((char*)(&binary_char),sizeof(binary_char));
	binary_long = (long)(*nelm);
	fem_msh_file->write((char*)(&binary_long),sizeof(binary_long));
	int binary_int=0;
	for(int i=0; i<(*nelm); i++){
		fem_msh_file->write((char*)(&i),sizeof(long));
		binary_int = idm[i] -1;
		fem_msh_file->write((char*)(&binary_int),sizeof(int));
		binary_long = -1;
		fem_msh_file->write((char*)(&binary_long),sizeof(binary_long));
		fem_msh_file->write((char*)("tri"),sizeof(char[4]));
		for(int j=0;j<3;j++){
			long node_id = *(mtj+KTE*j+i)-1;
			fem_msh_file->write((char*)(&node_id),sizeof(long));
		}
	}

	fem_msh_file->flush();
	fem_msh_file->close();

	return -1;
}

/**************************************************************************/
/* MSH - Function: ExecuteDelaunay2D
                                                                          */
/* Task: Interface function for thread  
/* History:
   03/2007     NW        Implementation
**************************************************************************/
unsigned __stdcall ExecuteDelaunay2DThread(void* outfilepath)
{
	returnCodeForThread = ExecuteDelaunay2DProcess((char*)outfilepath);
	return 0;
}

/**************************************************************************/
/* MSH - Function: ExecuteDelaunay2DProcess
                                                                          */
/* Task: Starts and Executes the 2D Delaunay mesh generator   
   char* outfilepath:  pointer of output path

/* fdelaun2d input data	                                                  
   ktj:		max number of node
   kbd:     max number of boundary
   kcm:
   nex:		number of closed boundary
   ibex:	1d array, number of vertex on closed boudary
   nidm:    1d array, domain id of closed boundary
   ibno:	2d array, vertex id list of closed boundary
   nib:		1d array, domain id list of closed boundary
   nbk:		number of open boundary
   ibreak:	1d array, number of vertex on open boudary
   nbreak:	2d array, vertex id list of open boundary
   nob:		number of vertex on boundary(both closed, open)
   nib:		number of special points 
   px,py:   2d-coordinates(X,Y) array of points
   pd:      point density
   dpp:     minimum distance between points (ex: 1.0E-7)
   stl:     standard element size
 
/* fdelaun2d output data	                                                  
   nelm:    number of elements
   mtj:     2d array, element-node relation array
   jac:     2d array, neighbor element relation array
   idm:     1d array, 
   node:    number of node
   px,py:   coordinate of points

/* History:
   02/2007     NW        1st Coding
**************************************************************************/
int ExecuteDelaunay2DProcess(char* outfilepath)
{
	int nex, nbk, nob, nib, node, nelm, ierrcode;
	int *ibex, **ibno, *nidm, *ibreak, **nbreak;
	int *mtj, *jac, *idm;
	double *px, *py, *pd;
	double dpp,stl;
	int ktj,kbd,kcm;

	nex = 0;
	nbk = 0;
	nob = 0;
	nib = 0;
	node = 0;
	nelm = 0;
	ierrcode = 0;
	ibex = (int*) malloc(sizeof(int)*KBD);
	ibno = (int**) malloc(sizeof(int)*KBD*KTJ);
	nidm = (int*) malloc(sizeof(int)*KBD);
	ibreak = (int*) malloc(sizeof(int)*KBD);
	nbreak = (int**) malloc(sizeof(int)*KBD*KTJ);
	mtj = (int*) malloc(sizeof(int)*KTE*3);
	jac = (int*) malloc(sizeof(int)*KTE*3);
	idm = (int*) malloc(sizeof(int)*KTE);
	px = (double*) malloc(sizeof(double)*(KTJ+3));
	py = (double*) malloc(sizeof(double)*(KTJ+3));
	pd = (double*) malloc(sizeof(double)*(KTJ+3));
	dpp = 1.0e-7;
	stl = 1.0;
	ktj = KTJ;
	kbd = KBD;
	kcm = KCM;
	for (int i=0; i<KBD; i++) {
		*(ibex+i) = 0;
		*(nidm+i) = 0;
		*(ibreak+i) = 0;
	}
	for (int i=0; i<KBD*KTJ; i++) {
		*(ibno+i) = 0;
		*(nbreak+i) = 0;
	}
	for (int i=0; i<KTE; i++) {
		*(idm+i) = 0;
	}
	for (int i=0; i<KTE*3; i++) {
		*(mtj+i) = 0;
		*(jac+i) = 0;
	}
	for (int i=0; i<KTJ+3; i++) {
		*(px+i) = 0.0;
		*(py+i) = 0.0;
		*(pd+i) = 0.0;
	}


	if (ibex==NULL || ibno==NULL || nidm==NULL || ibreak==NULL || nbreak==NULL 
		|| mtj==NULL || jac==NULL || idm==NULL || px==NULL || py==NULL || pd==NULL) {
		printf("ERROR: Fail to malloc()");
		return 999;
	}
	if(ReadInputData(&nex, &nbk, &nob, &nib, ibex, (int*)ibno, nidm, ibreak, (int*)nbreak, px, py, pd)) {
		triangulate(&kbd, &ktj, &kcm, &nex, ibex, ibno, nidm, &nbk, ibreak, nbreak, &nob, &nib, 
					px, py, pd, &dpp, &stl, &node, &nelm, mtj, jac, idm, &ierrcode);

		if (ierrcode == 0) {
			if (nelm < 1E6) {
				WriteOutputFile(outfilepath, &node,px,py,&nelm,(int*)mtj,idm);
			} else {
				char bin_out_file_path[256];
				strncpy(bin_out_file_path, outfilepath, sizeof(bin_out_file_path));
				char* retbuf = strstr(bin_out_file_path, ".msh");
				strcpy(retbuf, "_binary.msh");
				WriteOutputFileBIN(bin_out_file_path, &node,px,py,&nelm,(int*)mtj,idm);
			}
		}
	}

	free(ibex);
	free(ibno);
	free(nidm);
	free(ibreak);
	free(nbreak);
	free(mtj);
	free(jac);
	free(idm);
	free(px);
	free(py);
	free(pd);

	return ierrcode;
}


/**************************************************************************/
/* MSH - Function: ExecuteDelaunay2D
                                                                          */
/* Task: Interface for outside   

/* History:
   03/2007     NW        Implementation
**************************************************************************/
int ExecuteDelaunay2D(char* outfilepath)
{
#ifdef MFC
	HANDLE hThread;
    unsigned int dwThreadId;

	returnCodeForThread = 0;

    //Begin thread for mesh generation process
    hThread = (HANDLE)_beginthreadex(
							NULL,/* Security data structure */ 
							DELAUN2D_STACK_SIZE,/* Stack size */ 
							ExecuteDelaunay2DThread, /* start function address */
							(void*)outfilepath,/* arglist */
							0,/* initflag */
							&dwThreadId /* thread id */
							);

    WaitForSingleObject(hThread,INFINITE);  
    CloseHandle(hThread);

	return returnCodeForThread;
#else
	return ExecuteDelaunay2DProcess(outfilepath);
#endif

}


