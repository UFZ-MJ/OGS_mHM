/*

			dtmstdio.h
			>input and output

			last update : 2003.11.01		by t.manabe

*/


#ifndef DTM_STDIO_H
#define DTM_STDIO_H	5


#include<iostream>
#include<string>
#include<cstdio>

#include<fstream>

#include"dtm_crowd.h"
#include"dtm_tetgen.h"
#include"dtm_surface.h"
using namespace std;

namespace dtm{

// * Dtm * ... for DTM files(*.3dtri,*.3dtet,*.3dpt)
// * Avs * ... for MicroAvs(*.inp)
// * Rfi * ... for RFI files(*.rfi)

	//input file
	bool inpDtm3dPoint(char *fname,Crowd *crowd);
	bool inpDtm3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd);
	bool inpDtm3dTriangle(char *fname,Surface *surface,Crowd *crowd);
	bool inpDtm3dTriangleWithDensity(char *fname,Surface *surface,Crowd *crowd);

	//output file
	bool outRfi3dTriangle(char *fname,Surface *surface,Crowd *crowd,char *version);
	bool outMSH3dTetra(char *fname,Tetgen* tetgen,Crowd* crowd);

	bool outDli3dTriangle(char *fname,Surface *surface,Crowd *crowd);
	bool outDli3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd);

	bool outDtm3dTriangle(char *fname,Surface* surface,Crowd *crowd);
	bool outDtm3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd);

	bool outAvs3dTriangle(char *fname,Surface *surface,Crowd *crowd);
	bool outAvs3dTetraWithType(char *fname,Tetgen *tetgen,Crowd *crowd);
	bool outAvs3dTetra(char *fname,Tetgen *tetgen,Crowd *crowd);
	

}



#endif

//////////////////////////////////////////////////////////////////EOF