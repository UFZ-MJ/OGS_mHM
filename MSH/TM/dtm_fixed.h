/*

			dtm_fixed.h
			>fixed parameter

			last update : 2003.11.28		by t.manabe

*/



#ifndef DTM_DTMFIXED_H
#define DTM_DTMFIXED_H

namespace dtm{

#define LIMIT		0.000000000001
#define LIMIT_10	0.0000000001
#define LIMIT_8		0.00000001
#define LIMIT_6		0.000001

#define FRAT_ANGLE	0.99800
//#define FRAT_ANGLE	0.99900


#define NUM_FG_TYPE		2
#define FG_TYPE_TRI		0
#define FG_TYPE_TET		1
//#define FG_TYPE_HEX		2	//NUM_FG_TYPE 3
//#define FG_TYPE_LNE		3	//NUM_FG_TYPE 4

#define ND_NUM_TRI		3
#define ND_NUM_TET		4
//#define ND_NUM_HEX		8
//#define ND_NUM_LNE		2

//#define NB_NUM_TRI		3
//#define NB_NUM_TET		4
//#define NB_NUM_HEX		6


#define MAX_LAPLAS_TURN		100

#define LIMIT_NUMBER_OF_ELEMENTS		2000000
#define LIMIT_NUMBER_OF_NODES			500000
#define LIMIT_NUMBER_OF_TRIANGLES		1000000
#define LIMIT_NUMBER_OF_TETRAHEDRAS		2000000

//define file type
#define FILE_NON			0

#define FILE_DTM_PT			10
#define FILE_DTM_PTWD		15
#define FILE_DTM_TRI		12
#define FILE_DTM_TRIWD		17
#define FILE_DTM_TET		13
#define FILE_DTM_TETWD		18


#define FILE_RFI_TRI		21
#define FILE_RFI_TET		22

#define FILE_AVS_TRI		31
#define FILE_AVS_TET		32

/* 
	DLI FILE NUMBER
	40...44		: non density
	45...49		: with density
	N%5 == 0	: point
	N%5 == 1	: not use
	N%5 == 2	: triangle
	N%5 == 3	: tetrahedra
	N%5 == 4	: triangle and tetrahedra
	*/
#define FILE_DLI_PT			40
#define FILE_DLI_PTWD		45
#define FILE_DLI_TRI		42
#define FILE_DLI_TRIWD		47
#define FILE_DLI_TET		43
#define FILE_DLI_TETWD		48
#define FILE_DLI_TRITET		44
#define FILE_DLI_TRITETWD	49

}

#endif

//////////////////////////////////////////////////////////////////EOF