#include "tricommon.h"


static bool bStartTriangle_TRI( 
		std::vector<Cp_dbl2>	*Pnt2List,		// [i] list of node
		std::vector<Cp_dtri> *TriList,		// [i] list of triangles
		int				iUseConvexHul )	// [i] 1:YES,0:NO)
{	
	bool bRet, bRetCode = false;
	int iChPoints, i, j, k, opoints, npoints, nInPoints;
	int *iBelongCount = NULL;
	int *iBoxSorted = NULL;
	int *iNbrCount = NULL;

	Cp_dbl2 ep;
	CDT_Matrix<int> iNbrList;
	CDT_Matrix<int> iBlgList;
	std::vector<Cp_dtri> TriVector;
	std::vector<int> Cnv2List;
	std::vector<int> PntTable;

	TriList->clear( );
	
	// --- ptables: point number in triangulation <-> original number table. ---

	npoints = (int)Pnt2List->size( );
	nInPoints = npoints;

	// --- temporarily remove points with identical X and Y coordinates ---

	PntTable.clear( );
	for ( i = 0; i < npoints + iPsuedoPointCount_TRI( ); i++ ) {
		PntTable.push_back( i );
	}

	fprintf( stdout, "\rSorting %d Points....", npoints );
	fflush( stdout );

	for ( i = 0; i < npoints - 1; i++ ) {
		for ( j = i + 1; j < npoints; j++ ) {
			if( Pnt2List->at( i ).x != Pnt2List->at( j ).x ) continue;
			if( Pnt2List->at( i ).y != Pnt2List->at( j ).y ) continue;

			Pnt2List->erase( Pnt2List->begin( ) + j );
			PntTable.erase( PntTable.begin( ) + j );
			
			npoints--;
			j--;
		}
	}
	fprintf( stdout, "Done.\n" );
	

	i = (int)Pnt2List->size( );
	for( i = 0 ; i < npoints ; i++ ) {
		ep = Pnt2List->at( i );
	}

	if( 3 > npoints ) {
		printf( "Less than 3 unique points! %d\n ", npoints );
		goto PIX_EXIT;
	}

	/*------------------------------------------------------------------
	#### memory allocation for the neighbouring points list.
	#### 
	#### opoints	
	#### number of original points, excluding the duplicated points ( if any ).
	#### 
	#### npoints
	#### number of points including the pseudo points.
	#### 
	#### iNbrCount[i]	 
	#### number of points connected to point i including "i" itself.
	------------------------------------------------------------------*/
	opoints = npoints;
	npoints += iPsuedoPointCount_TRI( );

	iNbrCount = new int[npoints];
	if( NULL == iNbrCount ) {
		printf( "\n\7Failed in Memory Allocation.\n" );
		goto PIX_EXIT;
	}

	// --- initilize neibur list ---

	if( ! iNbrList.Alloc( npoints, npoints ) ) {
		printf( "Memory allocation Failed!" );
		goto PIX_EXIT;
	}

	for ( i = 0; i < npoints; i++ ) iNbrCount[i] = 0;	
	if( !bInitNeighbourList( npoints, &iNbrList ) ) {
		goto PIX_EXIT;
	}
	
	// --- iBelongCount[i]	number of triangles to which node i belongs ---
	// --- iBlgList[i][j]	j-th triangle to which node i belongs ---
	
	// --- initialize belong list ---

	if( ! iBlgList.Alloc( npoints, npoints ) ) {
		printf( "Memory allocation Failed!" );
		goto PIX_EXIT;
	}
	
	iBelongCount = new int[npoints];
	if( NULL == iBelongCount ) {
		printf( "\n\7Failed in Memory Allocation.\n" );
		goto PIX_EXIT;
	}

	for ( i = 0; i < npoints; i++ ) iBelongCount[i] = 0; 
	if( !bInitBelongList_TRI( npoints, &iBlgList ) ) {
		goto PIX_EXIT;
	}

	// --- generate pseudo points‚ð”­¶‚·‚é ---
	
	bRet = bSetPseudoPoints_TRI( opoints, Pnt2List );

	// --- BOX-sort points ---
	
	iBoxSorted = new int[npoints];
	if( NULL == iBoxSorted ) {
		printf( "Failed in Memory Allocation.\n" );
		goto PIX_EXIT;
	}

	bRet = bBoxSortPoint_TRI( opoints, Pnt2List, iBoxSorted );

	// --- START!! constructing neighbour and triangle list. ---

	TriVector.clear( );
	for( i = 0; i < opoints; i++ ) {
		fprintf( stdout, "\rProcessing Point %d/%d....", i, npoints );
		fflush( stdout );
		k = iBoxSorted[i];
		// --- get triangle and neigbour list around point "p". ---
		iNbrCount[k] = iListNeighbour_TRI( k, opoints, Pnt2List, 
			iBelongCount, &iNbrList, &iBlgList, &TriVector );
		if( iNbrCount[k] < 0 ) { 
			printf( "\7\n" ); 
			break; 
		}
	}
	fprintf( stdout, "Done.\n" );


	// --- process points on ConvexHull ---

	if( iUseConvexHul ) {	
		
		// --- generate convex hull ---	
		Cnv2List.clear( );
		iChPoints = CreateConvexHull( opoints, Pnt2List, &Cnv2List );
		iChPoints = (int)Cnv2List.size( );

		// --- process points on ConvexHull ---
		for( i = 0; i < iChPoints; i++ ) {
			fprintf( stdout, "\rProcessing Convex Hull %d/%d....", i, iChPoints );
			iListConvex_TRI( i, opoints, Pnt2List, &Cnv2List,  
				iNbrCount, &iNbrList, &TriVector );
		}
		fprintf( stdout, "Done.\n" );
	}

	// --- reconstruct triangle list ---

	fprintf( stdout, "\rWriting %d Triangles....", (int)TriVector.size( ) );
	bRet = bResetTriangleList_TRI( nInPoints, &PntTable, &TriVector, TriList );
	fprintf( stdout, "Done.\n" );

	// --- Done ---

	bRetCode = true;

PIX_EXIT:

	if( iBelongCount ) delete[] iBelongCount;
	if( iBoxSorted ) delete[] iBoxSorted;
	if( iNbrCount ) delete[] iNbrCount;
	
	iNbrList.Free( );
	iBlgList.Free( );
	Cnv2List.clear( );
	TriVector.clear( );
	PntTable.clear( );

	return bRetCode;
}

 bool bLoadStart_TRI(     // Function to be called!!!
		char	*szNodeFile,	// [i] input node file
		char	*szTriFile,		// [i] triangle file to be created
		int		iUseConvexHul )	// [i] 1:YES,0:NO
{
	bool bRet, bRetCode = false;
	int nTri;
	std::vector<Cp_dbl2> Pnt2List;	
	std::vector<Cp_dbl3> Pnt3List;	
	std::vector<Cp_dtri> TriList;	

	// --- load points ---
	bLoadNodeFile_TRI( szNodeFile, &Pnt3List );
	ConvertPntList( &Pnt3List, &Pnt2List );

	// --- start triangulation ---
	
	bStartTriangle_TRI( &Pnt2List, &TriList, iUseConvexHul );
	nTri = (int)TriList.size( );

	// --- save triangles ---

	bRet = bSaveTriFile_TRI( szTriFile, &TriList, &Pnt3List );
	if( !bRet ) {
		goto PIX_EXIT;
	}
		
	// --- Done ---

	bRetCode = true;

PIX_EXIT:

	Pnt2List.clear( );
	Pnt3List.clear( );
	TriList.clear( );

	return bRetCode;
}


void usage( char *argv[] )
{
	fprintf( stdout, "Usage: %s ", argv[0] );
	fprintf( stdout, "[Node File] " );
	fprintf( stdout, "(Use Convex Hull? Y/N ) " );
	fprintf( stdout, "\n" );
	fprintf( stdout, "===> Notes:\n" );
	fprintf( stdout, "===> [Node File] must be a name of existing Node file.\n" );
	fprintf( stdout, "===> [Node File] has a list of (X,Y) or (X,Y,Z) coordiinated \n" );
	fprintf( stdout, "===> delineated by comma(,). Output triangle file is created \n" );
	fprintf( stdout, "===> in the same directory with an extention 'tri'.\n" );
	fprintf( stdout, "===> Option (Use Convex Hull?) determins whether all area in \n" );
	fprintf( stdout, "===> the convex hull are triangulated. This option must be Y or N. \n" );
	fprintf( stdout, "===> If this option is omitted, convex hull is not used.\n" );
}

#ifdef WIN32
void main_program_process( int argc, char *argv[] )
{
	char szTriFile[FILE_NAME_SIZE], szNodeFile[FILE_NAME_SIZE];
	char szFileName[FILE_NAME_SIZE], szDir[FILE_NAME_SIZE];
	char szDrive[FILE_NAME_SIZE], szExt[FILE_NAME_SIZE];
	int iUseConvexHul;
	
	if( 2 == argc ) {
		strcpy( szNodeFile, argv[1] );
		iUseConvexHul = 0;
	} else if( 3 == argc ) {
		strcpy( szNodeFile, argv[1] );
		if( 'Y' == argv[2][0] || 'y'  == argv[2][0] ) {
			iUseConvexHul = 1;
		} else if( 'N' == argv[2][0] || 'n'  == argv[2][0] ) {
			iUseConvexHul = 0;
		} else {
			usage( argv );
			return;
		}
	} else {
		usage( argv );
		return;
	}

	_splitpath( szNodeFile, szDrive, szDir, szFileName, szExt );
	_makepath( szTriFile, szDrive, szDir, szFileName, ".tri" );

	if( bLoadStart_TRI( szNodeFile, szTriFile, iUseConvexHul ) ) {
		fprintf( stdout, "Triangle File %s Created.\n", szTriFile );
	} else {
		fprintf( stdout, "Triangulation Failed!\n" );
	}

	return;
}
#endif // WIN32

