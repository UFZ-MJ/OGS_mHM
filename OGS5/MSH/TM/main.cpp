/*

			main.cpp
			>main file of DTMesh

			last update : 2003.12.09		by t.manabe

*/


#include "stdafx.h" /* MFC */


#include"dtm_fixed.h"
#include"dtmesh.h"


void main()
{
	using dtm::Dtmesh;
	Dtmesh* dtmesh;
	char fname[256];//	char version[100];
	int flg;
	int type;
	int i,n;
	double density;


	cout << "+------------------------------------------------+\n";
	cout << "|                                                |\n";
	cout << "|         PROGRAM :dtmesh(cpp version 040)       |\n";
	cout << "|                                                |\n";
	cout << "|         last update : 2003.12.09               |\n";
	cout << "|         by t.manabe                            |\n";
	cout << "|                                                |\n";
	cout << "+------------------------------------------------+\n";
	cout << "\n";

	dtmesh = new Dtmesh;

	cout << "select mesh type\n";
	cout << "0: cancel\n";
	cout << "1 : triangle(3dtri) to tetrahedra\n";
	cout << "2 : triangle(dli) to tetrahedra\n";
	cout << "3 : triangle(rfi) to tetrahedra\n";
	cout << "4 : point(3dpt) to tetrahedra\n";
	cout << "5 : refinement triangle(3dtri)\n";
	cout << "6 : refinement triangle(dli)\n";
	cout << "7 : refinement triangle(rfi)\n";
	cout << "8 : refinement tetrahedra(3dtet)\n";
	cout << "9 : refinement tetrahedra(dli)\n";
	cout << "10 : refinement tetrahedra(rfi)\n";
	cout << "\n";

	cin >> type;


	//try main routin
	try{

		switch(type) {
		case 1:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputDtm(fname);
			if(flg == FILE_DTM_TRI) {
				density = dtmesh->getPolygon()->getMinlength();
				n = dtmesh->getCrowd()->getSize();
				for(i=0;i<n;i++) {
					dtmesh->getCrowd()->getNode(i)->setDensity(density);
				}
			} else if(flg != FILE_DTM_TRIWD) {
				cout << "cannot open input file\n";
				exit(0);
			}
	
			dtmesh->Mesh();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTetra(fname,"0000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 2:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputDli(fname);

			if(flg == FILE_DLI_TRI) {
				density = dtmesh->getPolygon()->getMinlength();
				n = dtmesh->getCrowd()->getSize();
				for(i=0;i<n;i++) {
					dtmesh->getCrowd()->getNode(i)->setDensity(density);
				}
			} else if(flg != FILE_DLI_TRIWD) {
				cout << "cannot open input file\n";
				exit(0);
			}
	
			dtmesh->Mesh();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTetra(fname,"0000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 3:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputRfi(fname);
			if(flg != FILE_RFI_TRI) {
				cout << "cannot open input file\n";
				exit(0);
			}
			density = dtmesh->getPolygon()->getMinlength();
			dtm::Crowd* cr;
			dtm::Node *nd;
			cr = dtmesh->getCrowd();
			n = cr->getSize();
			for(i=0;i<n;i++) {
				nd = cr->getNode(i);
				nd->setDensity(density);
			}
			dtmesh->Mesh();
	
		/*	dtmesh->meshTriangleToTetra();*/
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTetra(fname,"0000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 4:
			cout << "input file name\n";
			cin >> fname;
			if(!dtmesh->inputDtm(fname)) {
				cout << "cannot open input file\n";
				exit(0);
			}
	
			dtmesh->meshPointToTetra();
		
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTetra(fname,"0000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 5:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputDtm(fname);
			if(flg == FILE_DTM_TRI) {
				density = dtmesh->getPolygon()->getMinlength();
				n = dtmesh->getCrowd()->getSize();
				for(i=0;i<n;i++) {
					dtmesh->getCrowd()->getNode(i)->setDensity(density);
				}
			} else if(flg != FILE_DTM_TRIWD) {
				cout << "cannot open input file\n";
				exit(0);
			}
	
			dtmesh->refineSurface();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTriangle(fname,"00000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 6:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputDli(fname);
			if(flg == FILE_DLI_TRI) {
				density = dtmesh->getPolygon()->getMinlength();
				n = dtmesh->getCrowd()->getSize();
				for(i=0;i<n;i++) {
					dtmesh->getCrowd()->getNode(i)->setDensity(density);
				}
			} else if(flg != FILE_DLI_TRIWD) {
				cout << "cannot open input file\n";
				exit(0);
			}
	
			dtmesh->refineSurface();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTriangle(fname,"00000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 7:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputRfi(fname);
			if(flg != FILE_RFI_TRI) {
				cout << "cannot open input file\n";
				exit(0);
			}
			density = dtmesh->getPolygon()->getMinlength();
			n = dtmesh->getCrowd()->getSize();
			for(i=0;i<n;i++) {
				dtmesh->getCrowd()->getNode(i)->setDensity(density);
			}
	
			dtmesh->refineSurface();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTriangle(fname,"00000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
	
			break;
		case 8:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputDtm(fname);
			if(flg == FILE_DTM_TET) {
				n = dtmesh->getCrowd()->getSize();
				density = dtmesh->getCrowd()->getCurrentMinWidth();
				density = density*2/n;
				for(i=0;i<n;i++) {
					dtmesh->getCrowd()->getNode(i)->setDensity(density);
				}
			} else if(flg != FILE_DTM_TETWD) {
				cout << "cannot open input file\n";
				exit(0);
			}
	
			dtmesh->refineTetgen();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputAvsTetra(fname)) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 9:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputDli(fname);
			if(flg == FILE_DLI_TET) {
				n = dtmesh->getCrowd()->getSize();
				density = dtmesh->getCrowd()->getCurrentMinWidth();
				density = density*2/n;
				for(i=0;i<n;i++) {
					dtmesh->getCrowd()->getNode(i)->setDensity(density);
				}
			} else if(flg != FILE_DLI_TETWD) {
				cout << "cannot open input file\n";
				exit(0);
			}
			dtmesh->refineTetgen();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTetra(fname,"0000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 10:
			cout << "input file name\n";
			cin >> fname;
			flg = dtmesh->inputRfi(fname);
			if(flg != FILE_RFI_TET) {
				cout << "cannot open input file\n";
				exit(0);
			}
			n = dtmesh->getCrowd()->getSize();
			density = dtmesh->getCrowd()->getCurrentMinWidth();
			density = density*2/n;
			for(i=0;i<n;i++) {
				dtmesh->getCrowd()->getNode(i)->setDensity(density);
			}
	
			dtmesh->refineTetgen();
	
			cout << "\n";
			cout << "output file ? yes = 1, no = 0\n";
			cin >> flg;
			if(flg){
				cout << "output file name\n";
				cin >> fname;
				if(!dtmesh->outputRfiTetra(fname,"0000")) {
					cout << "cannot open output file\n";
					exit(0);
				}
			}
			break;
		case 0:
		default:
			break;
		}
	

	}
	catch(dtm::Errors err) {
		cout << "******ERROR******\n";
		cout << err.getFunctionName() << "\n";
		cout << err.getEfect() << "\n";
		cout << "exit program";

		delete dtmesh;

		exit(0);
	}

	delete dtmesh;


	return;
}




/* INPUT FILE FORMAT*/

/*///// *.3dtri /////////////////////////////////////////////

  //number_of_nodes[int]
  n

  //id x y z density_of_node[int,double,double,double,double]
  1    *.**  *.**  *.**  *.**
  ***  *.**  *.**  *.**  *.**
  n-1  *.**  *.**  *.**  *.**

  //number_of_triangle[int]
  m

  //id point1 point2 point3[int,int,int,int]
  1    ***  ***  ***
  ***  ***  ***  ***
  m-1  ***  ***  ***

  /////////////////////////////////////////////////////////*/












//////////////////////////////////////////////////////////////////EOF