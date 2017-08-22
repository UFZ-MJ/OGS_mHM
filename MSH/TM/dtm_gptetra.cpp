/*

			dtm_gptetra.h
			>GpTetra Class

			last update : 2003.12.08			by t.manabe

  */


#include "stdafx.h" /* MFC */


#include"dtm_gptetra.h"

namespace dtm{

	//GpTetra
	///03.12.08

	GpTetra::GpTetra()
		:Group(LIMIT_NUMBER_OF_TETRAHEDRAS)
	{


	}


	//~GpTetra
	///03.12.08

	GpTetra::~GpTetra()
	{
		clear();
	}


	//getSize
	///03.12.08

	int GpTetra::getSize()
	{
		return (int)tetra.size();
	}


	//erase
	///03.12.08

	int GpTetra::erase(int id)
	{
		int i,j;
		int size;
		Tetra *crtet,*nbtet,*endtet;

		crtet = tetra[id-1];
		for(i=0;i<4;i++) {
			nbtet = crtet->getNeighborElement(i);
			if(nbtet) {
				j = nbtet->getConection(crtet);
				if(j != -1)nbtet->resetNeighborElement(j,NULL);
			}
		}

		size = (int)tetra.size();
		if(id == size) {
			tetra.pop_back();
			delete crtet;
		} else {
			endtet = tetra[tetra.size()-1];
			tetra[id-1] = endtet;
			endtet->setId(id);
			tetra.pop_back();
			delete crtet;
		}

		return (int)tetra.size();
	}


	//clear
	///03.12.08

	void GpTetra::clear()
	{
		int n = (int)tetra.size();
		for(int i=0;i<n;i++) {
			if(tetra[i]) delete tetra[i];
		}
	}


	//getTetra
	///03.12.08

	Tetra* GpTetra::getTetra(int index)
	{
		int n = (int)tetra.size();
		if(index < 0 || index >= n) {
			Errors err(0,
				"Tetra* Tetgen::getItem(int index)",
				"illegal index");
			throw err;
		}

		return tetra[index];
	}


	//makeNewTetra
	///03.12.08

	Tetra* GpTetra::makeNewTetra(Node *nd1,Node *nd2,Node *nd3,Node *nd4)
	{
		Tetra* ptet;

		int n = (int)tetra.size();
		if(n >= getMaxSize()) {
			Errors err(1,
				"Tetra* GpTetra::makeNewTetra(Node *nd1,Node *nd2,Node *nd3,Node *nd4)",
				"number os tetrahedras is max");
			throw err;
		}

		int id = (int)tetra.size()+1;
		ptet = new Tetra(id,nd1,nd2,nd3,nd4);
		if(!ptet) {
			Errors err(-1,
				"Tetra* Tetgen::makeNewItem(Node *p0,Node *p1,Node *p2,Node *p3)",
				"failed making new object");
			throw err;
		}

		tetra.push_back(ptet);

		return ptet;
	}




}

//////////////////////////////////////////////////////////////////EOF