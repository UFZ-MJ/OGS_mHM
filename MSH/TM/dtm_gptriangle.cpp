/*

			dtm_gptriangle.h
			>GpTriangle Class

			last update  : 2003.12.08			by t.manabe

  */


#include "stdafx.h" /* MFC */


#include"dtm_gptriangle.h"

namespace dtm{



	//GpTriangle
	///03.12.08

	GpTriangle::GpTriangle()
		:Group(LIMIT_NUMBER_OF_TRIANGLES)
	{

	}


	//~GpTriangle
	///03.12.08

	GpTriangle::~GpTriangle()
	{
		clear();
	}


	//getSize
	///03.12.08

	int GpTriangle::getSize()
	{
		return (int)triangle.size();
	}


	//clear
	///03.12.08

	void GpTriangle::clear()
	{
		int n = (int)triangle.size();
		for(int i=0;i<n;i++) {
			delete triangle[i];
		}
		triangle.clear();
	}
	

	//erase
	///03.12.08
	
	int GpTriangle::erase(int id)
	{
		int edge;
		Triangle *ertri,*swtri,*nbtri;
		int size;

		int n = (int)triangle.size();
		if(id < 1 || id > n) {
			Errors err(0,
				"int GpTriangle::erase(int id)",
				"illegal id");
			throw err;
		}
		
		ertri = triangle[id-1];

		for(int i=0;i<3;i++) {
			int n = ertri->getNumberOfNeighborElement(i);
			for(int j=0;j<n;j++) {
				nbtri = ertri->getNeighborElement(i,j);
				edge = nbtri->getConection(ertri);
				if(edge == -1) continue;
				nbtri->eraseNeighborElement(edge,ertri);
			}
		}

		delete ertri;

		size = (int)triangle.size();
		if(id == size) {
			triangle.pop_back();
		} else {
			swtri = triangle[triangle.size()-1];
			swtri->setId(id);
			triangle[id-1] = swtri;
			triangle.pop_back();
		}

		return (int)triangle.size();
	}


	//getTriangle
	///03.12.08

	Triangle* GpTriangle::getTriangle(int index)
	{
		int n = (int)triangle.size();
		if(index < 0 || index >= n) {
			Errors err(0,
				"Triangle* GpTriangle::getTriangle(int index)",
				"illegal index");
			throw err;
		}

		return triangle[index];
	}


	//makeNewTriangle
	///03.12.08

	Triangle* GpTriangle::makeNewTriangle(Node *nd1,Node *nd2,Node *nd3)
	{
		Triangle* ptri;

		int n = (int)triangle.size();
		if(n >= getMaxSize()) {
			Errors err(1,
				"Triangle* GpTriangle::makeNewTriangle(Node *nd1,Node *nd2,Node *nd3)",
				"number of triangles is max");
			throw err;
		}

		int id = (int)triangle.size()+1;
		ptri = new Triangle(id,nd1,nd2,nd3);
		if(!ptri) {
			Errors err(-1,
				"Triangle* GpTriangle::makeNewTriangle(Node *nd1,Node *nd2,Node *nd3)",
				"failed making new object");
			throw err;
		}
		triangle.push_back(ptri);

		return ptri;
	}



	//view
	///03.12.08

	void GpTriangle::view()
	{
		cout << "POLYGON SIZE = " << (int)triangle.size() << "\n";
		int n = (int)triangle.size();
		for(int i=0;i<n;i++) {
			triangle[i]->view();
		}
	}











}

//////////////////////////////////////////////////////////////////EOF