/*

			dtm_groud.cpp
			>Group Class

			last update : 2003.12.08			by t.manabe

  */


#include "stdafx.h" /* MFC */


#include"dtm_gpnode.h"

namespace dtm{

	//GpNode
	///03.12.08

	GpNode::GpNode()
		:Group(LIMIT_NUMBER_OF_NODES)
	{

	}


	//~GpNode
	///03.12.08

	GpNode::~GpNode()
	{
		clear();
	}


	//getSize
	///03.12.08

	int GpNode::getSize()
	{
		return (int)node.size();
	}


	//erase
	///03.12.08

	int GpNode::erase(int id)
	{
		int n = (int)node.size();
		if(id < 1 || id > n) {
			return -1;
		}

		Node *tmp_nd = node[id-1];
		node[id-1] = node[node.size()-1];
		node[id-1]->setId(id);
		delete tmp_nd;
		node.pop_back();

		return (int)node.size();
	}


	//clear
	///03.12.08

	void GpNode::clear()
	{
		int n = (int)node.size();
		for(int i=0;i<n;i++) {
			delete node[i];
		}
		node.clear();

		return;
	}

	//getNode
	///03.12.08

	Node* GpNode::getNode(int index)
	{
		int n = (int)node.size();
		if(index < 0 || index >= n) {
			Errors err(0,
				"Node* GpNode::getNode(int id)",
				"illegal index");
			throw err;
		}
		return node[index];
	}


	//makeNewNode
	///03.12.08

	Node *GpNode::makeNewNode(double x,double y,double z,double density)
	{
		Node *pnd;

		int n = (int)node.size(); 
		if(n >= getMaxSize()) {
			Errors err(1,
				"Node *GpNode::makeNewNode(double x,double y,double z,double density)",
				"number of nodes is max");
			throw err;
		}

		int id = (int)node.size()+1;
		pnd = new Node(id,x,y,z,density);

		if(!pnd) {
			Errors err(-1,
				"Node *GpNode::makeNewNode(double x,double y,double z,double density)",
				"failed make new object");
			throw err;
		}

		node.push_back(pnd);

		return pnd;
	}


	//view
	///03.12.08

	void GpNode::view()
	{
		cout << (int)node.size() << "\n";

		int n = (int)node.size();
		for(int i=0;i<n;i++) {
			node[i]->view();
		}

		return;
	}




}



//////////////////////////////////////////////////////////////////EOF