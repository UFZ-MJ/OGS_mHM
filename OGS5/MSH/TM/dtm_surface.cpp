/*

			dtm_surface.cpp
			>Surface class

			lst update : 2003.12.08		by t.manabe

*/

#include "stdafx.h" /* MFC */


#include"dtm_surface.h"

#include"dtm_timer.h"
#include"dtm_stdio.h"

namespace dtm{

const static int ij1[3] = {1,2,0};
const static int ij2[3] = {2,0,1};


// Surface
///03.11.23

Surface::Surface()
:flg_lg(0) , cr_crowd(NULL),
laplas_flg(1),rate_density(0.9),
rate_minarea(0.05),stack_maxlong(1000)
{

}


//Surface
///03.12.08

Surface::Surface(Crowd *crowd)
:flg_lg(0) , cr_crowd(crowd),
laplas_flg(1),rate_density(0.9),
rate_minarea(0.05),stack_maxlong(1000)
{

}


// ~Surface
///03.11.23

Surface::~Surface()
{

}


//setCurrentCrowd
//03.12.08
void Surface::setCurrentCrowd(Crowd *crowd)
{
	cr_crowd = crowd;

	return;
}



//initialize
///03.11.23

void Surface::setSurfaceParameter(int lplsflg,double ratednst,
								  double ratearea,int stklong)
{
	laplas_flg = lplsflg;
	rate_density = ratednst;
	rate_minarea = ratearea;
	stack_maxlong = stklong;
}


// setNeighbor
///03.11.23
//determine neighbors of each triangle
//must be setted hanger of triangle

void Surface::setNeighbor()
{

	int i,j,k,l;
	Node *cr_nd1,*cr_nd2;
	Triangle *crtri,*optri;
	int hgsize1,hgsize2;
	int flg;
	int n;
	const int ij[4] = {-1,2,1,0};
	int conect;

	n = getSize();
	for(i=0;i<n;i++) 
		for(j=0;j<3;j++) getTriangle(i)->resetNeighborElement(j,NULL);

	for(i=0;i<n;i++) {
		crtri = getTriangle(i);
		for(j=0;j<3;j++) {
			/////
			//if(crtri->getTopNeighbor(j))continue;
			/////
			cr_nd1 = crtri->getNode((j+1)%3);
			cr_nd2 = crtri->getNode((j+2)%3);
			hgsize1 = cr_nd1->getHangerSize(FG_TYPE_TRI);
			hgsize2 = cr_nd2->getHangerSize(FG_TYPE_TRI);
			flg = 0;
			if(hgsize1 <= hgsize2) {
				for(k=0;k<hgsize1;k++) {
					optri = (Triangle*)cr_nd1->getHanger(k,FG_TYPE_TRI);
					if(*crtri == *optri) continue;
					for(l=0;l<3;l++) {
						if(optri->getNode(l) == cr_nd2) {
							crtri->setNewNeighborElement(j,optri);
							conect = optri->getNodeNumber(cr_nd1);
							optri->setNewNeighborElement(ij[l+conect],crtri);
						//	flg = 1;
							break;
						}
					}
				//	if(flg) break;
				}
			} else {
				for(k=0;k<hgsize2;k++) {
					optri = (Triangle*)cr_nd2->getHanger(k,FG_TYPE_TRI);
					if(*crtri == *optri) continue;
					for(l=0;l<3;l++) {
						if(optri->getNode(l) == cr_nd1) {
							crtri->setNewNeighborElement(j,optri);
							//optri->setNeighbor((l+1)%3,crtri);
							conect = optri->getNodeNumber(cr_nd2);
							optri->setNewNeighborElement(ij[l+conect],crtri);
						//	flg = 1;
							break;
						}
					}
				//	if(flg) break;
				}
			}
		}
	}


	return;
}



// setLength
///03.11.23
//set max and min length

void Surface::setLength()
{
	int i;
	double maxl,minl;
	Triangle *crtri;
	int n;

	if(!getSize()) return;

	crtri = getTriangle(0);
	maxlength = crtri->getMaxlength();
	minlength = crtri->getMinlength();
	maxaverage = maxlength;
	minaverage = minlength;

	n = getSize();
	for(i=1;i<n;i++) {
		crtri = getTriangle(i);
		maxl = crtri->getMaxlength();
		minl = crtri->getMinlength();
		if(maxl > maxlength) maxlength = maxl;
		if(minl < minlength) minlength = minl;
		maxaverage = (maxaverage + maxl)*0.50;
		minaverage = (minaverage + minl)*0.50;
	}

	flg_lg = 1;
}


// getMaxlength
///03.11.23

double Surface::getMaxlength()
{
	if(!flg_lg) setLength();

	return maxlength;
}


// getMinlength
///03.11.23

double Surface::getMinlength()
{
	if(!flg_lg) setLength();

	return minlength;
}


// getMaxaverage
///03.11.23

double Surface::getMaxaverage()
{
	if(!flg_lg) setLength();

	return maxaverage;
}


// getMinaverage
///03.11.23

double Surface::getMinaverage()
{
	if(!flg_lg) setLength();

	return minaverage;
}




// refineOnValue
///03.12.08
///refinement of triangle using triangle value
///value = minl / maxl


void Surface::refineOnValue()
{
	int i,n;
	deque<Triangle*>istack;
	deque<Triangle*>::iterator pp;
	Node *n0,*n1,*n2;
	Node *add_nd;
	Triangle *itri;
	double s,minarea,mindens;
	int ls,lm,le;
	int edge;


	//set minimum area
	mindens = cr_crowd->getCurrentMinDensity();
	minarea = mindens*mindens*rate_minarea;

	//pick up triangle which must be refinement
	while(1) {
		pushStack(istack,minarea);

		if(!istack.size()) break;

		//sort stack using value
		qSortValue(0,(int)istack.size()-1,istack);
		//cut stack over stack_maxlong
		int nn = (int)istack.size();
		while(nn > stack_maxlong) {
			istack.pop_back();
			nn = (int)istack.size();
		}

		n = getSize();
		for(i=0;i<n;i++) {
			itri = getTriangle(i);
			itri->setState(0);
		}

		//refinement triangle
		while(istack.size()) {
			itri = istack[0];
			edge = itri->getEdge();
			istack.pop_front();

			add_nd = createNode(itri,edge);

			//remeshing
			reMesh(add_nd,itri,edge);

			//laplasian
			if(laplas_flg) {
				if(add_nd->getType() == 2) laplasInTri(add_nd);
			}
			
			clearnStack(istack,add_nd);

			n = add_nd->getHangerSize(FG_TYPE_TRI);
			for(i=0;i<n;i++) {
				itri = (Triangle*)(add_nd->getHanger(i,FG_TYPE_TRI));
				if(!itri) continue;

				//itri->setState(0);
				edge = itri->getEdge();
				n0 = itri->getNode((edge+1)%3);
				n1 = itri->getNode((edge+2)%3);
				n2 = itri->getNode(edge);

				if(itri->getMaxlength() <
					(n0->getDensity()+n1->getDensity())*rate_density)continue;
				
				s = getTriangleArea(n2,n0,n1);
				if(s < minarea) continue;

				if(istack.size()) {
					if(itri->getValue() < istack[0]->getValue()) {
						istack.push_front(itri);
					} else if(itri->getValue() 
						< istack[istack.size()-1]->getValue()) {
						istack.push_back(itri);
					} else {
						ls = 0;
						le = (int)istack.size()-1;
						while((le-ls) > 1) {
							lm = (int)((ls+le)*0.5);
							if(itri->getValue() < istack[lm]->getValue()) {
								le = lm;
							} else {
								ls = lm;
							}
						}
						pp = istack.begin()+le;
						istack.insert(pp,itri);
					}
				} else {
					istack.push_back(itri);
				}
			}
			nn = (int)istack.size();
			while(nn > stack_maxlong) {
				istack.pop_back();
				nn = (int)istack.size();
			}
		}
	}

	return;
}


//refineOnRate
///03.12.08
///refinement of triangle using triangle rate
///rate = minl / (maxl * maxl)

void Surface::refineOnRate()
{
	int i,n;
	deque<Triangle*>istack;
	deque<Triangle*>::iterator pp;
	Node *n0,*n1,*n2;
	Node* add_nd;
	Triangle *itri;
	double s,minarea,mindens;
	int ls,lm,le;
	int edge;

	mindens = cr_crowd->getCurrentMinDensity();
	minarea = mindens*mindens*rate_minarea;



	while(1) {

		pushStack(istack,minarea);

		if(!istack.size()) break;
	
		qSortRate(0,(int)istack.size()-1,istack);

		int nn = (int)istack.size();
		while(nn > stack_maxlong) {
			istack.pop_back();
			nn = (int)istack.size();
		}

		n = getSize();
		for(i=0;i<n;i++) {
			getTriangle(i)->setState(0);
		}
		while(istack.size()) {
			itri = istack[0];
			edge = itri->getEdge();
			istack.pop_front();

			add_nd = createNode(itri,edge);

			//remeshing
			reMesh(add_nd,itri,edge);


			//laplasian
			if(laplas_flg) {
				if(add_nd->getType() == 2)	laplasInTri(add_nd);
			}/**/

			clearnStack(istack,add_nd);

			n = add_nd->getHangerSize(FG_TYPE_TRI);
			for(i=0;i<n;i++) {
				itri = (Triangle*)(add_nd->getHanger(i,FG_TYPE_TRI));
				if(!itri) continue;
				itri->setState(0);

				edge = itri->getEdge();
				n0 = itri->getNode((edge+1)%3);
				n1 = itri->getNode((edge+2)%3);
				n2 = itri->getNode(edge);

				if(itri->getMaxlength() <
					(n0->getDensity()+n1->getDensity())*rate_density)continue;

				s = getTriangleArea(n2,n0,n1);
				if(s < minarea) continue;

				if(istack.size()) {
					if(itri->getRate() < istack[0]->getRate()) {
						istack.push_front(itri);
					} else if(itri->getRate() 
						< istack[istack.size()-1]->getRate()) {
						istack.push_back(itri);
					} else {
						ls = 0;
						le = (int)istack.size()-1;
						while((le-ls) > 1) {
							lm = (int)((ls+le)*0.5);
							if(itri->getRate() < istack[lm]->getRate()) {
								le = lm;
							} else {
								ls = lm;
							}
						}
						pp = istack.begin()+le;

						istack.insert(pp,itri);
					}
				} else {
					istack.push_back(itri);
				}
			}
			nn = (int)istack.size();
			while(nn > stack_maxlong) {
				istack.pop_back();
				nn = (int)istack.size();
			}
		}

	}


	return;

}


//refineOnFase
///03.12.08
void Surface::refineOnFast()
{
	int i;
	int n;
	Node *n0,*n1,*n2;
	Node *add_nd;
	Triangle *itri;
	double s;
	double minarea;
	double mindens;
	int edge;
	int flg;

	Point ppt;
	deque<Triangle*>que;


	n = getSize();
	for(i=0;i<n;i++) {
		que.push_back(getTriangle(i));
	}
	mindens = cr_crowd->getCurrentMinDensity();
	minarea = mindens*mindens*rate_minarea*2;

	flg = 0;
	while(que.size()) {
		itri = que[0];
		que.pop_front();

	//	if(!itri->getType()) continue;

		edge = itri->getEdge();
		n0 = itri->getNode((edge+1)%3);
		n1 = itri->getNode((edge+2)%3);
		n2 = itri->getNode(edge);
		if(itri->getMaxlength() < 
			(n0->getDensity()+n1->getDensity())*rate_density)continue;

		s = getTriangleArea(n2,n0,n1);
		if(s < minarea) continue;

		add_nd = createNode(itri,edge);

		//remeshing
		reMesh(add_nd,itri,edge);

		//laplasian
		if(laplas_flg) {
			if(add_nd->getType() == 2) laplasInTri(add_nd);
		}
		n = add_nd->getHangerSize(FG_TYPE_TRI);
		for(i=0;i<n;i++) {
			que.push_back((Triangle*)(add_nd->getHanger(i,FG_TYPE_TRI)));
		}
	}

	return;
}


//pushStack
///03.12.08

int Surface::pushStack(deque<Triangle*>&istack,double minarea)
{
	Node *n0,*n1,*n2;
	Triangle *crtri;
	int edge;
	double area;

	int n = getSize();
	for(int i=0;i<n;i++) {
		crtri = getTriangle(i);
		//if(!crtri->getType()) continue;

		edge = crtri->getEdge();
		n0 = crtri->getNode((edge+1)%3);
		n1 = crtri->getNode((edge+2)%3);
		n2 = crtri->getNode(edge);
		if(crtri->getMaxlength() < 
			(n0->getDensity()+n1->getDensity())*rate_density)continue;
		
		area = getTriangleArea(n2,n0,n1);
		if(area < minarea) continue;

		istack.push_back(crtri);
	}

	return (int)istack.size();
}


//clearnStack
///03.12.08

int Surface::clearnStack(deque<Triangle*>&istack,Node *add_nd)
{
	int i,n;
	Triangle *crtri;
	deque<Triangle*>::iterator itr;

	n = add_nd->getHangerSize(FG_TYPE_TRI);
	for(i=0;i<n;i++) {
		crtri = (Triangle*)add_nd->getHanger(i,FG_TYPE_TRI);
		crtri->setState(1);
	}
	n = (int)istack.size();
	for(i=0;i<n;i++) {
		crtri = istack[i];
		if(crtri->getState() == 1) {
			itr = istack.begin()+i;
			istack.erase(itr);
			i--;
			n--;
			continue;
		}
	}
	n = add_nd->getHangerSize(FG_TYPE_TRI);
	for(i=0;i<n;i++) {
		crtri = (Triangle*)add_nd->getHanger(i,FG_TYPE_TRI);
		crtri->setState(0);
	}


	return (int)istack.size();
}


//qSortValue
///03.11.23
//quick sorting
//used only fineValue

void Surface::qSortValue(int f,int e,deque<Triangle*>&array)
{
	int i,k;
	double pv,key;
	Triangle* tmptri;
	int flg;

	if(f >= e) return;

	flg = 1; 
	pv = array[f]->getValue();
	for(i=f+1;i<=e;i++) {
		key = array[i]->getValue();
		if(key > pv) {
			pv = key;
			flg = 0;
			break;
		} else if(key < pv) {
			flg = 0;
			break;
		}
	}
	if(flg)return;

	k = e;

	for(i=f;i<=k;i++) {
		key = array[i]->getValue();
		if(key >= pv) {
			tmptri = array[i];
			array[i] = array[k];
			array[k] = tmptri;
			k--;
			i--;
		}
	}

	qSortValue(f,i-1,array);
	qSortValue(i,e,array);

	return;
}


//qSortRate
///03.11.23
//quick sorting
//used only fineRate

void Surface::qSortRate(int f,int e,deque<Triangle*>&array)
{
	int i,k;
	double pv,key;
	Triangle* tmptri;
	int flg;

	if(f >= e) return;

	flg = 1; 
	pv = array[f]->getRate();
	for(i=f+1;i<=e;i++) {
		key = array[i]->getRate();
		if(key > pv) {
			pv = key;
			flg = 0;
			break;
		} else if(key < pv) {
			flg = 0;
			break;
		}
	}
	if(flg)return;

	k = e;

	for(i=f;i<=k;i++) {
		key = array[i]->getRate();
		if(key >= pv) {
			tmptri = array[i];
			array[i] = array[k];
			array[k] = tmptri;
			k--;
			i--;
		}
	}

	qSortRate(f,i-1,array);
	qSortRate(i,e,array);

	return;
}


//createNode
///03.12.08

Node* Surface::createNode(Triangle *tri,int edge)
{
	Node *add_nd;
	Node *n0,*n1;
	double ix,iy,iz,iden;
	Triangle *ptri;

	n0 = tri->getNode((edge+1)%3);
	n1 = tri->getNode((edge+2)%3);
	ix = (n0->getX() + n1->getX())*0.5;
	iy = (n0->getY() + n1->getY())*0.5;
	iz = (n0->getZ() + n1->getZ())*0.5;
	iden = (n0->getDensity() + n1->getDensity())*0.5;

	add_nd = cr_crowd->makeNewNode(ix,iy,iz,iden);

	int num = tri->getNumberOfNeighborElement(edge);
	if(num == 1) {
		ptri = tri->getNeighborElement(edge,0);
		if(tri->getType() != ptri->getType()) {
			add_nd->setType(1);		//this node is on border of surface(edge)
		} else if(tri->getAngle(ptri) < FRAT_ANGLE ){
			add_nd->setType(1);		//this node is on edge
		} else {
			add_nd->setType(2);		//this node is on surface(not edge) 
		}
	} else {
		add_nd->setType(1);
	}

	return add_nd;
}


//reMesh
///03.12.08

void Surface::reMesh(Node *add_nd,Triangle *f_cr_tri,int edge)
{
	int i,nn;
	Triangle *f_nw_tri; //new triangle created in f_cr_tri
	Triangle *nb_tri;
	deque<Triangle*>ntri;
	deque<Triangle*>istack;
	Node *tp_nd;	// top node
	Node *rt_nd;	// right node
	Node *lt_nd;	// left node
	int nb_edge;
	int type;
	double density;
	double minarea;
	int domain;
	Neighbor neib;	// neighbor triangles which must be remeshed(=op_nb) 
	Neighbor op_nb;	// opposite side neighor(opposite top node)
	Neighbor ls_nb;	// left side neighor(opposite right node)
	Neighbor rs_nb;	// right side neighor(opposite left node)


	density = add_nd->getDensity();
	minarea = density*density*rate_minarea*0.5;

	neib = f_cr_tri->getNeighbor(edge);

	tp_nd = f_cr_tri->getNode(edge);
	rt_nd = f_cr_tri->getNode((edge+1)%3);
	lt_nd = f_cr_tri->getNode((edge+2)%3);
	op_nb = f_cr_tri->getNeighbor(edge);
	ls_nb = f_cr_tri->getNeighbor((edge+1)%3);
	rs_nb = f_cr_tri->getNeighbor((edge+2)%3);

	//renew element f_cr_tri
	type = f_cr_tri->getType();
	domain = f_cr_tri->getDomain();
	f_cr_tri->resetNode(add_nd,tp_nd,rt_nd);
	f_cr_tri->resetNeighbor(0,rs_nb);
	f_cr_tri->resetNeighbor(1,op_nb);
	f_cr_tri->setType(type);
	f_cr_tri->setDomain(domain);

	//make new element nelm3+2
	f_nw_tri = makeNewTriangle(add_nd,lt_nd,tp_nd);
	f_nw_tri->resetNeighbor(0,ls_nb);
	f_nw_tri->setNewNeighborElement(1,f_cr_tri);
	f_nw_tri->setType(type);
	f_nw_tri->setDomain(domain);

	f_cr_tri->resetNeighborElement(2,f_nw_tri);

	nn = ls_nb.getNumberOfElement();
	for(i=0;i<nn;i++) {
		nb_tri = (Triangle*)ls_nb.getElement(i);
		nb_edge = nb_tri->getConection(f_cr_tri);
		nb_tri->eraseNeighborElement(nb_edge,f_cr_tri);
		nb_tri->setNewNeighborElement(nb_edge,f_nw_tri);
	}

	if(ls_nb.getNumberOfElement() == 1) {
		//push stack for swap
		istack.push_back(f_nw_tri);
	} 
	if(rs_nb.getNumberOfElement() == 1) {
		//push stack for swap
		istack.push_back(f_cr_tri);
	} 

	remeshOppositeSide(neib,f_cr_tri,f_nw_tri,tp_nd,rt_nd,lt_nd,istack,add_nd);

	swapTriangle(istack,add_nd,minarea);


	return;
}


//remeshOppositeSide
///03.12.08

int Surface::remeshOppositeSide(Neighbor &neib,
								Triangle *f_cr_tri,Triangle *f_nw_tri,
								Node *tp_nd,Node *rt_nd,Node *lt_nd,
								deque<Triangle*>&istack,
								Node *add_nd)
{
	deque<Triangle*>nb_tri_buff;
	Triangle *op_cr_tri; //neighbor triangle of f_cr_tri
	Triangle *op_nw_tri; //new triangle created in op_cr_tri
	int op_tp_edge,op_rt_edge,op_lt_edge;
	Node *op_tp_nd;		// neighbor`s top node
	Neighbor op_op_nb;	//	neighbor`s opposite side neighbor
	Neighbor op_rs_nb;	//	neighbor`s right side neighbor
	Neighbor op_ls_nb;	//	neighbor`s left side neighbor
	int type,domain;
	Triangle *nb_tri;	//triangle of one of neighbors
	Triangle *nn_tri;	//triangle of neighbor`s neighbor
	int nnedge;
	long tkdummy;

	if (tp_nd == NULL) tkdummy=0;

	int num = neib.getNumberOfElement();
	for(int i=0;i<num;i++) {
		op_cr_tri = neib.getElement(i);
		op_tp_edge = op_cr_tri->getConection(f_cr_tri);
		op_rt_edge = op_cr_tri->getNodeNumber(rt_nd);
		op_lt_edge = op_cr_tri->getNodeNumber(lt_nd);
		op_tp_nd = op_cr_tri->getNode(op_tp_edge);

		op_op_nb = op_cr_tri->getNeighbor(op_tp_edge);
		op_rs_nb = op_cr_tri->getNeighbor(op_rt_edge);
		op_ls_nb = op_cr_tri->getNeighbor(op_lt_edge);

		//renew opposite tririangle
		type = op_cr_tri->getType();
		domain = op_cr_tri->getDomain();
		op_cr_tri->resetNode(add_nd,rt_nd,op_tp_nd);
		op_cr_tri->resetNeighbor(0,op_ls_nb);
		op_cr_tri->resetNeighbor(2,op_op_nb);
		op_cr_tri->setType(type);
		op_cr_tri->setDomain(domain);

		//make new triangle in opposite
		op_nw_tri = makeNewTriangle(add_nd,op_tp_nd,lt_nd);
		op_nw_tri->resetNeighbor(0,op_rs_nb);
		op_nw_tri->resetNeighborElement(2,op_cr_tri);
		op_nw_tri->setType(type);
		op_nw_tri->setDomain(domain);

		//conect new triangle
		op_cr_tri->resetNeighborElement(1,op_nw_tri);
		f_nw_tri->setNewNeighborElement(2,op_nw_tri);
		op_nw_tri->setNewNeighborElement(1,f_nw_tri);

		int n = (int)nb_tri_buff.size();
		for(int j=0;j<n;j++) {
			nb_tri = nb_tri_buff[j];
			nb_tri->setNewNeighborElement(1,op_nw_tri);
			op_nw_tri->setNewNeighborElement(1,nb_tri);
		}
		nb_tri_buff.push_back(op_nw_tri);

		int nn = op_rs_nb.getNumberOfElement();
		for(int j=0;j<nn;j++) {
			nn_tri = op_rs_nb.getElement(j);
			nnedge = nn_tri->getConection(op_cr_tri);
			nn_tri->eraseNeighborElement(nnedge,op_cr_tri);
			nn_tri->setNewNeighborElement(nnedge,op_nw_tri);
		}
		if(op_rs_nb.getNumberOfElement() == 1) {
			//push stack for swap
			istack.push_back(op_nw_tri);
		} 
		if(op_ls_nb.getNumberOfElement() == 1) {
			//push stack for swap
			istack.push_back(op_cr_tri);
		} 
	}

	return 0;
}


//swapTriangle
///03.12.08

int Surface::swapTriangle(deque<Triangle*>&istack,Node *add_nd,double minarea)
{
	int i,nn;
	Triangle *lt_tri;	// this side triangle(left triangle after swap)
	Triangle *rt_tri;	// opposite side triangle(right triangle after swap)
	Triangle *nb_tri;
	int op_tp_edge;		// top node of opposite triangle
	int op_lt_edge;		// left node of opposite triangle
	int op_rt_edge;		// right node of opposite triangle
	int th_lt_edge;		// left node of this side
	int th_rt_edge;		// right node of this side
	int nb_edge;
	Neighbor op_ls_nb;
	Neighbor op_rs_nb;
	Neighbor th_ls_nb;
	Neighbor th_rs_nb;
	int type,domain;
	Node *tp_nd,*lt_nd,*rt_nd;
	double hori_dist,vert_dist;
	double angl;
	double prev_area1,prev_area2,new_area1,new_area2;

	while(istack.size()) {
		//pop triangle from stack
		lt_tri = istack[istack.size()-1];
		rt_tri = lt_tri->getNeighborElement(0,0);
		istack.pop_back();

		op_tp_edge = rt_tri->getConection(lt_tri);
		op_lt_edge = (op_tp_edge+1)%3;
		op_rt_edge = (op_tp_edge+2)%3;
		tp_nd = rt_tri->getNode(op_tp_edge);
		lt_nd = rt_tri->getNode(op_lt_edge);
		rt_nd = rt_tri->getNode(op_rt_edge);

		//compare horizontal line with vertical line
		hori_dist = getDistanceSquare(lt_nd,rt_nd);
		vert_dist = getDistanceSquare(tp_nd,add_nd);
		if(hori_dist < vert_dist) continue;

		//check angle if flat
		angl = lt_tri->getAngle(rt_tri);
		if(angl < FRAT_ANGLE) continue;

		//check area if same
		prev_area1 = getTriangleArea(add_nd,rt_nd,lt_nd);
		prev_area2 = getTriangleArea(tp_nd,lt_nd,rt_nd);
		new_area1 = getTriangleArea(add_nd,tp_nd,lt_nd);
		new_area2 = getTriangleArea(add_nd,rt_nd,tp_nd);
		if(new_area1 < minarea || new_area2 < minarea)continue;
		if( fabs((new_area1+new_area2)-(prev_area1+prev_area2)) > LIMIT )continue;

		//preparation of swapping
		th_lt_edge = lt_tri->getNodeNumber(lt_nd);
		th_rt_edge = lt_tri->getNodeNumber(rt_nd);

		th_ls_nb = lt_tri->getNeighbor(th_lt_edge);
		th_rs_nb = lt_tri->getNeighbor(th_rt_edge);
		op_ls_nb = rt_tri->getNeighbor(op_lt_edge);
		op_rs_nb = rt_tri->getNeighbor(op_rt_edge);

		//swap lt_tri(this side to left)
		type = lt_tri->getType();
		domain = lt_tri->getDomain();
		lt_tri->resetNode(add_nd,tp_nd,lt_nd);
		lt_tri->resetNeighbor(0,op_rs_nb);
		lt_tri->resetNeighbor(1,th_rs_nb);
		lt_tri->resetNeighborElement(2,rt_tri);
		lt_tri->setType(type);
		lt_tri->setDomain(domain);

		//swap rt_tri(opposit to right)
		type = rt_tri->getType();
		domain = rt_tri->getDomain();
		rt_tri->resetNode(add_nd,rt_nd,tp_nd);
		rt_tri->resetNeighbor(0,op_ls_nb);
		rt_tri->resetNeighborElement(1,lt_tri);
		rt_tri->resetNeighbor(2,th_ls_nb);
		rt_tri->setType(type);
		rt_tri->setDomain(domain);


		//re-conect neighbors
		nn = op_rs_nb.getNumberOfElement();
		for(i=0;i<nn;i++) {
			nb_tri = op_rs_nb.getElement(i);
			nb_edge = nb_tri->getConection(rt_tri);
			nb_tri->eraseNeighborElement(nb_edge,rt_tri);
			nb_tri->setNewNeighborElement(nb_edge,lt_tri);
		}
		nn = th_ls_nb.getNumberOfElement();
		for(i=0;i<nn;i++) {
			nb_tri = th_ls_nb.getElement(i);
			nb_edge = nb_tri->getConection(lt_tri);
			nb_tri->eraseNeighborElement(nb_edge,lt_tri);
			nb_tri->setNewNeighborElement(nb_edge,rt_tri);
		}
		if(op_ls_nb.getNumberOfElement() == 1) {
			//push stack next triangle
			istack.push_back(rt_tri);
		}
		if(op_rs_nb.getNumberOfElement() == 1) {
			//push stack next triangle
			istack.push_back(lt_tri);
		}
	}
		
	return 0;
}


//pickBoundary
///03.12.08

void Surface::pickBoundary()
{
	int i,n;

	n = getSize();
	for(i=0;i<n;i++) {
		if(getTriangle(i)->getType() == 1) {
			getTriangle(i)->setVisible();
		} else {
			getTriangle(i)->setInvisible();
		}
	}

	return;
}


//pickDomain
///03.12.08

void Surface::pickDomain(vector<int>&area)
{
	int n = getSize();
	int nn = (int)area.size();

	for(int i=0;i<n;i++) {
		getTriangle(i)->setInvisible();
		for(int j=0;j<nn;j++) {
			if(getTriangle(i)->getDomain() == area[j]) {
				getTriangle(i)->setVisible();
				break;
			}
		}
	}


	return;
}


//sortDirection
///03.11.27
///check suraface direction using ray trace method

void Surface::sortDirection()
{
	int i;
	Triangle *hit,*itri;
	Point hitpt,ipt;
	double iarea,parea,aarea,barea,carea;
	double dist,ans,dans,a,b,c,d;
	double ix,iy,iz;
	double xm,ym,zm;
	deque<double>stkdist;
	deque<Triangle*>stktri;
	int flg = 0;
	Point lmtpt(-1000.,-1000.,-1000.);	//start point of ray
	int num;


	///cr_crowd->normalize();

	//divide the section
	divideSection();

	//hit the ray to represent point of each section
	int n = getSize();
	while(hittri.size()) {
		hit = hittri[0];
		hittri.pop_front();
		hitpt = getTriangleCenter(hit);
		xm = hitpt.getX() - lmtpt.getX();
		ym = hitpt.getY() - lmtpt.getY();
		zm = hitpt.getZ() - lmtpt.getZ();
		stkdist.clear();
		stktri.clear();
		for(i=0;i<n;i++) {
			itri = getTriangle(i);
			if(!itri->isVisible()) continue;
			a = itri->getPlaneVector(0);
			b = itri->getPlaneVector(1);
			c = itri->getPlaneVector(2);
			d = itri->getPlaneVector(3);
			ans = a*xm + b*ym + c*zm;
			dans = itri->getPlaneAnswer(&lmtpt);

			if(ans < LIMIT && ans > -LIMIT)continue;
			dist = -1.*(dans/ans);
			ix = xm*dist + lmtpt.getX();
			iy = ym*dist + lmtpt.getY();
			iz = zm*dist + lmtpt.getZ();
			ipt.setCoordinates(ix,iy,iz);

			parea = itri->getArea();
			aarea = getTriangleArea(&ipt,itri->getNode(1),itri->getNode(2));
			barea = getTriangleArea(&ipt,itri->getNode(2),itri->getNode(0));
			carea = getTriangleArea(&ipt,itri->getNode(0),itri->getNode(1));
			iarea = aarea + barea + carea;

			if(fabs(parea - iarea) > LIMIT)continue;

			stkdist.push_back(dist);
			stktri.push_back(itri);
			itri->setState(2);
		}

		//sort to small to big
		if(stkdist.size() > 1) {
			qSortDist(0,(int)stkdist.size()-1,stkdist,stktri);
		}

		checkSamePlane(stktri,stkdist,lmtpt);

		//check direction
		num = (int)stktri.size();
		for(i=0;i<num;i++) {
			itri = stktri[i];
			dans = itri->getPlaneAnswer(&lmtpt);
			if(i%2) {
				if(dans < 0) {
					itri->setState(1);
				}else {
					flg = 1;
					itri->setState(-1);
				}
			} else {
				if(dans < 0) {
					flg = 1;
					itri->setState(-1);
				}else {
					itri->setState(1);
				}
			}
		}
	}

	if(flg) {
		setDirection();
	}else {
		n = getSize();
		for(i=0;i<n;i++) getTriangle(i)->setState(0);
	}


	//cr_crowd->restore();

	return;
}


//checkSamePlane
///03.12.08

int Surface::checkSamePlane(deque<Triangle*>&stktri,deque<double>&stkdist,
							Point &lmtpt)
{
	Triangle *itri,*ptri;
	double dans,dd;

	//check the error of number
	int num = (int)stktri.size();
	for(int i=0;i<num-1;i++) {
		itri = stktri[i];
		dans = itri->getPlaneAnswer(&lmtpt);
		for(int j=0;j<3;j++) {
			int nn = itri->getNumberOfNeighborElement(j);
			for(int l=0;l<nn;l++) {
				ptri = itri->getNeighborElement(i,j);
				if(!ptri)continue;
				if(ptri->getState() != 2)continue;
				if(!ptri->isVisible())continue;
				dd = ptri->getPlaneAnswer(&lmtpt);
				if(dans*dd < LIMIT)continue;
				for(int k=i+1;k<num;k++) {
					if(ptri == stktri[k]) {
						deque<double>::iterator dp = stkdist.begin() + k;
						deque<Triangle*>::iterator dt = stktri.begin() + k;
						ptri->setState(0);
						stkdist.erase(dp);
						stktri.erase(dt);
						num = (int)stktri.size();
						break;
					}
				}
			}
		}
	}

	return (int)stktri.size();
}



//divideSection
///03.11.27
///divide the suraface to each section(non conection suraface)
///arrange the suraface direction in same direction

void Surface::divideSection()
{
	int i;
	Triangle *itri;
	deque<Triangle*>que;
	int flg;
	int nsec;
	int start;


	start = 0;
	hittri.clear();
	if(!getSize()) return;
	int n = getSize();
	for(i=0;i<n;i++) { 
		if(!(getTriangle(i)->isVisible())) continue;
		start = i;
		break;
	}

	nsec =1;
	itri = getTriangle(start);
	itri->setState(1);
	que.push_back(itri);
	hittri.push_back(itri);

	flg = 1;
	while(flg) {
		flg = 0;

		correctTriangleDirection(que);

		n = getSize();
		for(i=0;i<n;i++) {
			if(!(getTriangle(i)->isVisible())) continue;
			if(!(getTriangle(i)->getState())) {
				nsec++;
				que.push_back(getTriangle(i));
				hittri.push_back(getTriangle(i));
				flg = 1;
				break;
			}
		}
	}

	n = getSize();
	for(i=0;i<n;i++) getTriangle(i)->setState(0);


	return;
}


//setDirection
///03.11.27
///correct the surface direction

void Surface::setDirection()
{
	int i;
	Triangle *itri;
	deque<Triangle*>que;

	int n = getSize();
	for(i=0;i<n;i++) {
		itri = getTriangle(i);
		if(itri->getState() == -1) {
			itri->reverseDirection();
			itri->setState(1);
			que.push_back(itri);
			continue;
		} else if(itri->getState() == 1) {
			que.push_back(itri);
			continue;
		}
		itri->setState(0);
	}

	correctTriangleDirection(que);

	n = getSize();
	for(i=0;i<n;i++) getTriangle(i)->setState(0);

	return;
}


//correctTriangleDirection
///03.12.08

int Surface::correctTriangleDirection(deque<Triangle*>&que)
{
	Triangle *itri,*jtri;
	int conect;
	Node *in1,*in2,*jn1,*jn2;

	while(que.size()) {
		itri = que[0];
		que.pop_front();
		for(int i=0;i<3;i++) {
			int num = itri->getNumberOfNeighborElement(i);
			for(int j=0;j<num;j++) {
				jtri = itri->getNeighborElement(i,j);
				if(!jtri) continue;
				if(!jtri->isVisible()) continue;
				if(jtri->getState()) continue;

				jtri->setState(1);
				que.push_back(jtri);

				conect = jtri->getConection(itri);
				in1 = itri->getNode((i+1)%3);
				in2 = itri->getNode((i+2)%3);
				jn1 = jtri->getNode((conect+1)%3);
				jn2 = jtri->getNode((conect+2)%3);
				if( (*in1 == *jn2) && (*in2 == *jn1) ) {
					continue;
				} else if( (*in1 == *jn1) && (*in2 == *jn2) ){
					jtri->reverseDirection();
				} else {
					Errors err(0,
						"void Surface::correctDirection()",
						"conection is false");
					throw err;
				}
			}
		}
	}

	return 0;
}


//selectDomain
///03.12.08
///set triangle state 
///in domain = 1, out domain = -1 , other domain = 0

void Surface::selectDomain(Tetgen *tetgen)
{
	int flg;
	int judge;
	Tetra *itet;
	Point pt;
	Triangle itri;

	int n = tetgen->getSize();
	for(int i=0;i<n;i++) {

		flg = 0;
		itet = tetgen->getTetra(i);
		focus.clear();
		center.clear();

		flg = pickTriangle(itet);
		if(!flg) {
			flg = pickLine(itet);
		}
		if(!flg) {
			flg = pickPoint(itet);
		}
		if(!flg) {
			itet->setState(0);
			continue;
		}
		int nn = (int)focus.size();
		for(int j=0;j<nn;j++) {
			pt = getTriangleCenter(focus[j]);
			center.push_back(pt);
		}

		judge = judgeInOut(itet);

		if(judge == 0) {
			itet->setState(-1);
		} else {
			itet->setState(1);
		}

		nn = (int)focus.size();
		for(int j=0;j<nn;j++) {
			focus[j]->setState(0);
		}
	}

	return;
}


//pickTriangle
///03.12.08

int Surface::pickTriangle(Tetra *tet)
{

	int i=0;
	Node *p0,*p1,*p2;
	Triangle *t0,*t1,*t2;
	int i0=0,i1=0,i2=0;
	int id0=0,id1=0,id2=0;
	int ct[4][3] = {1,3,2,2,3,0,3,1,0,0,1,2};
	int n0=0,n1=0,n2=0;


	for(i=0;i<4;i++) {
		p0 = tet->getNode(ct[i][0]);
		p1 = tet->getNode(ct[i][1]);
		p2 = tet->getNode(ct[i][2]);

		n0 = p0->getHangerSize(FG_TYPE_TRI);
		n1 = p1->getHangerSize(FG_TYPE_TRI);
		n2 = p2->getHangerSize(FG_TYPE_TRI);


		i0 = i1 = i2 = 0;
		while(1) {
			if(i0 >= n0)break;
			if(i1 >= n1)break;
			if(i2 >= n2)break;
			t0 = (Triangle*)p0->getHanger(i0,FG_TYPE_TRI);
			if(!t0->isVisible()) {
				i0++;
				continue;
			}
			id0 = t0->getId();
			while(i1 != n1) {
				t1 = (Triangle*)p1->getHanger(i1,FG_TYPE_TRI);
				if(!t1->isVisible()) {
					i1++;
					continue;
				}
				id1 = t1->getId();
				if(id1 >= id0) break;
				i1++;
			}
			if(i1 == n1)break;

			while(i2 != n2) {
				t2 = (Triangle*)p2->getHanger(i2,FG_TYPE_TRI);
				if(!t2->isVisible()) {
					i2++;
					continue;
				}
				id2 = t2->getId();
				if(id2 >= id1) break;
				i2++;
			}
			if(i2 == n2)break;
			while(i0 != n0) {
				t0 = (Triangle*)p0->getHanger(i0,FG_TYPE_TRI);
				if(!t0->isVisible()) {
					i0++;
					continue;
				}
				id0 = t0->getId();
				if(id0 >= id2) break;
				i0++;
			}
			if(i0 == n0)break;

			if( (id0 == id1) && (id0 == id2) ){
				focus.push_back((Triangle*)p0->getHanger(i0,FG_TYPE_TRI));
			//	i0++;
			//	i1++;
			//	i2++;
				break;
			}
		}
	}


	return (int)focus.size();

}


// pickLine
///03.11.23

int Surface::pickLine(Tetra *tet)
{
	int i;
	Node *p0,*p1;
	Triangle *t0,*t1;
	int i0,i1;
	int id0=0,id1=0;
	int n0,n1;
	const int j0[] = {0,0,0,1,1,2};
	const int j1[] = {1,2,3,2,3,3};

	for(i=0;i<6;i++) {
		p0 = tet->getNode(j0[i]);
		p1 = tet->getNode(j1[i]);
		n0 = p0->getHangerSize(FG_TYPE_TRI);
		n1 = p1->getHangerSize(FG_TYPE_TRI);
		i0 = i1 = 0;
		while(1) {
			if(i0 >= n0)break;
			if(i1 >= n1)break;
			t0 = (Triangle*)p0->getHanger(i0,FG_TYPE_TRI);
			if(!t0->isVisible()) {
				i0++;
				continue;
			}
			id0 = t0->getId();
			while(i1 != n1) {
				t1 = (Triangle*)p1->getHanger(i1,FG_TYPE_TRI);
				if(!t1->isVisible()) {
					i1++;
					continue;
				}
				id1 = t1->getId();
				if(id1 >= id0) break;
				i1++;
			}
			if(i1 == n1)break;
			while(i0 != n0) {
				t0 = (Triangle*)p0->getHanger(i0,FG_TYPE_TRI);
				if(!t0->isVisible()) {
					i0++;
					continue;
				}
				id0 = t0->getId();
				if(id0 >= id1) break;
				i0++;
			}
			if(i0 == n0)break;
			if(id0 == id1) {
				focus.push_back(t0);
				i0++;
				i1++;
			}
		}
	}

	return (int)focus.size();
}


// pickPoint
///03.11.23

int Surface::pickPoint(Tetra *tet)
{
	int i,j;
	Node *p0;
	int n;
	Triangle *tri;

	for(i=0;i<4;i++) {
		p0 = tet->getNode(i);
		n = p0->getHangerSize(FG_TYPE_TRI);
		for(j=0;j<n;j++){
			tri = (Triangle*)p0->getHanger(j,FG_TYPE_TRI);
			if(!tri->isVisible()) continue;
			focus.push_back(tri);
		}
	}

	return (int)focus.size();

}

// judgeInOut
///03.11.23
/// return 1 : inside , return 0 : outside

int Surface::judgeInOut(Tetra *tet)
{
	int i,j;
	int flgin,flgout;
	Point cp;
	Triangle *itri,*jtri;
	Plane plane;
	double ans;
	double ave;
	double param;
	int fnum;


	fnum =(int)focus.size();
	ave = 0.;
	for(i=0;i<fnum;i++) {
		ave += focus[i]->getMinlength();
	}
	ave = ave / fnum;
	param = ave*LIMIT_10;

	flgin = flgout = 0;
	cp = getTetraCenter(tet);
	for(i=0;i<fnum;i++) {
		itri = focus[i];
		itri->setState(0);
		plane = itri->getPlane();
		ans = plane.getAnswer(&cp);

		if(ans > param) {
			flgout = 1;
			itri->setState(-1);
		} else {
			flgin = 1;
			itri->setState(1);
		}
	}
	
	if(flgout && !flgin) return 0;
	if(!flgout && flgin) return 1;

	for(i=0;i<fnum;i++) {
		itri = focus[i];
		plane = itri->getPlane();
		if(itri->getState() != 1) continue;
		for(j=0;j<fnum;j++) {
			jtri = focus[j];
			if(jtri->getState() != -1) continue;
			cp = center[j];
			ans = plane.getAnswer(&cp);
			if(ans > param) {
				jtri->setState(0);
			}
		}
	}
	for(i=0;i<fnum;i++) if(focus[i]->getState() == -1) return 0;
	if(!focus.size()) return 0;

	return 1;
}


//cleanBadTriangle
///03.11.23

void Surface::cleanBadTriangle(double min_value)
{
	int edge;
	Node *nd;
	Triangle *itri;
	double mindens,minarea;


	mindens = cr_crowd->getCurrentMinDensity();
	minarea = mindens*mindens*rate_minarea;
	
	int n = getSize();
	for(int i=0;i<n;i++) {
		itri = getTriangle(i);

		if(itri->getArea() > minarea &&
			itri->getValue() > min_value) continue;
		edge = itri->getEdge();
		for(int j=0;j<3;j++) {
			nd = itri->getNode((edge+j)%3);
			if(nd->getType() == 0 )continue;
			if(eraseNode(nd))break;
		}
		n = getSize();
	}

	return;
}


//eraseNode
///03.11.23
///remove node from triangle

Node* Surface::eraseNode(Node* er_nd)
{
	int i,j,k;
	int n,nn;
	Triangle *itri;
	Node *pnd,*tmpnd;
	double lg;
	deque<double>lg_buff;
	deque<Node*>po_nd_buff;
	deque<double>::iterator it_d;
	deque<Node*>::iterator it_n;
	int flg;
	int type;
	int tp_flg = 0;

	type = er_nd->getType();
	if(type == 2) {
		tp_flg = 0;
	} else if(type == 1) {
		tp_flg = 1;
	} else {
		return NULL;
	}

	n = er_nd->getHangerSize(FG_TYPE_TRI);
	for(i=0;i<n;i++) {
		itri = (Triangle*)er_nd->getHanger(i,FG_TYPE_TRI);
		for(j=0;j<3;j++) {
			pnd = itri->getNode(j);
			if(*pnd == *er_nd)continue;
			if(tp_flg) {
				if(pnd->getType() == 2)continue;
			}
			flg = 0;
			nn = (int)po_nd_buff.size();
			for(k=0;k<nn;k++) {
				tmpnd = po_nd_buff[k];
				if(*tmpnd == *pnd) {
					flg = 1;
					break;
				}
			}
			if(flg)continue;
			lg = getDistanceSquare(pnd,er_nd);
			for(k=0;k<nn;k++) {
				if(lg < lg_buff[k])break;
			}
			if(k == 0) {
				po_nd_buff.push_front(pnd);
				lg_buff.push_front(lg);
			} else if(k == nn){
				po_nd_buff.push_back(pnd);
				lg_buff.push_back(lg);
			} else {
				it_n = po_nd_buff.begin() + k;
				po_nd_buff.insert(it_n,pnd);
				it_d = lg_buff.begin() + k;
				lg_buff.insert(it_d,lg);
			}
		}
	}

	//try merge in order of length
	n = (int)po_nd_buff.size();
	for(i=0;i<n;i++) {
		pnd = po_nd_buff[i];
		if(!mergeNode(pnd,er_nd)) continue;
		return pnd;
	}

	return NULL;

}


//mergeNode
///03.12.10
///merge node ne_nd to po_nd(erase ne_nd)
///if succeeded return po_nd,if not return NULL

Node* Surface::mergeNode(Node *po_nd,Node *ne_nd)
{
	int i,j,k,l;
	Triangle *potri = NULL,*netri = NULL,*nbtri = NULL;
	deque<Triangle*>pohang;
	deque<Triangle*>nehang;
	Neighbor neib,nw_nb;
	int edge,ndedge;
	int flg;

	//set pohang ang nehang
	//triangle of pohang state = 1
	//triangle of nehang state = -1
	divHanger(po_nd,ne_nd,pohang,nehang);

	//check area, if can merge or not
	if(!isPossibleMerge(po_nd,ne_nd,pohang,nehang)) {
		int nn = (int)pohang.size();
		for(i=0;i<nn;i++) {
			pohang[i]->setState(0);
		}
		nn = (int)nehang.size();
		for(i=0;i<nn;i++) {
			nehang[i]->setState(0);
		}
		return NULL;
	}



	int n = (int)pohang.size();
	for(i=0;i<n;i++) {
		potri = pohang[i];
		flg = 0;
		for(j=0;j<3;j++) {
			neib = potri->getNeighbor(j);
			int nn = neib.getNumberOfElement();
			for(k=0;k<nn;k++) {
				netri = neib.getElement(k);
				if(netri->getState() == -1) {
					flg = 1;
					break;
				}
			}
			if(flg) {
				//conect nw_nb to potri
				potri->eraseNeighborElement(j,netri);
				edge = netri->getNodeNumber(ne_nd);
				nw_nb = netri->getNeighbor(edge);
				potri->mergeNeighbor(j,nw_nb);
	
				//conect potri to nw_nb
				int nnn = nw_nb.getNumberOfElement();
				for(l=0;l<nnn;l++) {
					nbtri = nw_nb.getElement(l);
	
					edge = nbtri->getConection(netri);
					nbtri->eraseNeighborElement(edge,netri);
					nbtri->setNewNeighborElement(edge,potri);

				}
				flg = 0;
			}
		}
		//change node to po_nd
		ndedge = potri->getNodeNumber(ne_nd);
		potri->setNode(ndedge,po_nd);
		///potri->setState(0);
	}
	n = (int)nehang.size();
	for(i=0;i<n;i++) {
		netri = nehang[i];
		erase(netri->getId());
	}
	cr_crowd->erase(ne_nd->getId());

	//swapping
	swapLawson(pohang);

	n = (int)pohang.size();
	for(i=0;i<n;i++) {
		potri = pohang[i];
		potri->setState(0);
	}


	return po_nd;
}


//divHanger
///03.12.10
///divide hanger to remain trianles and erase triangles

void Surface::divHanger(Node *po_nd,Node *ne_nd,deque<Triangle*>&pohang,
						deque<Triangle*>&nehang)
{
	int flg;
	Node *pnd;
	Triangle *hgtri;

	int n = ne_nd->getHangerSize(FG_TYPE_TRI);
	for(int i=0;i<n;i++) {
		hgtri = (Triangle*)ne_nd->getHanger(i,FG_TYPE_TRI);
		flg = 0;
		for(int j=0;j<3;j++) {
			pnd = hgtri->getNode(j);
			if(*pnd == *po_nd) {
				flg = 1;
				break;
			}
		}
		if(flg) {
			hgtri->setState(-1);
			nehang.push_back(hgtri);	//erase triangles
		} else {
			hgtri->setState(1);
			pohang.push_back(hgtri);	//remain triangles
		}
	}

	return;
}


//isPossibleMerge
///03.12.10

int Surface::isPossibleMerge(Node *po_nd,Node *ne_nd,deque<Triangle*>&pohang,
							deque<Triangle*>&nehang)
{
	double area,pre_area,pst_area;
	Triangle *ptri;
	Node *nd1,*nd2;
	int edge;
	double mindens,minarea;

	mindens = cr_crowd->getCurrentMinDensity();
	minarea = mindens*mindens*rate_minarea*0.5;

	pre_area = pst_area = 0.;
	int n = (int)pohang.size();
	for(int i=0;i<n;i++) {
		ptri = pohang[i];
		pre_area += ptri->getArea();
		edge = ptri->getNodeNumber(ne_nd);
		nd1 = ptri->getNode((edge+1)%3);
		nd2 = ptri->getNode((edge+2)%3);
		area = getTriangleArea(po_nd,nd1,nd2);

		//check area
		//if new area is too small( area < minarea ), do not merge
		//but it possible to merge, if you want to do
		if(area < minarea) return 0;
		pst_area += area;
	}
	n = (int)nehang.size();
	for(int i=0;i<n;i++) {
		ptri = nehang[i];
		pre_area += ptri->getArea();
	}

	//if pst_area > pre_area(or pst_area < pre_area ), cannot merge
	if(fabs(pre_area - pst_area) > LIMIT) return 0;


	return 1;
}


//swapLawson
///03.12.10

void Surface::swapLawson(deque<Triangle*>istack)
{
	int i,nn;
	Triangle *lt_tri = NULL;	// this side triangle(left triangle after swap)
	Triangle *rt_tri = NULL;	// opposite side triangle(right triangle after swap)
	Triangle *nb_tri = NULL;
	int op_tp_edge=0;		// top node of opposite triangle
	int op_lt_edge=0;		// left node of opposite triangle
	int op_rt_edge=0;		// right node of opposite triangle
	int th_lt_edge=0;		// left node of this side
	int th_rt_edge=0;		// right node of this side
	int nb_edge;
	Neighbor op_ls_nb;
	Neighbor op_rs_nb;
	Neighbor th_ls_nb;
	Neighbor th_rs_nb;
	int type,domain;
	Node *bt_nd = NULL,*tp_nd = NULL,*lt_nd = NULL,*rt_nd = NULL;
	double hori_dist,vert_dist;
	double angl;
	double prev_area1,prev_area2,new_area1,new_area2;
	int flg=0;
	double mindens,minarea;

	mindens = cr_crowd->getCurrentMinDensity();
	minarea = mindens*mindens*rate_minarea*0.5;

	while(istack.size()) {

		//pop triangle from stack
		lt_tri = istack[istack.size()-1];
		istack.pop_back();

		if(lt_tri->getState() > 10)continue;
		for(i=0;i<3;i++) {
			flg = 1;
			if(lt_tri->getNumberOfNeighborElement(i) != 1)continue;

			/////////////////////////////////////////
			rt_tri = lt_tri->getNeighborElement(i,0);

			//////////////if(rt_tri == NULL) getch();

			///////////////
			if(rt_tri->getState() == 0) continue;

			op_tp_edge = rt_tri->getConection(lt_tri);
			op_lt_edge = (op_tp_edge+1)%3;
			op_rt_edge = (op_tp_edge+2)%3;
			
			bt_nd = lt_tri->getNode(i);

			tp_nd = rt_tri->getNode(op_tp_edge);
			lt_nd = rt_tri->getNode(op_lt_edge);
			rt_nd = rt_tri->getNode(op_rt_edge);

			//compare horizontal line with vertical line
			hori_dist = getDistanceSquare(lt_nd,rt_nd);
			vert_dist = getDistanceSquare(tp_nd,bt_nd);

			if(hori_dist < vert_dist) continue;

			//check angle if flat
			angl = lt_tri->getAngle(rt_tri);


			if(angl < FRAT_ANGLE) continue;

			//check area if same
			prev_area1 = getTriangleArea(bt_nd,rt_nd,lt_nd);
			prev_area2 = getTriangleArea(tp_nd,lt_nd,rt_nd);
			new_area1 = getTriangleArea(bt_nd,tp_nd,lt_nd);
			new_area2 = getTriangleArea(bt_nd,rt_nd,tp_nd);

			if(new_area1 < minarea || new_area2 < minarea)continue;
			if( fabs((new_area1+new_area2)-(prev_area1+prev_area2)) > LIMIT_10 )continue;
			
			flg = 0;
			break;
		}
		if(flg) continue;

		//preparation of swapping

		th_lt_edge = lt_tri->getNodeNumber(lt_nd);
		th_rt_edge = lt_tri->getNodeNumber(rt_nd);

		/////////
		//printf("\nedge lt%5d rt%5d\n\n",th_lt_edge,th_rt_edge);
		
		th_ls_nb = lt_tri->getNeighbor(th_lt_edge);
		th_rs_nb = lt_tri->getNeighbor(th_rt_edge);
		op_ls_nb = rt_tri->getNeighbor(op_lt_edge);
		op_rs_nb = rt_tri->getNeighbor(op_rt_edge);


		//swap lt_tri(this side to left)
		type = lt_tri->getType();
		domain = lt_tri->getDomain();
		lt_tri->resetNode(bt_nd,tp_nd,lt_nd);
		lt_tri->resetNeighbor(0,op_rs_nb);
		lt_tri->resetNeighbor(1,th_rs_nb);
		lt_tri->resetNeighborElement(2,rt_tri);
		lt_tri->setType(type);
		lt_tri->setDomain(domain);

		//swap rt_tri(opposit to right)
		type = rt_tri->getType();
		domain = rt_tri->getDomain();
		rt_tri->resetNode(bt_nd,rt_nd,tp_nd);
		rt_tri->resetNeighbor(0,op_ls_nb);
		rt_tri->resetNeighborElement(1,lt_tri);
		rt_tri->resetNeighbor(2,th_ls_nb);
		rt_tri->setType(type);
		rt_tri->setDomain(domain);


		//re-conect neighbors
		nn = op_rs_nb.getNumberOfElement();
		for(i=0;i<nn;i++) {
			nb_tri = op_rs_nb.getElement(i);

			nb_edge = nb_tri->getConection(rt_tri);
			nb_tri->eraseNeighborElement(nb_edge,rt_tri);
			nb_tri->setNewNeighborElement(nb_edge,lt_tri);
		}

		nn = th_ls_nb.getNumberOfElement();
		for(i=0;i<nn;i++) {
			nb_tri = th_ls_nb.getElement(i);

			nb_edge = nb_tri->getConection(lt_tri);
			nb_tri->eraseNeighborElement(nb_edge,lt_tri);
			nb_tri->setNewNeighborElement(nb_edge,rt_tri);

		}

		/*
		if(op_ls_nb.getNumberOfElement() == 1) {
			//push stack next triangle
			istack.push_back(rt_tri);
		}
		if(op_rs_nb.getNumberOfElement() == 1) {
			//push stack next triangle
			istack.push_back(lt_tri);
		}*/

		//push
		lt_tri->setState( lt_tri->getState() + 1 );
		istack.push_back(lt_tri);
	}
		


	return;
}


//qSortDist
///03.11.23
//quick sorting
//used only fineRate

void Surface::qSortDist(int f,int e,deque<double>&dist,deque<Triangle*>&array)
{
	int i,k;
	double pv,key;
	double tmpd;
	Triangle* tmptri;
	int flg;

	if(f >= e) return;

	flg = 1; 
	pv = dist[f];
	for(i=f+1;i<=e;i++) {
		key = dist[i];
		if(key > pv) {
			pv = key;
			flg = 0;
			break;
		} else if(key < pv) {
			flg = 0;
			break;
		}
	}
	if(flg)return;

	k = e;

	for(i=f;i<=k;i++) {
		key = dist[i];
		if(key >= pv) {
			tmpd = dist[i];
			tmptri = array[i];
			dist[i] = dist[k];
			array[i] = array[k];
			dist[k] = tmpd;
			array[k] = tmptri;
			k--;
			i--;
		}
	}

	qSortDist(f,i-1,dist,array);
	qSortDist(i,e,dist,array);

	return;
}




}

//////////////////////////////////////////////////////////////////EOF