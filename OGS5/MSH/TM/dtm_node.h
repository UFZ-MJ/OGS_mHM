/*

			dtm_node.h
			>Node class
			>>based Element and Point Class

			last update : 2003.12.03			by t.manabe

*/

#pragma warning (disable : 4786) //for bug of Visual C++ 6.0

#ifndef DTM_NODE_H
#define DTM_NODE_H	1

#include<iostream>
#include<deque>
#include<vector>

#include"dtm_fixed.h"
#include"dtm_element.h"
#include"dtm_point.h"
#include"dtm_figure.h"

using namespace std;

namespace dtm{


	class Node : public Element , public Point	//Node class
	{
	private:
		double cr_density;	//density of node
		int type;			//node type
		/* input node(on surface edge) = 0  on surface(edge) = 1,
		   on surface(not edge) = 2  inside = 3 */

		vector< vector<Figure*> >hanger; //figures which have this node

	public:
		/* constructor and destructor */

		Node();
		Node(int n);
		Node(int n,double x,double y,double z);
		Node(int n,double x,double y,double z,double density);
		Node(int n,Point ppt,double density);
		Node(int n,Point ppt);
		~Node(){}


		/* override of Point methods */

		void setCoordinates(double x,double y,double z){
			Point::setCoordinates(x,y,z);
			renew();
		}
		void setX(double x){ Point::setX(x);	renew();}
		void setY(double y){ Point::setY(y);	renew();}
		void setZ(double z){ Point::setZ(z);	renew();}
		void moveX(double x){	Point::moveX(x);	renew();}
		void moveY(double y){ 	Point::moveY(y);	renew();}
		void moveZ(double z){ 	Point::moveZ(z);	renew();}
		void move(Point mp)	{	Point::move(mp);	renew();}
	
		void reduce(double rate){
			Point::reduce(rate);
			cr_density = cr_density/rate;
			renew();
		}
		void expand(double rate){
			Point::expand(rate);
			cr_density = cr_density*rate;
			renew();
		}

		/* original methods */

		void setDensity(double density)
		{ cr_density = density; /*flg_density = true;*/ }
		double getDensity(){ return cr_density; }
		void setType(int n){ type = n; }
		int getType(){ return type; }
	
		Figure* getHanger(int index,int fgType);
		int getHangerSize(int fgType){ return (int)hanger[fgType].size(); }
		void insertHanger(Figure *hang,int fgType);
		void eraseHanger(int id,int fgType);
		void clearHanger(int fgType = -1);

		void renew(){
			int n = (int)hanger.size();
			for(int i=0;i<n;i++) {
				int nn = (int)hanger[i].size();
				for(int j=0;j<nn;j++) {
					hanger[i][j]->renew();
				}
			}
		}	//renew parameter of hanger objects
		

		/* friend function */

		friend bool operator==(Node &op1,Node &op2) {
			return ( op1.getId() == op2.getId() );
		}
	

		//for debug
		virtual void view();
		virtual void viewHanger(int fgType);
		///
	};


}



#endif

//////////////////////////////////////////////////////////////////EOF