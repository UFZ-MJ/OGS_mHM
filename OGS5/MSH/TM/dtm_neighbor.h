/*

			dtm_neighbor.h
			>Neighbor class for Triangle

			last update : 2003.12.03			by t.manabe

  */


#ifndef DTM_NEIGHBOR_H
#define DTM_NEIGHBOR_H	2

#include"dtm_error.h"
//#include"dtm_figure.h"
//#include"dtm_triangle.h"

namespace dtm{

	class Triangle;


	class Neighbor 
	{
		private:
			unsigned int array_size;
			unsigned int number_of_element;
			Triangle **element_array;

		public:
			/* constructor and destructor */

			Neighbor();
			Neighbor(Triangle *elm);
			Neighbor(const Neighbor& op);
			~Neighbor();


			/* original methods */

			int getNumberOfElement() const { return number_of_element; }
			int getArraySize() const { return array_size; }

			Triangle* getElement(int index) const;
			int resetElement(int index,Triangle *elm);
			int pushElement(Triangle *elm);
			int eraseElement(int index);
			int eraseElement(Triangle *elm);
			void clearElement();

			Neighbor& operator=(const Neighbor& op);

			//for debug
			//void view();
	};




}

#endif

//////////////////////////////////////////////////////////////////EOF