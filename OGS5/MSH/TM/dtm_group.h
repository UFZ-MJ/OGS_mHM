/*

			dtm_group.h
			>Group Class
			>>group of elements

			last update : 2003.12.05			by t.manabe

  */

#ifndef DTM_GROUP_H
#define DTM_GROUP_H	0

#include"dtm_fixed.h"

namespace dtm{


	//group of element
	class Group
	{
	private:
		int max_number_of_elements;

	public:
		/* constructor and destructor */

		Group():max_number_of_elements(LIMIT_NUMBER_OF_ELEMENTS){}
		Group(unsigned int max_size):max_number_of_elements(max_size){}
		~Group(){}


		/* original methods */

		virtual int getSize() = 0;
		virtual void clear() = 0;
		virtual int erase(int id) = 0;
		virtual void setMaxSize(unsigned int n){ max_number_of_elements = n; }
		virtual int getMaxSize(){ return max_number_of_elements; }

	};

}

#endif

//////////////////////////////////////////////////////////////////EOF