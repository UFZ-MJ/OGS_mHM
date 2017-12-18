/*

			dtm_neighbor.cpp
			>Neighbor Class

			last update : 2003.12.03			by t.manabe

  */

#include "stdafx.h" /* MFC */


#include"dtm_neighbor.h"

namespace dtm{


//Neighbor
///03.12.03

Neighbor::Neighbor():array_size(1), number_of_element(0)
{
	element_array = new Triangle*[array_size];
	element_array[0] = NULL;
}


//Neighbor
///03.12.03

Neighbor::Neighbor(Triangle *elm):array_size(2), number_of_element(0)
{
	if(elm) {
		element_array = new Triangle*[array_size];
		element_array[0] = elm;
		element_array[1] = NULL;
		number_of_element++;
	} else {
		element_array = new Triangle*[array_size];
		element_array[0] = NULL;
		element_array[1] = NULL;
	}

}


//Neighbor
///03.12.03

Neighbor::Neighbor(const Neighbor& op)
{
	clearElement();
	int n = op.getNumberOfElement();
	for(int i=0;i<n;i++) {
		pushElement(op.getElement(i));
	}
}


//~Neighbor
///03.12.03

Neighbor::~Neighbor()
{
	delete [] element_array;
}


//getElement
///03.12.03

Triangle* Neighbor::getElement(int index) const
{
	int n = number_of_element;
	if(index >= n || index < 0 ) return NULL;
	else return element_array[index];
}


//resetElement
///03.12.03

int Neighbor::resetElement(int index,Triangle *elm)
{
	int n = number_of_element;
	if(index < 0 ){
		Errors err(0,
			"int Neighbor::resetElement(int index,Triangle *elm)",
			"illegal index");
		throw err;
	} else if(index >= n) {
		pushElement(elm);
	}
	if(elm) {
		element_array[index] = elm;
	} else {
		eraseElement(index);
	}
	return number_of_element;
}


//pushElement
///03.12.03

int Neighbor::pushElement(Triangle *elm)
{
	int n = number_of_element;
	if(!elm) return number_of_element;
	for(int i=0;i<n;i++) {
		if(element_array[i] == elm) {
			return number_of_element;
		}
	}
	int nn = array_size;
	if(n >= nn) {
		Triangle** tmp;
		tmp = new Triangle*[array_size*2];
		if(!tmp) {
			Errors err(-1,
				"int Neighbor::pushElement(Triangle *elm)",
				"failed making new object");
			throw err;
		}
		for(int i=0;i<n;i++) {
			tmp[i] = element_array[i];
		}
		for(int i=nn;i<nn*2;i++) {
			tmp[i] = NULL;
		}
		delete [] element_array;
		element_array = tmp;
		array_size *= 2;
	}
	element_array[number_of_element] = elm;
	number_of_element++;

	return number_of_element;
}


//eraseElement
///03.12.03

int Neighbor::eraseElement(int index)
{
	int n = number_of_element;
	if(index >= n || index < 0 ){
		/*Errors err(0,
			"int Neighbor::eraseElement(int index)",
			"illegal index");
		throw err;*/
		return -1;
	}
	element_array[index] = element_array[number_of_element-1];
	element_array[number_of_element-1] = NULL;
	number_of_element--;
	return number_of_element;
}


//eraseElement
///03.12.03

int Neighbor::eraseElement(Triangle *elm)
{
	int index = -1;
	int n = number_of_element;
	for(int i=0;i<n;i++) {
		if(element_array[i] == elm) {
			index = i;
			break;
		}
	}
	if(index == -1) return -1;
	element_array[index] = element_array[number_of_element-1];
	element_array[number_of_element-1] = NULL;
	number_of_element--;
	return number_of_element;
}


//clearElement
///03.12.03

void Neighbor::clearElement()
{
	int n = number_of_element;
	for(int i=0;i<n;i++) {
		element_array[i] = NULL;
	}
	number_of_element = 0;
}


//operator=
///03.12.03

Neighbor& Neighbor::operator =(const Neighbor& op)
{
	clearElement();
	int n = op.getNumberOfElement();
	for(int i=0;i<n;i++) {
		pushElement(op.getElement(i));
	}

	return *this;
}


/*
void Neighbor::view()
{
	int n = op.getNumberOfElement();
	for(int i=0;i<n;i++) {
		element_array[i]
}
*/

}

//////////////////////////////////////////////////////////////////EOF