/*

			element.h
			>Element class
			>>base class of all objects(node,figure,triangle,tetra,cube)

			last update : 2003.11.23		by t.manabe

*/


#ifndef DTM_ELEMENT_H
#define DTM_ELEMENT_H	0

namespace dtm{
	
	class Element //base class of all objects
	{
	private:
		int element_id;	//ID of element

	public:
		/* constructor and destructor */

		Element():element_id(0){}
		Element(int id):element_id(id){}
		virtual ~Element(){}


		/* original methods */

		virtual void setId(int id){element_id = id;}
		virtual int getId() const {return element_id;}


		/* friend function */

		friend bool operator==(Element &op1,Element &op2){
			return (op1.element_id == op2.element_id);
		}
	};


}






#endif


//////////////////////////////////////////////////////////////////EOF