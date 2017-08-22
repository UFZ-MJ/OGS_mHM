/*

			dtm_error.h
			>Errors Class
			>>throwed when error happend

			last update : 2003.09.24

*/


#ifndef DTM_ERROR_H
#define DTM_ERROR_H	0


#include<iostream>
#include<string>
#include<conio.h>


using namespace std;

namespace dtm{


	class Errors		//throwed each error
	{
	private:
		int ernum;		//error number
		string func;	//function name
		string efect;	//efect of error

	public:
		/* constructor and destructor */

		Errors():ernum(0),func(""),efect("") { }
		Errors(int n,string nf,string ef): ernum(n),func(nf),efect(ef){ }


		/*original methos */

		int getErrorNumber(){ return ernum; }
		string getFunctionName(){ return func; }
		string getEfect(){ return efect; }
	};

//	void errorView(string where,string how);
//	void errorExit();

}


#endif

//////////////////////////////////////////////////////////////////EOF