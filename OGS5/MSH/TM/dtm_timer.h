/*







*/


#ifndef DTM_TIMER_H
#define DTM_TIMER_H	1

#include<time.h>
#include<iostream>
#include<cstring>

using namespace std;

namespace dtm{

	class Timer
	{
	private:

		clock_t start;
		clock_t stop;
		clock_t total;
		int count;

		char name[128];

	public:
		Timer(): total(0) , count(0){ strcpy(name,""); }
		Timer(char *nm): total(0) , count(0){ strcpy(name,nm); }
		~Timer(){ }

		virtual clock_t Start(){
			total = 0;
			start = clock();
			return start;
		}
		virtual clock_t Stop(){
			stop = clock();
			total += stop - start;
			count++;
			return stop;
		}
		virtual void reStart(){
			start = clock();
		}
		virtual void clear(){
			total = 0;
			count = 0;
		}
		virtual clock_t getTime(){ return total; }
		virtual clock_t getStart(){ return start; }
		virtual clock_t getStop(){ return stop; }
		virtual void view(){
			std::cout << name << " TIME : " << total << " (" 
				<< total/CLOCKS_PER_SEC << "[s])\n";
		}
		virtual clock_t getAverage(){
			if(!count) return 0;
			return total/count;
		}
		virtual void viewAverage(){
			if(!count) return;
			std::cout << name << "TIME average : " << total/count
				<< " count : " << count << "\n";
		}



	};


}





#endif


//////////////////////////////////////////////////////////////////EOF