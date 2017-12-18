#ifndef _BENCHTIMER_H
#define _BENCHTIMER_H

#ifndef _WIN32
#include <sys/time.h>
#endif


class BenchTimer {
public:
  
  BenchTimer();
  ~BenchTimer();
  
  void start();
  void stop();
  
  bool running();
  
  double time_s();
  double time_ms();
  
  double precision_s();
  double precision_ms();
  
  
protected:

#ifndef _WIN32
  struct timeval start_val;
  struct timeval end_val;
#endif

  double prec_s;
  bool is_running;
  
};

#endif
