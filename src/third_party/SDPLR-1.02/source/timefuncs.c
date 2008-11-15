#include "sdplrlib.h"


#ifdef __PC

#include <time.h>

double current_time(clock_t timeorig)
{
  return (double)(clock() - timeorig)/CLOCKS_PER_SEC;
}


#else


#include <sys/times.h>
#include <unistd.h>

double current_time(double timeorig)
{
  int     clocks_per_second = sysconf( _SC_CLK_TCK );
  struct  tms timearr;

  times(&timearr);

  return (timearr.tms_utime - timeorig)/clocks_per_second;
}


#endif

