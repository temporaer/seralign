#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifdef __PC
  #include <time.h>
#else
  #include <sys/times.h>
  #include <unistd.h>
#endif

#ifdef __USER
  #include "user.h"
#endif

