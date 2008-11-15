#include "sdplrlib.h"

double linesearch(problemdata* data, double* R, double* D, double max, double* funcval, int update)
{
  int i;

  double ss=0.0;
  double quartic[5], cubic[4], roots[3], f0, fmax, f1, f2, f3, minval;

  extern double *UVt;
  extern double *ARD;
  extern double *ADD;

  // This code assumes that the violation for R is stored in data->vio,
  // i.e., function() has not been called since function(R) was
  // evaluated. Same for RRt.

  ZEROVEC(ARD, data->m); ARD[0] = 0.0;
  ZEROVEC(ADD, data->m); ADD[0] = 0.0;

  Aoper(data,R,D,UVt,0,1,ARD);
  EASYDSCAL(data->m,2.0,ARD); ARD[0] *= 2.0;

  Aoper(data,D,D,UVt,1,1,ADD);

  // Calculate constant for quartic
  quartic[0] = data->vio[0] -
               EASYDDOT(data->m, data->lambda, data->vio) + 
               0.5*data->sigma*pow(EASYDNRM2(data->m, data->vio), 2.0);
//   quartic[0] = *funcval;

  // Calculate linear for quartic
  quartic[1] = ARD[0] - EASYDDOT(data->m, data->lambda, ARD) +
               data->sigma*EASYDDOT(data->m, data->vio, ARD);

  // Calculate quadratic for quartic
  quartic[2] = ADD[0] - EASYDDOT(data->m, data->lambda, ADD) +
               data->sigma*(EASYDDOT(data->m, data->vio, ADD) +
               0.5*pow(EASYDNRM2(data->m, ARD), 2.0));

  // Calculate cubic for quartic
  quartic[3] = data->sigma*EASYDDOT(data->m, ARD, ADD);

  // Calculate quartic for quartic
  quartic[4] = 0.5*data->sigma*pow(EASYDNRM2(data->m, ADD), 2.0);

  // Calculate cubic
  cubic[0] = 1.0*quartic[1]; if(cubic[0] > 0.0) {
    printf("Problem!  %f should be less then 0.\n", cubic[0]); return 0.0; }
  cubic[1] = 2.0*quartic[2];
  cubic[2] = 3.0*quartic[3];
  cubic[3] = 4.0*quartic[4];

  if(cubic[3] == 0.0) {
    printf("Surprise! Got a quadratic function!\n");
    exit(0);
  }

  roots[0] = roots[1] = roots[2] = 1.0e10;
  gsl_poly_solve_cubic(cubic[2]/cubic[3], cubic[1]/cubic[3], cubic[0]/cubic[3], roots, roots+1, roots+2);
  //printf("%f    %f    %f\n", roots[0], roots[1], roots[2]);

  f0 = quartic[0];

  fmax = gsl_poly_eval(quartic, 5, max);

  if(roots[0] == 1.0e10 || roots[0] <= 0.0 || roots[0] > max) f1 = 1.0e20;
  else f1 = gsl_poly_eval(quartic, 5, roots[0]);

  if(roots[1] == 1.0e10 || roots[1] <= 0.0 || roots[1] > max) f2 = 1.0e20;
  else f2 = gsl_poly_eval(quartic, 5, roots[1]);

  if(roots[2] == 1.0e10 || roots[2] <= 0.0 || roots[2] > max) f3 = 1.0e20;
  else f3 = gsl_poly_eval(quartic, 5, roots[2]);

  minval = 1.0e20;
  minval = mymin(f0,minval);
  minval = mymin(fmax,minval);
  minval = mymin(f1,minval);
  minval = mymin(f2,minval);
  minval = mymin(f3,minval);

  if(f0 == minval) ss = 0.0;
  if(fmax == minval) ss = max;
  if(f1 == minval) ss = roots[0];
  if(f2 == minval) ss = roots[1];
  if(f3 == minval) ss = roots[2];

  // Setup next values

  *funcval = minval;
  if(update) {
    for(i = 1; i <= data->m; i++)
      data->vio[i] += ss*(ARD[i] + ss*ADD[i]);
    data->vio[0] += ss*(ARD[0] + ss*ADD[0]);
  }

  return ss;
}

