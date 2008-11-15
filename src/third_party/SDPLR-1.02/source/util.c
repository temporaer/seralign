#include "sdplrlib.h"

#define mymax(A, B) ((A) > (B) ? (A) : (B))
#define mymin(A, B) ((A) < (B) ? (A) : (B))

int    daxpy_();
int    dcopy_();
double ddot_();
double dnrm2_();
int    dscal_();

int copyscaledvectovec(double* dy, double da, double* dx, int n)
{
#ifndef __OPT_MYBLAS
  EASYDCOPY(n, dx, dy);
  EASYDSCAL(n, da, dy);
  return 0;
#elif __OPT_UNROLLEDLOOP
  int i=0, blocklimit;
  blocklimit = (n/8)*8; 
  while(i < blocklimit) {
    dy[i+1] = da*dx[i+1];
    dy[i+2] = da*dx[i+2];
    dy[i+3] = da*dx[i+3];
    dy[i+4] = da*dx[i+4];
    dy[i+5] = da*dx[i+5];
    dy[i+6] = da*dx[i+6];
    dy[i+7] = da*dx[i+7];
    dy[i+8] = da*dx[i+8];
    i += 8; 
  }
  if(i < n) { 
    switch(n - i) { 
      case 7 : dy[i+1] = da*dx[i+1]; i++; 
      case 6 : dy[i+1] = da*dx[i+1]; i++; 
      case 5 : dy[i+1] = da*dx[i+1]; i++; 
      case 4 : dy[i+1] = da*dx[i+1]; i++; 
      case 3 : dy[i+1] = da*dx[i+1]; i++; 
      case 2 : dy[i+1] = da*dx[i+1]; i++; 
      case 1 : dy[i+1] = da*dx[i+1]; 
    }
  }
  return 0;
#else
  int i;
  for(i = 1; i <= n; i++)
    dy[i] = da*dx[i];
  return 0;
#endif
}

int move_in_dir(double* dy, double* dx, double da, double* dir, int n)
{
#ifndef __OPT_MYBLAS
  if(dy == dx) EASYDAXPY(n, da, dir, dy);
  else if(dy == dir) {
    EASYDSCAL(n, da, dy);
    EASYDAXPY(n, 1.0, dx, dy);
  }
  else {
    EASYDCOPY(n, dx, dy);
    EASYDAXPY(n, da, dir, dy);
  }
  return 1;
#elif __OPT_UNROLLEDLOOP
  int i=0, blocklimit;
  blocklimit = (n/8)*8; 
  while(i < blocklimit) {
    dy[i+1] = dx[i+1] + da*dir[i+1];
    dy[i+2] = dx[i+2] + da*dir[i+2];
    dy[i+3] = dx[i+3] + da*dir[i+3];
    dy[i+4] = dx[i+4] + da*dir[i+4];
    dy[i+5] = dx[i+5] + da*dir[i+5];
    dy[i+6] = dx[i+6] + da*dir[i+6];
    dy[i+7] = dx[i+7] + da*dir[i+7];
    dy[i+8] = dx[i+8] + da*dir[i+8];
    i += 8; 
  }
  if(i < n) { 
    switch(n - i) { 
      case 7 : dy[i+1] = dx[i+1] + da*dir[i+1]; i++; 
      case 6 : dy[i+1] = dx[i+1] + da*dir[i+1]; i++; 
      case 5 : dy[i+1] = dx[i+1] + da*dir[i+1]; i++; 
      case 4 : dy[i+1] = dx[i+1] + da*dir[i+1]; i++; 
      case 3 : dy[i+1] = dx[i+1] + da*dir[i+1]; i++; 
      case 2 : dy[i+1] = dx[i+1] + da*dir[i+1]; i++; 
      case 1 : dy[i+1] = dx[i+1] + da*dir[i+1]; 
    }
  }
  return 0;
#else
  int i;
  for(i = 1; i <= n; i++)
    dy[i] = dx[i] + da*dir[i];
  return 0;
#endif
}

int mydaxpy(int n, double da, double* dx, int incx, double* dy, int incy)
{
#ifndef __OPT_MYBLAS
  return daxpy_(&n,&da,dx,&incx,dy,&incy);
#elif __OPT_UNROLLEDLOOP
  int i=0, blocklimit, ix=0, iy=0;
  if(incx == 1 && incy == 1) {
    blocklimit = (n/8)*8; 
    while(i < blocklimit) {
      dy[i] += da*dx[i];
      dy[i+1] += da*dx[i+1];
      dy[i+2] += da*dx[i+2];
      dy[i+3] += da*dx[i+3];
      dy[i+4] += da*dx[i+4];
      dy[i+5] += da*dx[i+5];
      dy[i+6] += da*dx[i+6];
      dy[i+7] += da*dx[i+7];
      i += 8; 
    }
    if(i < n) { 
      switch(n - i) { 
        case 7 : dy[i] += da*dx[i]; i++; 
        case 6 : dy[i] += da*dx[i]; i++; 
        case 5 : dy[i] += da*dx[i]; i++; 
        case 4 : dy[i] += da*dx[i]; i++; 
        case 3 : dy[i] += da*dx[i]; i++; 
        case 2 : dy[i] += da*dx[i]; i++; 
        case 1 : dy[i] += da*dx[i]; 
      }
    }
  }
  else {
    for(i = n; i--; ) {
      dy[iy] += da*dx[ix];
      ix += incx;
      iy += incy;
    }
  }
  return 0;
#else
  int i, ix=0, iy=0;
  if(incx == 1 && incy == 1) {
    for(i = n; i--; )
      dy[i] += da*dx[i];
  }
  else {
    for(i = n; i--; ) {
      dy[iy] += da*dx[ix];
      ix += incx;
      iy += incy;
    }
  }
  return 0;
#endif
}

int mydcopy(int n, double* dx, int incx, double* dy, int incy)
{
#ifndef __OPT_MYBLAS
  return dcopy_(&n,dx,&incx,dy,&incy);
#elif __OPT_UNROLLEDLOOP
  int i=0, blocklimit, ix=0, iy=0;
  if(incx == 1 && incy == 1) {
    blocklimit = (n/8)*8; 
    while(i < blocklimit) {
      dy[i] = dx[i];
      dy[i+1] = dx[i+1];
      dy[i+2] = dx[i+2];
      dy[i+3] = dx[i+3];
      dy[i+4] = dx[i+4];
      dy[i+5] = dx[i+5];
      dy[i+6] = dx[i+6];
      dy[i+7] = dx[i+7];
      i += 8; 
    }
    if(i < n) { 
      switch(n - i) { 
        case 7 : dy[i] = dx[i]; i++; 
        case 6 : dy[i] = dx[i]; i++; 
        case 5 : dy[i] = dx[i]; i++; 
        case 4 : dy[i] = dx[i]; i++; 
        case 3 : dy[i] = dx[i]; i++; 
        case 2 : dy[i] = dx[i]; i++; 
        case 1 : dy[i] = dx[i]; 
      }
    }
  }
  else {
    for(i = n; i--; ) {
      dy[iy] = dx[ix];
      ix += incx;
      iy += incy;
    }
  }
  return 0;
#else
  int i, ix=0, iy=0;
  if(incx == 1 && incy == 1) {
    for(i = n; i--; )
      dy[i] = dx[i];
  }
  else {
    for(i = n; i--; ) {
      dy[iy] = dx[ix];
      ix += incx;
      iy += incy;
    }
  }
  return 0;
#endif
}

double myddot(int n, double* dx, int incx, double* dy, int incy)
{
/*
#ifndef __OPT_MYBLAS
*/
  return ddot_(&n, dx, &incx, dy, &incy);
/*
#elif __OPT_UNROLLEDLOOP
  int i=0, blocklimit, ix=0, iy=0;
  double sum=0.0;
  if(incx == 1 && incy == 1) {
    blocklimit = (n/8)*8; 
    while(i < blocklimit) {
      sum += dy[i]*dx[i] + dy[i+1]*dx[i+1] + dy[i+2]*dx[i+2] +
             dy[i+3]*dx[i+3] + dy[i+4]*dx[i+4] + dy[i+5]*dx[i+5] +
             dy[i+6]*dx[i+6] + dy[i+7]*dx[i+7];
      i += 8; 
    }
    if(i < n) { 
      switch(n - i) { 
        case 7 : sum += dy[i]*dx[i]; i++; 
        case 6 : sum += dy[i]*dx[i]; i++; 
        case 5 : sum += dy[i]*dx[i]; i++; 
        case 4 : sum += dy[i]*dx[i]; i++; 
        case 3 : sum += dy[i]*dx[i]; i++; 
        case 2 : sum += dy[i]*dx[i]; i++; 
        case 1 : sum += dy[i]*dx[i]; 
      }
    }
  }
  else {
    for(i = n; i--; ) {
      sum += dx[ix]*dy[iy];
      ix += incx;
      iy += incy;
    }
  }
  return sum;
#else
  int i, ix=0, iy=0;
  double sum=0.0;
  if(incx == 1 && incy == 1) {
    for(i = n; i--; )
      sum += dx[i]*dy[i];
  }
  else {
    for(i = n; i--; ) {
      sum += dx[ix]*dy[iy];
      ix += incx;
      iy += incy;
    }
  }
  return sum;
#endif
*/
}

double mydnrm2(int n, double* dx, int incx)
{
  return dnrm2_(&n, dx, &incx);
}

int mydscal(int n, double da, double* dx, int incx)
{
#ifndef __OPT_MYBLAS
  return dscal_(&n, &da, dx, &incx);
#elif __OPT_UNROLLEDLOOP
  int i=0, blocklimit, ix=0;
  if(incx == 1) {
    blocklimit = (n/8)*8; 
    while(i < blocklimit) {
      dx[i] *= da;
      dx[i+1] *= da;
      dx[i+2] *= da;
      dx[i+3] *= da;
      dx[i+4] *= da;
      dx[i+5] *= da;
      dx[i+6] *= da;
      dx[i+7] *= da;
      i += 8; 
    }
    if(i < n) { 
      switch(n - i) { 
        case 7 : dx[i] *= da; i++; 
        case 6 : dx[i] *= da; i++; 
        case 5 : dx[i] *= da; i++; 
        case 4 : dx[i] *= da; i++; 
        case 3 : dx[i] *= da; i++; 
        case 2 : dx[i] *= da; i++; 
        case 1 : dx[i] *= da; 
      }
    }
  }
  else {
    for(i = n; i--; ) {
      dx[ix] *= da;
      ix += incx;
    }
  }
  return 0;
#else
  int i, ix=0;
  if(incx == 1) {
    for(i = n; i--; )
      dx[i] *= da;
  }
  else {
    for(i = n; i--; ) {
      dx[ix] *= da;
      ix += incx;
    }
  }
  return 0;
#endif
}


int createlowrankmat(problemdata* data, lowrankmat** passedR, int ncol, int nrow)
{
  lowrankmat *R;

  MYCALLOC(R, lowrankmat, 1);
  
  R->ncol = ncol;
  R->nrow = nrow;
  
  MYCALLOC(R->d, double, ncol + 1);
  MYCALLOC(R->ent, double, nrow*ncol + 1);

  *passedR = R;

  return 1;
}

int destroylowrankmat(lowrankmat* R)
{
  MYFREE(R->d);
  MYFREE(R->ent);
  MYFREE(R);

  return 1;
}

int createsparsesymmmat(sparsesymmmat** passedS, int nnz)
{
  sparsesymmmat *S;

  MYCALLOC(S, sparsesymmmat, 1);
  MYCALLOC(S->row, int, nnz+1);
  MYCALLOC(S->col, int, nnz+1);
  S->nnz = nnz;
  MYCALLOC(S->ent, double, nnz+1);
//   MYCALLOC(S->inagg, int, nnz+1);

  *passedS = S;

  return 1;

}

int destroysparsesymmmat(sparsesymmmat* S)
{
  MYFREE(S->row);
  MYFREE(S->col);
  MYFREE(S->ent);
//   MYFREE(S->inagg);
  MYFREE(S);

  return 1;
}

int creatediagmat(diagmat** passedD, int nnz)
{
  diagmat *D;

  MYCALLOC(D, diagmat, 1);
  MYCALLOC(D->ind, int, nnz+1);
  D->nnz = nnz;
  MYCALLOC(D->ent, double, nnz+1);

  *passedD = D;

  return 1;

}

int destroydiagmat(diagmat* D)
{
  MYFREE(D->ind);
  MYFREE(D->ent);
  MYFREE(D);

  return 1;
}


int createdatamat(problemdata* data, datamat** passedA, char type, int ncol_or_nnz, int dim, char* label)
{
  datamat *A;

  MYCALLOC(A, datamat, 1);
  A->type = type;
  MYCALLOC(A->label, char, 30);
  strcpy(A->label, label);

  if(type == 'l') createlowrankmat(data, &(A->lr), ncol_or_nnz, dim);
  
  if(type == 's') {
    createsparsesymmmat(&(A->sp), ncol_or_nnz);
    MYCALLOC(A->sp->XS_in, int, ncol_or_nnz + 1);
  }

  if(type == 'd') {
    creatediagmat(&(A->diag), ncol_or_nnz);
    MYCALLOC(A->diag->XS_in, int, ncol_or_nnz + 1);
  }

  // if type = 'u' then do nothing

  *passedA = A;

  return 1;
}

int destroydatamat(datamat* A)
{
  if(A->type == 'l') destroylowrankmat(A->lr);
  
  if(A->type == 's') {
    destroysparsesymmmat(A->sp);
    MYFREE(A->sp->XS_in);
  }

  if(A->type == 'd') {
    destroydiagmat(A->diag);
    MYFREE(A->diag->XS_in);
  }

  MYFREE(A->label);
  MYFREE(A);

  return 1;
}

