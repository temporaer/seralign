#include "sdplrlib.h"

#ifdef __ARPACK

int my_arpack(int (*matvec)(double*, double*), int n, int nev, int ncv, char* which, int maxitr, int printlevel, double* evals, double* evecs, int* nconv, int* nummatvec)
{
  // Variables taken from ARPACK examples
  double *v, *workl, *workd, *d, *resid;
  int *select, iparam[11], ipntr[11];
  char bmat[2], all[4];
  int ldv, ido, lworkl, info, ierr;
  int mode, ishfts, rvec;
  double tol, sigma;
  double zero;

  int i;
  double minevalest;

  // Extra variables
  int tempint;

  // Initialize counter for number of matvec operations
  *nummatvec = 0;

  // Setup different parameters
  ldv = n;
  zero = 0.0;
  strcpy(bmat, "I\0");
  lworkl = ncv*(ncv+8);
  //tol = zero;
  tol = (double)1.0/n;
  info = 0;
  ido = 0;
  ishfts = 1;
  mode = 1;
  iparam[1-1] = ishfts;
  iparam[3-1] = maxitr;
  iparam[7-1] = mode;
  rvec = 1;
  strcpy(all, "All\0");

  // Allocate memory
  v = (double*)calloc(ldv*ncv, sizeof(double));
  workl = (double*)calloc(ncv*(ncv+8), sizeof(double));
  workd = (double*)calloc(3*n, sizeof(double));
  d = (double*)calloc(ncv*2, sizeof(double));
  resid = (double*)calloc(n, sizeof(double));
  select = (int*)calloc(ncv, sizeof(int));

  // Do ARPACK loop

  do {

    dsaupd_(&ido, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &info);

    if(ido == -1 || ido == 1) {

      // This line assumes matvec numbers from 0 (not 1)
      matvec(workd + ipntr[2-1] - 1, workd + ipntr[1-1] - 1);
      *nummatvec = *nummatvec + 1;

    }
    else if(info < 0) {

      if(printlevel > 0) printf("ARPACK error (info = %d).\n", info);
      goto END_ROUTINE;

    }

  } while(ido != 99);

  // Post process

  dseupd_(&rvec, all, select, d, v, &ldv, &sigma, bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &ierr);

  if(ierr != 0) {
    if(printlevel > 0) printf("ARPACK error (ierr = %d).\n", ierr);
    info = ierr;
    goto END_ROUTINE;
  }

  *nconv = iparam[5-1];

  if(printlevel > 0 && *nconv < nev) {
    printf("Warning: %d out of %d evals converged.", *nconv, nev);
    if(nev == 1) printf(" Best estimate returned.\n");
    else printf("\n");
  }

  // Save eigen pairs in user's memory (even unreliable ones)
  mydcopy(nev, d, 1, evals, 1);
  tempint = nev*n;
  mydcopy(tempint, v, 1, evecs, 1);

  if(nev == 1 && *nconv == 0) {
    minevalest = 1.0e10;
    for(i = 1; i <= ncv; i++)
      minevalest = mymin(minevalest, workl[ipntr[6-1]-2+i]);
    evals[0] = minevalest;
  }


END_ROUTINE:

  free(select);
  free(resid);
  free(d);
  free(workd);
  free(workl);
  free(v);


  if(*nconv < nev) return -1;
  return info;
}

#endif

int simple_getWWt(double* out, double* in)
{
  extern problemdata* global_data;

  // This routine is called by ARPACK, which assumes
  // out and in are numbered from 0, while getWWt assumes from 1
  getWWt(global_data, global_copyR, out - 1, in - 1);

  return 0;
}

int simple_Stimesvec_block(double* out, double* in)
{
  extern int          global_blk;
  extern problemdata *global_data;

  // This routine is called by ARPACK, which assumes
  // out and in are numbered from 0

  Stimesmat(global_data, global_data->S + global_data->XS_blkptr[global_blk] - 1, global_data->y, in - 1, out - 1, global_data->blksz[global_blk], 1, global_blk, 0);

  return 0;
}

int simple_getMv(double* out, double* in)
{
  extern problemdata *global_data;
  extern double *global_tempvec;

  // This routine is called by ARPACK, which assumes
  // out and in are numbered from 0

  getMv(global_data, global_copyR, in-1, global_data->S, out-1, global_tempvec);

  return 0;
}


int Smineval(problemdata* d, double* eval)
{
  int problem, ct, i, k;
  double *blkevals;

  *eval = 1.0e10;

  ct = 0;
  for(k = 1; k <= d->numblk; k++)
    if(d->blktype[k] == SDPBLK) ct++;
    else if(d->blktype[k] == DIAGBLK) ct += d->blksz[k];

  MYCALLOC(blkevals, double, ct+1);

  problem = Sblockmineval(d, blkevals);

  for(i = 1; i <= ct; i++)
    *eval = mymin(*eval, blkevals[i]);

  MYFREE(blkevals);

  return problem;
}

int Sblockmineval(problemdata* d, double* blkevals)
{
  int problem=0;
  int ct, i, k, ret, n, nev, ncv, maxitr, printlevel;
  int nummatvec, nconv;
  char which[3];
  double eval, *evec;
  int maxblksz, info, dimwork;
  char uplo, jobz;
  double *work, *Scopy, *evals;
  extern int          global_blk;
  extern problemdata *global_data;

  // Note this function overwrites S!

  global_data = d;
  nev = 1;
  strcpy(which, "SA\0");
  maxitr = 500;
  printlevel = 0;


  // If necessary, allocate workspace for dsyev_
  maxblksz = 0;
  for(k = 1; k <= d->numblk; k++)
    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE)
      maxblksz = mymax(maxblksz, d->blksz[k]);
  if(maxblksz > 0) {
    dimwork = mymax(1, 3*maxblksz - 1);
    MYCALLOC(Scopy, double, maxblksz*maxblksz + 1);
    MYCALLOC(work, double, dimwork + 1);
    MYCALLOC(evals, double, maxblksz + 1);
  }
  else { work = Scopy = evals = NULL; }

  ct = 0;

  for(k = 1; k <= d->numblk; k++) {

    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == SPARSE) {
      
#ifdef __ARPACK

      global_blk = k;
      n = d->blksz[k];
      ncv = mymin(2*mymax(2*nev,20),n); // This is Matlab's rule (times two).
      MYCALLOC(evec, double, n+1);

      // When finished, evals[1] gives first eigenvalue
      // and evecs[1] is start of eigenvectors
      ret = my_arpack(simple_Stimesvec_block, n, nev, ncv, which, maxitr, printlevel, &eval, evec+1, &nconv, &nummatvec);
      MYFREE(evec);
      if(ret == -1) problem = -1;

#else

      eval = -1.0e10;

#endif

      blkevals[++ct] = eval;

    }
    else if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE) {

        jobz = 'n'; uplo = 'l';
        EASYDCOPY(d->blksz[k]*d->blksz[k], d->S + d->XS_blkptr[k] - 1, Scopy);
        dsyev_(&jobz, &uplo, &(d->blksz[k]), Scopy + 1, &(d->blksz[k]), evals + 1, work + 1, &dimwork, &info);
        if(info != 0) { printf("Eigenvalue computation failed.\n"); exit(0); }
        blkevals[++ct] = evals[1];

    }
    else if(d->blktype[k] == DIAGBLK) {

      for(i = 1; i <= d->blksz[k]; i++)
        blkevals[++ct] = d->S[d->XS_blkptr[k] - 1 + i];

    }

  }

  // If necessary, deallocate workspace used for dsyev_
  if(maxblksz > 0) {
    MYFREE(Scopy);
    MYFREE(work);
    MYFREE(evals);
  }

  return problem;
}

// double dualbound(problemdata* data)
// {
//   int problem;
//   double eval;

//   EASYDCOPY(data->m, data->lambda, data->S[1]->y);
//   problem = Smineval(data, &eval);
//  
//   return EASYDDOT(data->m, data->b, data->S[1]->y) - data->evaladj*mymin(eval, 0.0);

// }


