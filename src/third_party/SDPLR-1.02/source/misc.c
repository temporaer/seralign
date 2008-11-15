#include "sdplrlib.h"
#include "version.h"

int print_notes(void)
{
  printf("\n");
  printf("Note: Readdata_sdplr() assumes only one block!\n");
  printf("Note: Dual bounds currently only partially enabled!\n");
  printf("Note: Currently, ARPACK tol = 1/n. Appropriate?\n");
  printf("Note: ARPACK does not print error when it fails to converge.\n");
  printf("Note: No PC'ing for diagonal blocks.\n");
  printf("Note: Scaling doesn't work b/c of low-rank data matrices.\n");
  printf("\n");

  return 0;
}

int printparams(problemdata* data)
{
  printf("rho_f        = %.1e\n", data->rho_f);
  printf("rho_c        = %.1e\n", data->rho_c);
  printf("method       = %d\n", data->method);
  printf("sigmafac     = %.1f\n", data->sigmafac);
  printf("rankreduce   = %d\n", data->rankreduce);
  printf("timelim      = %d\n", data->timelim);
  printf("dthresh_dim  = %d\n", data->dthresh_dim);
  printf("dthresh_dens = %.2f\n", data->dthresh_dens);

  printf("numbfgsvecs  = %d\n", data->numbfgsvecs);
  printf("mixedtol     = %.1e\n", data->mixedtol);
  printf("doAR         = %d\n", data->doAR);
  printf("rankredtol   = %.16e\n", data->rankredtol);

  printf("pc           = %d\n", data->pc);
  printf("gdens        = %.1f\n", data->gdens);
  printf("reorder      = %d\n", data->reorder);
  printf("gaptol       = %.1e\n", data->gaptol);
  printf("checkbd      = %d\n", data->checkbd);
  printf("typebd       = %d\n", data->typebd);

  return 0;
}


int print_dimacs_errors(problemdata* data, double* R)
{
  int problem;
  int i, j, k, one=1;
  double CdotRRt, bty;
  double Cnorm, bnorm;
  double tempval;
  double *X, *S;
  int *colptr, *rowind;
  double d1,d2,d3,d4,d5,d6;

  d2 = d3 = 0.0;

//   printf("DIMACS error measures: ");

  // This routine does not work if C is a low-rank
  // data matrix. (Have I fixed this?)

  Aoper(data,R,R,data->X,1,1,data->vio);
  EASYDAXPY(data->m, -1.0, data->b, data->vio);
  CdotRRt = data->vio[0];

  bty = EASYDDOT(data->m, data->b, data->lambda);
  Cnorm = C_normdatamat(data);
  bnorm = fabs(data->b[ idamax_(&(data->m), data->b + 1, &one) ]);

  d1 = EASYDNRM2(data->m, data->vio)/(1.0 + bnorm);

#ifdef __ARPACK
  // Begin: Make sure S is based on y only
  for(k = 1; k <= data->m; k++)
    data->y[k] = -data->lambda[k];
  AToper(data, data->y, data->S, 1);
  // End: Make sure S is based on y only
  problem = Smineval(data,&tempval);
  d4 = mymax(0.0, -tempval/(1.0 + Cnorm));
#else
  problem = Smineval(data,&tempval);
  if(tempval != -1.0e10) {
    d4 = mymax(0.0, -tempval/(1.0 + Cnorm));
  }
  else { d4 = -1.0e10; }
#endif

  d5 = (CdotRRt - bty)/(1.0 + fabs(CdotRRt) + fabs(bty));

  tempval = 2.0*EASYDDOT(data->XS_blkptr[data->numblk+1]-1, data->X, data->S);
  for(k = 1; k <= data->numblk; k++) {
    X = data->X + data->XS_blkptr[k] - 1;
    S = data->S + data->XS_blkptr[k] - 1;
    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == SPARSE) {
      colptr = data->XS_colptr[k];
      rowind = data->XS_rowind[k];
      for(j = 1; j <= data->blksz[k]; j++)
        for(i = colptr[j]; i <= colptr[j+1]-1; i++) // This is a bit exorbitant.
          if(rowind[i] == j)
            tempval -= X[i]*S[i];
    }
    else if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == DENSE) {
      for(j = 1; j <= data->blksz[k]; j++)
        tempval -= X[SMATIND(j,j,data->blksz[k])]*S[SMATIND(j,j,data->blksz[k])];
    }
    else if(data->blktype[k] == DIAGBLK) {
      for(j = 1; j <= data->blksz[k]; j++)
        tempval -= X[j]*S[j];
    }
  }

  d6 = tempval/(1.0 + fabs(CdotRRt) + fabs(bty));

  if(d4 != -1.0e10) printf("DIMACS error measures: %.2e %.2e %.2e %.2e %.2e %.2e\n", d1, d2, d3, d4, d5, d6);
  else printf("DIMACS error measures: %.2e %.2e %.2e (n/a) %.2e %.2e\n", d1, d2, d3, d5, d6);
  
  if(problem == -1)
    printf("Warning (ARPACK): Eigenvalue calculation failed to converge. Best estimate returned.\n");

  printf("\n");

  return 0;
}


// int scale(problemdata* data)
// {
//   int i, j, k;
//   double normsq, norm;
//   datamat *A;

//   for(i = 1; i <= data->m; i++) {
//     norm = new_normdatamat(data, i);
//     data->b[i] /= norm;
//     for(k = 1; k <= data->numblk; k++) {
//       if(data->blktype[k] == 's') {
//         if(data->A[i][k]->type == 's') {
//           for(j = 1; j <= data->A[i][k]->sp->nnz; j++) {
//             data->A[i][k]->sp->ent[j] /= norm;
//           }
//         }
//         else if(data->A[i][k]->type == 'l') 
//           EASYDSCAL(data->A[i][k]->lr->ncol, 1.0/norm, data->A[i][k]->lr->d);
//       }
//       else if(data->blktype[k] == 'd')
//         for(j = 1; j <= data->A[i][k]->diag->nnz; j++) 
//           data->A[i][k]->diag->ent[j] /= norm;
//     }
//   }

//   return 0;

// }

// int printstats(problemdata* d)
// {
//   //int k;
//   int i;
//   int ct1,ct2;
//   //datamat *A;

//   //printf("Detailed Problem Statistics:\n\n");

//   printf("m = %d\n", d->m);

//   //printf("\n");

//   printf("blk = %d  ", d->numblk);
//   ct1 = 0; ct2 = 0;
//   for(i = 1; i <= d->numblk; i++) {
//     if(d->blktype[i] == 's') ct1++;
//     if(d->blktype[i] == 'd') ct2++;
//   }
//   //printf("(sdp:  %d     lp:  %d)\n", ct1, ct2);
//   printf("\n");

//   //printf("\n");

//   for(i = 1; i <= d->numblk; i++) {
//     printf("( %c   %d   %d   ", d->blktype[i], d->blksz[i], d->rank[i]);
//     if(d->blktype[i] == 's' && d->XS_blksto[i] == 's') printf("%3.2f )\n", (double)200.0*d->S[i]->sp->nnz/(d->blksz[i]*(d->blksz[i]+1)));
//     else printf("%3.2f )\n", 100.0);
//   }

//   //printf("\n");

//   printf("nr = %d\n", d->nr);

//   printf("\n");


//   /*
//   for(i = 0; i <= d->m; i++) {
//     for(k = 1; k <= d->numblk; k++) {
//       if(i == 0) A = d->C[k];
//       else A = d->A[i][k];
//       if(A->type == 's')
//         printf("A(%d,%d) = %.2f   ", k, i, (double)200.0*A->sp->nnz/(d->blksz[k]*(d->blksz[k]+1)));
//       if(A->type == 'd')
//         printf("A(%d,%d) = %.2f   ", k, i, (double)100.0*A->diag->nnz/d->blksz[k]);
//     }
//     printf("\n");
//   }
//   */

//   //getchar();

//   return 0;
// }


int printheading(int start)
{
  if(start == 1) {
#ifdef __PC
    printf("\n");
    printf("                  ***   SDPLR %.2f   ***\n\n", VERSION);
    printf("=============================================================\n");
    printf(" major   minor         val         infeas       cgs    time  \n");
    printf("-------------------------------------------------------------\n");
#else
    printf("\n");
    printf("                 ***   SDPLR %.2f   ***\n\n", VERSION);
    printf("===========================================================\n");
    printf(" major   minor        val        infeas       cgs    time  \n");
    printf("-----------------------------------------------------------\n");
#endif
  }

  if(start == 0) {
#ifdef __PC
    printf("=============================================================\n\n");
#else
    printf("===========================================================\n\n");
#endif
  }

  fflush(stdout);

  return 1;
}



  

