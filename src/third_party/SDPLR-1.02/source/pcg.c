#include "sdplrlib.h"
#include "globals.h"

// We only need CGtype and R to get a direction. We assume that we
// already have relavent 'y' and 'sigma'.

int getPCG(problemdata *d, double* R, double* dR, double* x)
{
  /* if CGtype = 0, it's straight CG
        CGtype = 1, PCG with our preconditioning. */
  
  double accinit=0.0, acc, CG_tol=0.0, rtz, beta, oldrtz=0.0, ptHp, alpha;
  double *r, *rhs, *p, *z, *Hp, *temp;
  int CGS, auxCGS, iter, CG_maxiter;
  double *tempvec, *bestx;

  int first=1;
  double firstchange=0.0, val=0.0, oldval, change;

  double ss, nextval, bestnextval = 1.0e10;

  int i, oldestnew=1, ind1, ind2, newestnew;

  // Compute W
  if(d->doAR) AR(d,R);

  // Preprocess for preconditioning if applicable
  if(d->pc) precond_preprocess(d);


  // Begin CG (CG_tol will be set during first iteration of PCG)
  CG_maxiter = 100000;

  CGS  = 0; auxCGS = 0;
  iter = 0;

  r    = global_tempvec;
  rhs  = global_tempvec2;
  p    = global_tempvec3;
  z    = global_tempvec4;
  Hp   = global_tempvec5;
  temp = global_tempvec6;

  if(d->pc%10 == PC_LARGEVECS || d->pc%10 == PC_EIGENVECS) 
    MYCALLOC(tempvec, double, d->U_ncol + 1);
  else tempvec = NULL;

  MYCALLOC(bestx, double, d->nr+1);

  copyscaledvectovec(rhs,-1.0,dR,d->nr);
  ZEROVEC(x,d->nr);
  copyscaledvectovec(r,-1.0,dR,d->nr);

  while(iter <= CG_maxiter) {

    iter++;

    EASYDCOPY(d->nr,r,z);

    if(d->pc) dopc(d,z,tempvec);

//     if(iter == 1) {
//       EASYDCOPY(d->nr,z,x);
//       break;
//     }

    rtz = EASYDDOT(d->nr,r,z);

    if(iter == 1) {
      accinit = sqrt(rtz);
      CG_tol = mymin(0.1,accinit);
      acc = 1.0;
      EASYDCOPY(d->nr,z,p);
    }
    else {
      acc = sqrt(rtz/accinit); if(_isnan(acc)) break; // if(_isnan(acc)) { printf("Problem with acc.\n"); }
      beta = rtz/oldrtz;
      move_in_dir(p,z,beta,p,d->nr);
    }

    if(acc <= CG_tol) break;

    getMv(d, R, p, d->S, Hp, temp);
    CGS++;

    if(d->pc == PC_LBFGS) {
      EASYDCOPY(d->nr, x, (d->pcvecsnew + oldestnew)->s);
      EASYDCOPY(d->nr, r, (d->pcvecsnew + oldestnew)->y);
      oldestnew = oldestnew%d->numbfgsvecs + 1;
    }

    ptHp = EASYDDOT(d->nr,p,Hp);
    if(ptHp <= 0.0) {
      if(iter == 1) EASYDCOPY(d->nr,p,x);
      CGS = -CGS;
      break;
    }


    alpha = rtz/ptHp;
    EASYDAXPY(d->nr,alpha,p,x);
    if(iter%50 == 0) {
      getMv(d, R, x, d->S, Hp, temp);
      CGS++;
      move_in_dir(r,rhs,-1.0,Hp,d->nr);
    }
    else EASYDAXPY(d->nr,-alpha,Hp,r);

    oldrtz = rtz;

    if(iter%5 == 0 && EASYDDOT(d->nr,x,dR) < 0.0) {
      ss = linesearch(d,R,x,1.0,&nextval,0);
      if(first) {
        oldval = global_generic_value; val = nextval;
        firstchange = fabs(val - oldval)/(0.5*(fabs(oldval) + fabs(val)));
        first = 0;
      }
      else {
        oldval = val; val = nextval;
        change = fabs(val - oldval)/(0.5*(fabs(oldval) + fabs(val)));
        if(change < 0.1*firstchange || change < mymin(d->rho_f*d->rho_c,d->rho_f*d->rho_c/d->sigma)) iter = CG_maxiter + 1;
      }
      if(nextval > bestnextval) iter = CG_maxiter + 1;
      else {
        bestnextval = nextval;
        EASYDCOPY(d->nr,x,bestx);
      }
      if(iter == CG_maxiter + 1) EASYDCOPY(d->nr,bestx,x);
    }

  }

  if(d->pc%10 == PC_LARGEVECS || d->pc%10 == PC_EIGENVECS) 
    MYFREE(tempvec);

  if(d->pc == PC_LBFGS) {

    global_oldest = oldestnew;

    if(oldestnew > 1) newestnew = oldestnew - 1;
    else newestnew = d->numbfgsvecs;

    for(i = 1; i <= d->numbfgsvecs; i++) {

      ind1 = oldestnew + i - 1; if(ind1 > d->numbfgsvecs) ind1 = ind1%d->numbfgsvecs;
      ind2 = oldestnew + i;     if(ind2 > d->numbfgsvecs) ind2 = ind2%d->numbfgsvecs;

      if(ind1 != newestnew) {
        move_in_dir((d->pcvecs + ind1)->s, (d->pcvecsnew + ind2)->s, -1.0, (d->pcvecsnew + ind1)->s, d->nr);
        move_in_dir((d->pcvecs + ind1)->y, (d->pcvecsnew + ind2)->y, -1.0, (d->pcvecsnew + ind1)->y, d->nr);
      }
      else {
        move_in_dir((d->pcvecs + ind1)->s, x, -1.0, (d->pcvecsnew + ind1)->s, d->nr);
        move_in_dir((d->pcvecs + ind1)->y, r, -1.0, (d->pcvecsnew + ind1)->y, d->nr);
      }
      EASYDSCAL(d->nr,-1.0,(d->pcvecs + ind1)->y);
      (d->pcvecs + ind1)->rho = EASYDDOT(d->nr, (d->pcvecs + ind1)->y, (d->pcvecs + ind1)->s);
      if((d->pcvecs + ind1)->rho <= 0.0) {
        ZEROVEC((d->pcvecs + ind1)->s, d->nr);
        ZEROVEC((d->pcvecs + ind1)->y, d->nr);
        (d->pcvecs + ind1)->rho = 0.0;
      }
      else (d->pcvecs + ind1)->rho = 1.0/(d->pcvecs + ind1)->rho;


    }

  }

  MYFREE(bestx);

  return CGS;

}

int getMv(problemdata *d, double* R, double* vec, double *S, double* Mv, double* temp)
{
  StimesR(d, S, d->y, vec, Mv, 0);
  getWWt(d,R,temp,vec);
  EASYDAXPY(d->nr,2.0*d->sigma,temp,Mv);
  EASYDSCAL(d->nr,2.0,Mv);
  return 0;
}

int getWWt(problemdata* d, double* R, double* out, double* in)
{
  double *tv2;
  // Note there is no conflict if out = in

  // This assumes out and in are indexed from 1
  int h, i, k, ell, base;
  double *tv, tempval;
  int *ARindik, nnzr, blksz, rank;
  double *ARik;


  if(d->doAR == 0) {

    MYCALLOC(tv2, double, d->m + 1);
    ZEROVEC(out, d->nr);
    Aoper(d,R,in,UVt,0,0,tv2);
    AToper(d,tv2,UVt,0);
    StimesR(d,UVt,tv2,R,out,0);
    MYFREE(tv2);
    return 0;

  }


  MYCALLOC(tv, double, d->m + 1);

  for(i = 1; i <= d->m; i++) {

    base = 0;

    tempval = 0.0;

    for(k = 1; k <= d->numblk; k++) {

      ARik = d->AR[i][k];
      ARindik = d->ARind[i][k];
      nnzr = ARindik[0];
      blksz = d->blksz[k];
      rank = d->rank[k];

      for(h = 1; h <= nnzr; h++)
        for(ell = 1; ell <= rank; ell++)
          tempval += ARik[(ell-1)*nnzr + h]*in[base + (ell-1)*blksz + ARindik[h]];

      base += blksz*rank;

    }

    tv[i] = tempval;

  }

  ZEROVEC(out, d->nr);

  for(i = 1; i <= d->m; i++) {

    base = 0;

    for(k = 1; k <= d->numblk; k++) {

      ARik = d->AR[i][k];
      ARindik = d->ARind[i][k];
      nnzr = ARindik[0];
      blksz = d->blksz[k];
      rank = d->rank[k];

      for(h = 1; h <= nnzr; h++)
        mydaxpy(rank, tv[i], ARik + h, nnzr, out + base + ARindik[h], blksz);

      base += blksz*rank;

    }

  }

  MYFREE(tv);

  return 0;

}

int AR(problemdata* d, double* R)
{
  int h, i, k;
  int maxblksz, *inv;
  int base;
  int *ARindik, nnzr, rank, blksz;
  int row, col, ind;
  double *ARik, ent;

  maxblksz = -1;
  for(k = 1; k <= d->numblk; k++)
    maxblksz = mymax(maxblksz, d->blksz[k]);
  MYCALLOC(inv, int, maxblksz + 1);

  for(i = 1; i <= d->m; i++) {

    base = 0;

    for(k = 1; k <= d->numblk; k++) {

      ARindik = d->ARind[i][k];
      ARik = d->AR[i][k];
      rank = d->rank[k];
      blksz = d->blksz[k];
      nnzr = ARindik[0];

      ZEROVEC(ARik, nnzr*rank);

      for(h = 1; h <= blksz; h++) inv[h] = 0;
      for(h = 1; h <= nnzr; h++) inv[ ARindik[h] ] = h;

      if(d->blktype[k] == SDPBLK) {

        if(d->A[i][k]->type == 's') {

          for(h = 1; h <= d->A[i][k]->sp->nnz; h++) {

            ent = d->A[i][k]->sp->ent[h];
            row = d->A[i][k]->sp->row[h];
            col = d->A[i][k]->sp->col[h];

            mydaxpy(rank, ent, R + base + col, blksz, ARik + inv[row], nnzr);

            if(row != col)
              mydaxpy(rank, ent, R + base + row, blksz, ARik + inv[col], nnzr);
        
          }

        }

      }
      else if(d->blktype[k] == 'd') {

        for(h = 1; h <= d->A[i][k]->diag->nnz; h++) {
          ind = d->A[i][k]->diag->ind[h];
          ent = d->A[i][k]->diag->ent[h];
          ARik[ inv[ind] ] = ent*R[base + ind];
        }

      }

      base += d->blksz[k]*d->rank[k];

    }

  }

  MYFREE(inv);

  return 0;

}


int CG_alloc(problemdata *d)
{
  int h, i, k;
  int maxblksz, *cts;

  MYCALLOC(d->ARind, int**, d->m + 1);
  for(i = 1; i <= d->m; i++)
    MYCALLOC(d->ARind[i], int*, d->numblk + 1);

  MYCALLOC(d->AR, double**, d->m + 1);
  for(i = 1; i <= d->m; i++)
    MYCALLOC(d->AR[i], double*, d->numblk + 1);

  maxblksz = -1;
  for(k = 1; k <= d->numblk; k++)
    maxblksz = mymax(maxblksz, d->blksz[k]);

  MYCALLOC(cts, int, maxblksz + 1);

  for(i = 1; i <= d->m; i++)
    for(k = 1; k <= d->numblk; k++) {
      for(h = 1; h <= d->blksz[k]; h++) cts[h] = 0;
      if(d->blktype[k] == 's') {
        if(d->A[i][k]->type == 's') {
          for(h = 1; h <= d->A[i][k]->sp->nnz; h++) {
            cts[ d->A[i][k]->sp->row[h] ] = 1;
            cts[ d->A[i][k]->sp->col[h] ] = 1;
          }
        }
        else if(d->A[i][k]->type == 'l') {
          for(h = 1; h <= d->blksz[k]; h++) cts[h] = 1; // indicates completely dense
        }
      }
      else if(d->blktype[k] == 'd') {
        for(h = 1; h <= d->A[i][k]->diag->nnz; h++)
          cts[ d->A[i][k]->diag->ind[h] ] = 1;
      }
      cts[0] = 0;
      for(h = 1; h <= d->blksz[k]; h++) if(cts[h] == 1) cts[0]++;
      MYCALLOC(d->ARind[i][k], int, cts[0] + 1);
      MYCALLOC(d->AR[i][k], double, cts[0]*d->rank[k] + 1);
      d->ARind[i][k][0] = cts[0];
      //printf("%d %d  %f\n", i, k, (double)100*cts[0]/d->blksz[k]);
      cts[0] = 0;
      for(h = 1; h <= d->blksz[k]; h++) if(cts[h] == 1) {
        cts[0]++;
        d->ARind[i][k][cts[0]] = h;
      }
    }

  MYFREE(cts);

  return 0;

}

int CG_dealloc(problemdata *d)
{
  int i, k;

  for(i = 1; i <= d->m; i++)
    for(k = 1; k <= d->numblk; k++) {
      MYFREE(d->ARind[i][k]);
      MYFREE(d->AR[i][k]);
    }

  for(i = 1; i <= d->m; i++) MYFREE(d->AR[i]);
  MYFREE(d->AR);

  for(i = 1; i <= d->m; i++) MYFREE(d->ARind[i]);
  MYFREE(d->ARind);

  return 0;

}

