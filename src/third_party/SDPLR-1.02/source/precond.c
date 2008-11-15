#include "sdplrlib.h"

int simple_getLiWWtLit(double* out, double* in)
{
  extern problemdata* global_data;

  // This routine is called by ARPACK, which assumes
  // out and in are numbered from 0

  // Other routines assume from 1

  EASYDCOPY(global_data->nr, in-1, out-1);
  block_solve_LT(global_data, out-1);

  // Note: this routine overwrites input vector
  getWWt(global_data, global_copyR, out-1, out-1);

  block_solve_L(global_data, out - 1);

  return 0;
}


int block_solve_L(problemdata* d, double* vec)
{
  int base, k, r;
  double one = 1.0;
  char side, uplo, trans, diag;

  base = 0;
  for(k = 1; k <= d->numblk; k++) {
    if(d->blktype[k] == 's' && d->XS_blksto[k] == 's') {
      multrhs_solve_L_vec(vec + base, d->L[k], d->Gcolptr[k], d->Gcolind[k], d->blksz[k], d->rank[k]);
    }
    else if(d->blktype[k] == 's' && d->XS_blksto[k] == 'd') {
      side = 'l'; uplo = 'l'; trans = 'n'; diag = 'n';
      dtrsm_(&side, &uplo, &trans, &diag, &(d->blksz[k]), &(d->rank[k]), &one,
             d->L[k] + 1, &(d->blksz[k]), vec + base + 1, &(d->blksz[k]));
    }
    else if(d->blktype[k] == 'd') {
      for(r = 1; r <= d->blksz[k]; r++) vec[base + r] /= d->L[k][r];
    }
    base += d->blksz[k]*d->rank[k];
  }

  return 0;
}

int block_solve_LT(problemdata* d, double* vec)
{
  int base, k, r;
  double one = 1.0;
  char side, uplo, trans, diag;

  base = 0;
  for(k = 1; k <= d->numblk; k++) {
    if(d->blktype[k] == 's' && d->XS_blksto[k] == 's') {
      multrhs_solve_LT_vec(vec + base, d->L[k], d->Gcolptr[k], d->Gcolind[k], d->blksz[k], d->rank[k]);
    }
    else if(d->blktype[k] == 's' && d->XS_blksto[k] == 'd') {
      side = 'l'; uplo = 'l'; trans = 't'; diag = 'n';
      dtrsm_(&side, &uplo, &trans, &diag, &(d->blksz[k]), &(d->rank[k]), &one,
             d->L[k] + 1, &(d->blksz[k]), vec + base + 1, &(d->blksz[k]));
    }
    else if(d->blktype[k] == 'd') {
      for(r = 1; r <= d->blksz[k]; r++) vec[base + r] /= d->L[k][r];
    }
    base += d->blksz[k]*d->rank[k];
  }

  return 0;
}




int multrhs_solve_L_vec(double* X, double* L, int* Lcolptr, int* Lcolind, int n, int m)
{
  int i, j;

  // X is n x m

  for(j = 1; j <= n-1; j++) {
    mydscal(m, 1.0/L[j], X + j, n);
    for(i = Lcolptr[j]; i <= Lcolptr[j+1]-1; i++)
      mydaxpy(m, -L[i], X + j, n, X + Lcolind[i], n);
  }

  mydscal(m, 1.0/L[n], X + n, n);
  
  return 0;
}


int multrhs_solve_LT_vec(double* X, double* L, int* Lcolptr, int* Lcolind, int n, int m)
{
  int i, j;

  mydscal(m, 1.0/L[n], X + n, n);

  for(j = n-1; j >= 1; j--) {
    for(i = Lcolptr[j]; i <= Lcolptr[j+1]-1; i++)
      mydaxpy(m, -L[i], X + Lcolind[i], n, X + j, n);
    mydscal(m, 1.0/L[j], X + j, n);
  }
  
  return 0;
}


int copy_sparse(double* new, double* old, int* Ecolptr, int* Econversion, int* EinG, int n)
{
  int i, j;

  // new is of the form G
  // old is of the form E

  for(j = 1; j <= n; j++) {

    new[j] = old[j];

    for(i = Ecolptr[j]; i <= Ecolptr[j+1]-1; i++)
      new[ EinG[i] ] = old[ Econversion[i] ];

  }

  return 0;

}

int convert_sparse_format(int* Ecolptr, int* Ecolind, int* Econversion, int* cptr, int* rind, int n)
{

  int i, j, k, ct;
  int **adj;

  MYCALLOC(adj, int*, n+1);
  for(i = 1; i <= n; i++) MYCALLOC(adj[i], int, n+1);

  for(i = 1; i <= n; i++) for(j = 1; j <= n; j++) adj[i][j] = 0;

  k = 1;
  for(j = 1; j <= n; j++)
    for(i = cptr[j]; i <= cptr[j+1]-1; i++)
      adj[ rind[i] ][ j ] = adj[ j ][ rind[i] ] = k++;

  for(i = 1; i <= n; i++)
    if(adj[i][i] == 0) {
      printf("Problem! Not all diag elements included!\n");
      exit(0);
    }

  Ecolptr[0] = 1; ct = n;
  for(k = 1; k <= n; k++) {
    Ecolind[k] = k;
    Econversion[k] = adj[k][k];
  }
  for(j = 1; j <= n; j++) {
    Ecolptr[j] = k = Ecolptr[j-1] + ct;
    ct = 0;
    for(i = j+1; i <= n; i++) if(adj[i][j] > 0) {
      ct++;
      Ecolind[k] = i;
      Econversion[k++] = adj[i][j];
    }
  }
  Ecolptr[n+1] = Ecolptr[n] + ct;

  for(i = 1; i <= n; i++) MYFREE(adj[i]);
  MYFREE(adj);

  return 0;
}

double incomplete_Cholesky(double* L, int* Lcolptr, int* Lcolind, int n)
  {

    int i, j, k;
    double *vec, veck;

    MYCALLOC(vec, double, n+1);
    if(vec == NULL) {
      printf("Memory allocation problem!\n");
      exit(0);
    }

    // Loop through columns

    for(j = 1; j <= n; j++)
      {

        // Calculate diagonal element
    
        if(L[j] > 0.0)
          {
            L[j] = sqrt( L[j] );
          }
        else
          {
            return L[j];
          }

        // Divide off-diag elts of j-th column by diag elt

        for(i = Lcolptr[j]; i <= Lcolptr[j+1]-1; i++)
          {
            L[i] /= L[j];
          }

        // Store j-th col in temp vec

        ZEROVEC(vec, n);
        // ZEROVEC(vec + j - 1, n - j + 1); // alternate

        vec[j] = L[j];

        for(i = Lcolptr[j]; i <= Lcolptr[j+1]-1; i++)
          {
            vec[ Lcolind[i] ] = L[i];
          }

        // Update remaining columns (j+1 .. n) [index = k]

        for(k = j+1; k <= n; k++)
          {
            veck = vec[k];

            L[k] -= veck * veck;

            for(i = Lcolptr[k]; i <= Lcolptr[k+1]-1; i++)
              {
                L[i] -= vec[ Lcolind[i] ] * veck;
              }
          }

   
      }

    MYFREE(vec);

    return 1.0;

  }



int makeG(problemdata* data, double density)
{
  int i, j, k, s, dim, ret, *xadj, *adjncy;
  int Ennz, *Ecolptr, *Ecolind, Fnnz, *Fcolptr, *Fcolind, *temp_Fcolind;
  int *EinF, numtodelete;
  int Gnnz, *Gcolptr, *Gcolind, *EinG;
  int *tiv1, *tiv2;

  for(k = 1; k <= data->numblk; k++) {

    if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == SPARSE) {

      // Save dimension of current block
      dim = data->blksz[k];
  
      // Save sparse E info for current block
      // Will point to only strictly lower (upper?) part
      // which differs from default format
      MYCALLOC(Ecolptr, int, dim+2);
      for(i = 1; i <= dim+1; i++) Ecolptr[i] = data->Ecolptr[k][i] - dim;
      Ennz = data->Ecolptr[k][dim+1]-1 - dim;
      MYCALLOC(Ecolind, int, Ennz + 1);
      for(i = 1; i <= Ennz; i++) Ecolind[i] = data->Ecolind[k][i + dim];

      if(Ennz > 0) {

        // tiv1 and 2 will store sorted rows and columns
        // of symmetric version of off diag of E
        MYCALLOC(tiv1, int, 2*Ennz + 1);
        MYCALLOC(tiv2, int, 2*Ennz + 1);
        for(j = 1; j <= dim; j++)
          for(i = Ecolptr[j]; i <= Ecolptr[j+1]-1; i++) {
            tiv1[i] = tiv2[i + Ennz] = Ecolind[i];
            tiv2[i] = tiv1[i + Ennz] = j;
          }
        quicksort2(tiv1, tiv2, 1, 2*Ennz);

        // Make xadj and adjncy for symbolic factorization procedure
        MYCALLOC(xadj, int, dim+2);
        MYCALLOC(adjncy, int, 2*Ennz+1);
        s = 1; xadj[1] = 1; xadj[dim+1] = 2*Ennz + 1;
        for(i = 1; i <= 2*Ennz; i++) {
          while(s != tiv1[i]) xadj[++s] = i;
          adjncy[i] = tiv2[i];
        }
        while(s != dim) xadj[++s] = 2*Ennz + 1;

        // Do symbolic factorization procedure, store in Fcolptr and Fcolind
        MYCALLOC(Fcolptr, int, dim+2); temp_Fcolind = NULL;
        ret = Ennz;
        while(ret > 0) {
          MYCALLOC(temp_Fcolind, int, 5*ret+1);
          ret = symbfactor(dim, Fcolptr, &Fnnz, temp_Fcolind, xadj, adjncy, 5*ret);
          if(ret > 0) MYFREE(temp_Fcolind);
        }
        MYCALLOC(Fcolind, int, Fnnz+1);
        for(i = 1; i <= Fnnz; i++) Fcolind[i] = temp_Fcolind[i];
        // free up leftovers
        MYFREE(temp_Fcolind);
        MYFREE(xadj);
        MYFREE(adjncy);
        MYFREE(tiv1);
        MYFREE(tiv2);

        //
        // Now have both E and F in similar structures
        // (no diagonal part)
        //

        // Next goal is to get G in similar structure

        // First, mark which entries in F are also in E
        MYCALLOC(EinF, int, Fnnz+1);
        for(i = 1; i <= Fnnz; i++) EinF[i] = 0;
        for(j = 1; j <= dim; j++) {
          s = Fcolptr[j];
          for(i = Ecolptr[j]; i <= Ecolptr[j+1]-1; i++) {
            while(Fcolind[s] != Ecolind[i]) s++;
            EinF[s] = 1;
          }
        }

        // Second, mark which entries will be deleted
        numtodelete = (int)(1.0 - density)*(Fnnz - Ennz);
        for(j = dim; j >= 1; j--)
          for(i = Fcolptr[j+1]-1; i >= Fcolptr[j]; i--)
            if(EinF[i] == 0 && numtodelete > 0) {
              EinF[i] = -1;
              numtodelete--;
            }


        // Third, make G from the information gathered
    
        // Count nnz in G
        Gnnz = 0;
        for(j = 1; j <= dim; j++)
          for(i = Fcolptr[j]; i <= Fcolptr[j+1]-1; i++)
            if(EinF[i] != -1)
              Gnnz++;

        // Allocate space for Gcolptr and Gcolind
        MYCALLOC(Gcolptr, int, dim+2);
        MYCALLOC(Gcolind, int, Gnnz+1);

        // Setup Gcolptr and Gcolind
        Gcolptr[1] = 1;
        for(j = 1; j <= dim; j++) {
          s = 0;
          for(i = Fcolptr[j]; i <= Fcolptr[j+1]-1; i++)
            if(EinF[i] != -1) {
              Gcolind[Gcolptr[j] + s] = Fcolind[i];
              s++;
            }
          Gcolptr[j+1] = Gcolptr[j] + s;
        }

        // Finally, make ptr from E to G
        MYCALLOC(EinG, int, Ennz+1);
        for(i = 1; i <= Ennz; i++) EinG[i] = 0;
        for(j = 1; j <= dim; j++) {
          s = Gcolptr[j];
          for(i = Ecolptr[j]; i <= Ecolptr[j+1]-1; i++) {
            while(Gcolind[s] != Ecolind[i]) s++;
            EinG[i] = s;
          }
        }

        // Now that we have G and EinG, we can free memory
        MYFREE(EinF);
        MYFREE(Fcolptr);
        MYFREE(Fcolind);
        MYFREE(Ecolptr);
        MYFREE(Ecolind);

        // Copy over G structures to official data

        MYCALLOC(data->Gcolptr[k], int, dim + 2);
        MYCALLOC(data->Gcolind[k], int, dim + Gnnz + 1);
        MYCALLOC(data->EinG[k], int, dim + Ennz + 1);
        for(j = 1; j <= dim+1; j++) data->Gcolptr[k][j] = Gcolptr[j] + dim;
        for(i = 1; i <= dim; i++) data->Gcolind[k][i] = i;
        for(i = 1; i <= Gnnz; i++) data->Gcolind[k][i + dim] = Gcolind[i];
        for(i = 1; i <= dim; i++) data->EinG[k][i] = i;
        for(i = 1; i <= Ennz; i++) data->EinG[k][i + dim] = EinG[i] + dim;
        data->Gnnz[k] = dim + Gnnz;

        MYFREE(EinG);
        MYFREE(Gcolptr);
        MYFREE(Gcolind);

      }
      else {

        MYCALLOC(data->Gcolptr[k], int, dim + 2);
        MYCALLOC(data->Gcolind[k], int, dim + 1);
        MYCALLOC(data->EinG[k], int, dim + 1);
        for(j = 1; j <= dim+1; j++) data->Gcolptr[k][j] = 1 + dim;
        for(i = 1; i <= dim; i++) data->Gcolind[k][i] = i;
        for(i = 1; i <= dim; i++) data->EinG[k][i] = i;
        data->Gnnz[k] = dim;

      }

    } // end if
    else if(data->blktype[k] == 's' && data->XS_blksto[k] == 'd') {
      data->Gnnz[k] = data->blksz[k]*data->blksz[k];
    }
    else if(data->blktype[k] == 'd') {
      data->Gnnz[k] = data->blksz[k];
    }

  } // end for

  return 0;

}

int symbfactor(int n, int* Fcolptr, int* Fnnz, int* temp_Fcolind, int* xadj, int* adjncy, int maxsub)
{
  int flag, *mrglnk, *rchlnk, *tempxnzsub, maxsublocal, *marker;
  int i, j, k, jstrt, jstop;
  int *tempcolind, *perm, *invp;

  mrglnk = (int*)calloc(n+1, sizeof(int));
  rchlnk = (int*)calloc(n+1, sizeof(int));
  tempxnzsub = (int*)calloc(n+2, sizeof(int));
  marker = (int*)calloc(n+1, sizeof(int));

  maxsublocal = maxsub;
  
  perm = (int*)calloc(n+1, sizeof(int));
  invp = (int*)calloc(n+1, sizeof(int));

  for(i = 1; i <= n; i++) perm[i] = invp[i] = i;

  /* Perm and invp needed for subroutine. */

  smbfct_2(n, xadj, adjncy, &maxsublocal, Fcolptr, tempxnzsub, temp_Fcolind, Fnnz, &flag, mrglnk, rchlnk, marker);

  if(flag > 0 || *Fnnz > maxsub) return mymax(flag, (*Fnnz)/4);

  tempcolind = (int*)calloc(*Fnnz+1, sizeof(int));

  /* Copy temp_L_colind into tempcolind, uncompressing.
     Tempcolind is temporary storage. Copy back into temp_Fcolind. */

  k = 1;
  for(i = 1; i <= n; i++) 
  {
    jstrt = tempxnzsub[i];
    jstop = jstrt + Fcolptr[i+1] - 1 - Fcolptr[i];
    for(j = jstrt; j <= jstop; j++)
      tempcolind[k++] = temp_Fcolind[j];
  }
  for(i = 1; i <= *Fnnz; i++) temp_Fcolind[i] = tempcolind[i];


  free(tempxnzsub);
  free(mrglnk);
  free(rchlnk);
  free(marker);
  free(perm);
  free(invp);
  free(tempcolind);

  return 0;

}


int smbfct_2(int neqns, int* xadj, int *adjncy, int* maxsub, 
			 int* xlnz, int* xnzsub, int* nzsub, int* maxlnz, 
			 int* flag, int* mrglnk, int* rchlnk, int* marker)
{
  int i, inz, j, jstop, jstrt, k, knz, kxsub, mrgk;
  int lmax, m, nabor, node, np1, nzbeg, nzend, rchm, mrkflg;

  nzbeg = 1; nzend = 0; xlnz[1] = 1;
  for(k = 1; k <= neqns; k++) mrglnk[k] = marker[k] = 0;

  np1 = neqns + 1;
  for(k = 1; k <= neqns; k++) {

    knz = 0; mrgk = mrglnk[k]; mrkflg = 0; marker[k] = k;
    if(mrgk != 0) marker[k] = marker[mrgk];
    xnzsub[k] = nzend;
    node = k;
    jstrt = xadj[node]; jstop = xadj[node+1] - 1;
    if(jstrt > jstop) goto L1500;

    rchlnk[k] = np1;
    for(j = jstrt; j <= jstop; j++) {
      nabor = adjncy[j]; nabor = nabor;
      if(nabor <= k) goto L300;
      rchm = k;
L200:
      m = rchm;
      rchm = rchlnk[m];
      if(rchm <= nabor) goto L200;
      knz++;
      rchlnk[m] = nabor;
      rchlnk[nabor] = rchm;
      if(marker[nabor] != marker[k]) mrkflg = 1;
L300:;
    }

    lmax = 0;
    if(mrkflg != 0 || mrgk == 0) goto L350;
    if(mrglnk[mrgk] != 0) goto L350;
    xnzsub[k] = xnzsub[mrgk] + 1;
    knz = xlnz[mrgk+1] - xlnz[mrgk] - 1;
    goto L1400;

L350:
    i = k;
L400:
    i = mrglnk[i];
    if(i == 0) goto L800;
    inz = xlnz[i+1] - xlnz[i] - 1;
    jstrt = xnzsub[i] +1; jstop = xnzsub[i] + inz;
    if(inz <= lmax) goto L500;
    lmax = inz;
    xnzsub[k] = jstrt;

L500:
    rchm = k;
    for(j = jstrt; j <= jstop; j++) {
      nabor = nzsub[j];
L600:
      m = rchm;
      rchm = rchlnk[m];
      if(rchm < nabor) goto L600;
      if(rchm == nabor) goto L700;
      knz++;
      rchlnk[m] = nabor;
      rchlnk[nabor] = rchm;
      rchm = nabor;
L700:;
    }
    goto L400;

L800:
    if(knz == lmax) goto L1400;
    if(nzbeg > nzend) goto L1200;
    i = rchlnk[k];
    for(jstrt = nzbeg; jstrt <= nzend; jstrt++) {
      if(nzsub[jstrt] - i < 0) goto L900;
      if(nzsub[jstrt] - i == 0) goto L1000;
      if(nzsub[jstrt] - i > 0) goto L1200;
L900:;
    }
    goto L1200;
L1000:
    xnzsub[k] = jstrt;
    for(j = jstrt; j <= nzend; j++) {
      if(nzsub[j] != i) goto L1200;
      i = rchlnk[i];
      if(i > neqns) goto L1400;
    }
    nzend = jstrt - 1;

L1200:
    nzbeg = nzend + 1;
    nzend = nzend + knz;
    if(nzend > *maxsub) goto L1600;
    i = k;
    for(j = nzbeg; j <= nzend; j++) {
      i = rchlnk[i];
      nzsub[j] = i;
      marker[i] = k;
    }
    xnzsub[k] = nzbeg;
    marker[k] = k;

L1400:
    if(knz <= 1) goto L1500;
    kxsub = xnzsub[k];
    i = nzsub[kxsub];
    mrglnk[k] = mrglnk[i];
    mrglnk[i] = k;
L1500:
    xlnz[k+1] = xlnz[k] + knz;
  }

  *maxlnz = xlnz[neqns] - 1;
  *maxsub = xnzsub[neqns];
  xnzsub[neqns+1] = xnzsub[neqns];
  *flag = 0;
  return 1;

L1600:
  *flag = nzend;
  return 0;

}


int quicksort2(int* A1, int* A2, int p, int r)
{
   int q;

   if(p < r) 
   {
      q = partition2(A1, A2, p, r);
      quicksort2(A1, A2, p, q);
      quicksort2(A1, A2, q+1, r);
   }

   return 1;
}

int partition2(int* A1, int* A2, int p, int r)
{
   int k, m, i, j, t1, t2;

   k = A1[p]; m = A2[p];
   i = p-1;
   j = r+1;
   while(i < j) {
      do j--;
      while(A1[j] > k || (A1[j] == k && A2[j] > m) );
      do i++;
      while(A1[i] < k || (A1[i] == k && A2[i] < m) );
      if(i < j) {
         t1 = A1[j]; t2 = A2[j];
         A1[j] = A1[i]; A2[j] = A2[i];
         A1[i] = t1; A2[i] = t2;
      }
      else return j;
   }

   return 0;
}


int quicksort2d(double* A1, int* A2, int p, int r)
{
   int q;

   if(p < r) 
   {
      q = partition2d(A1, A2, p, r);
      quicksort2d(A1, A2, p, q);
      quicksort2d(A1, A2, q+1, r);
   }

   return 1;
}

int partition2d(double* A1, int* A2, int p, int r)
{
   int m, i, j, t2;
   double k, t1;

   k = A1[p]; m = A2[p];
   i = p-1;
   j = r+1;
   while(i < j) {
      do j--;
      while(A1[j] < k || (A1[j] == k && A2[j] < m) );
      do i++;
      while(A1[i] > k || (A1[i] == k && A2[i] > m) );
      if(i < j) {
         t1 = A1[j]; t2 = A2[j];
         A1[j] = A1[i]; A2[j] = A2[i];
         A1[i] = t1; A2[i] = t2;
      }
      else return j;
   }

   return 0;
}


int precond_preprocess(problemdata* d)
{
  int i, j, k, ell, base1, base2;
  double tempval, *tempvec=NULL;

  double *diagWWt=NULL, *diagS, *ARnorm=NULL, *normsave=NULL;
  int minindex;

  // Allocate necessary storage for this subroutine.
  if(d->pc == PC_S_DIAG + PC_DIAG || d->pc == PC_S_CHOL + PC_DIAG)
    MYCALLOC(diagWWt, double, d->nr + 1);

  MYCALLOC(diagS, double, d->n + 1);

  if(d->pc == PC_S_DIAG + PC_LARGEVECS || d->pc == PC_S_CHOL + PC_LARGEVECS)
    MYCALLOC(ARnorm, double, d->m + 1);

  if(d->pc == PC_S_FACT + PC_DIAG || d->pc == PC_S_FACT + PC_LARGEVECS)
    MYCALLOC(tempvec, double, d->nr + 1);

  if(d->pc == PC_S_FACT + PC_LARGEVECS)
    MYCALLOC(normsave, double, d->U_ncol+1);

  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////
  /////////////////////////////////////////////////////

  if(d->pc == PC_S_DIAG + PC_DIAG) {

    // pc = 1 is a simple diagonal preconditioner

    // Get diag of WW^T
    getdiagWWt(d,diagWWt);

    // Get diag of S
    getdiagS(d,diagS);

    // Store diag(2*S + 4*sigma*WW^T)
    EASYDSCAL(d->n, 2.0, diagS);
    EASYDSCAL(d->nr, 4.0*d->sigma, diagWWt);
    base1 = base2 = 0;
    for(j = 1; j <= d->numblk; j++) {
      for(ell = 1; ell <= d->rank[j]; ell++) {
        for(k = 1; k <= d->blksz[j]; k++) {
          d->D[base2 + k] = diagS[base1 + k] + diagWWt[base2 + k];
          if(d->D[base2 + k] < 0.0) d->D[base2 + k] *= -1.0; // adjust if necessary
        }
        base2 += d->blksz[j];
      }
      base1 += d->blksz[j];
    }

    goto END_PC_PREPROCESS;

  }

  if(d->pc == PC_S_DIAG + PC_LARGEVECS || d->pc == PC_S_DIAG + PC_EIGENVECS) {

    // pc = 2 approximates S + WW^T by Diag(S) + low-rank(WW^T) ...
    // based on large-norm approximation to WW^T

    // Get diag S
    getdiagS(d,diagS);

    // Figure out approximation to W
    if(d->pc == PC_S_DIAG + PC_LARGEVECS) getlargevecs(d,ARnorm);
    else {
#ifdef __ARPACK
      if(d->pc == d->pc == PC_S_DIAG + PC_EIGENVECS) geteigenvecs(d,simple_getWWt);
#else
      printf("Error (precond_preprocess): ARPACK not linked in, so pc=3 not valid.\n");
      exit(0);
#endif
    }

    // Store pos def approx to diag S
    // Could we get numerical error if diag entry is too small?

    // Currently our approximation is the idenity I (b/c of numerical difficulties?)

    base1 = base2 = 0;
    for(j = 1; j <= d->numblk; j++) {
      for(ell = 1; ell <= d->rank[j]; ell++) {
        for(k = 1; k <= d->blksz[j]; k++) {
          d->D[base2 + k] = 1.0; //mymax(sqrt(fabs(diagS[base1 + k])), 1.0); // adjust if necessary
        }
        base2 += d->blksz[j];
      }
      base1 += d->blksz[j];
    }

    // In order to integrate pc = 3 into the SMW calculations of other preconditioners,
    // we must get a system like D^{1/2} (I + UU^T) D^{1/2}...
    // Adjust U by D^{-1/2}
    for(j = 1; j <= d->nr; j++) 
      mydscal(d->U_ncol, 1.0/d->D[j], d->U + j, d->nr);

    // Scale so that we can forget about constants
    EASYDSCAL(d->nr, sqrt(2.0), d->D);
    EASYDSCAL(d->U_ncol*d->nr, sqrt(2.0*d->sigma), d->U);

    // In case we ever need it, copy W over for saving. It is
    // overwritten in SMW_preprocess
    //EASYDCOPY(d->U_ncol*d->nr, d->U, d->U_save);
    
    SMW_preprocess(d);

    goto END_PC_PREPROCESS;
  }

  if(d->pc == PC_S_CHOL + PC_DIAG) {

    // Pc = 11 approximates S by LL^T and WW^T by diag(WW^T). The unique
    // approach is, however, that diag(WW^T) is incorporated into L

    // Get diag of WW^T
    getdiagWWt(d,diagWWt);

    // Get diag of S (used in getting Chol fact)
    getdiagS(d,diagS);

    // Before going onto calculation of inc Chol fact
    // of S, alter diagS by adding (scaled)diagWWt.
    EASYDSCAL(d->nr, 2.0*d->sigma, diagWWt);
    base1 = base2 = 0;
    for(j = 1; j <= d->numblk; j++) {
      for(ell = 1; ell <= d->rank[j]; ell++) {
        for(k = 1; k <= d->blksz[j]; k++)
          diagS[base1 + k] += diagWWt[base2 + k];
        base2 += d->blksz[j];
      }
      base1 += d->blksz[j];
    }

    // Get Chol factorization
    getL(d,diagS);

    // Scale so we can forget about constants
    for(k = 1; k <= d->numblk; k++) EASYDSCAL(d->Gnnz[k], sqrt(2.0), d->L[k]);

    goto END_PC_PREPROCESS;

  }

  if(d->pc == PC_S_CHOL + PC_LARGEVECS || d->pc == PC_S_CHOL + PC_EIGENVECS) {

    // pc = 12 approximates Hess by LL^T + WW^T, where UU^T is
    // a low-rank spectral approximation of WW^T

    // Get diag of S (used in Chol fact)
    getdiagS(d,diagS);

    // Get Chol fact
    getL(d,diagS);

    // Get spectral decomposition
    if(d->pc == PC_S_CHOL + PC_LARGEVECS) getlargevecs(d,ARnorm);
    else {
#ifdef __ARPACK
      geteigenvecs(d,simple_getWWt);
#else
      printf("Error (precond_preprocess): ARPACK not linked in, so pc=13 not valid.\n");
      exit(0);
#endif

    }

    // For SMW, need system of the form L (I + UU^T) L^T...
    // Calculate U
    for(i = 1; i <= d->U_ncol; i++) block_solve_L(d, d->U + (i-1)*d->nr);

    // Scale so we can forget about constants
    for(k = 1; k <= d->numblk; k++) EASYDSCAL(d->Gnnz[k], sqrt(2.0), d->L[k]);
    EASYDSCAL(d->U_ncol*d->nr, sqrt(2.0*d->sigma), d->U);

    // In case we ever need it, copy W over for safe keeping
    //EASYDCOPY(d->U_ncol*d->nr, d->U, d->U_save);

    // Prepare SMW
    SMW_preprocess(d);

    goto END_PC_PREPROCESS;
  }

  if(d->pc == PC_S_FACT + PC_DIAG) {

    // Pc = 21 first approximates S by LL^T, then approximates
    // I + UU^T, where U = L^{-1}W by a simple diag preconditioner.
    // But this is not so simple because calculating diag(UU^T)
    // is quite involved.

    // Get diag of S (used in Chol fact)
    getdiagS(d,diagS);

    // Get Chol fact
    getL(d,diagS);

    // Calculate diag of UU^T
    ZEROVEC(d->D, d->nr);
    for(i = 1; i <= d->m; i++) {
      getLiAiR(d,tempvec,i);
      for(j = 1; j <= d->nr; j++)  d->D[j] += pow(tempvec[j], 2.0);
    }

    // Scale so we can forget about constants
    EASYDSCAL(d->nr, 2.0*d->sigma, d->D);

    // Add in identity
    for(j = 1; j <= d->nr; j++) d->D[j] += 1.0;

    // Scale so we can forget about constants
    for(k = 1; k <= d->numblk; k++) EASYDSCAL(d->Gnnz[k], sqrt(2.0), d->L[k]);

    goto END_PC_PREPROCESS;

  }

  if(d->pc == PC_S_FACT + PC_LARGEVECS) {

    getdiagS(d,diagS);
    getL(d,diagS);

    // Get large-norm vecs of U
    minindex = 1;
    for(j = 1; j <= d->U_ncol; j++) normsave[j] = -1.0e10;

    for(i = 1; i <= d->m; i++) {
      getLiAiR(d,tempvec,i);
      tempval = pow(EASYDNRM2(d->nr, tempvec), 2.0);
      if(tempval > normsave[minindex]) {
        EASYDCOPY(d->nr, tempvec, d->U + (minindex - 1)*d->nr);
        normsave[minindex] = tempval;
        for(j = 1; j <= d->U_ncol; j++)
          if(normsave[j] < tempval) {
            tempval = normsave[j];
            minindex = j;
          }
      }
    }

    for(k = 1; k <= d->numblk; k++) EASYDSCAL(d->Gnnz[k], sqrt(2.0), d->L[k]);
    EASYDSCAL(d->U_ncol*d->nr, sqrt(2.0*d->sigma), d->U);
    //EASYDCOPY(d->U_ncol*d->nr, d->U, d->U_save);
    SMW_preprocess(d);

    goto END_PC_PREPROCESS;

  }

  if(d->pc == PC_S_FACT + PC_EIGENVECS) {

    getdiagS(d,diagS);
    getL(d,diagS);
#ifdef __ARPACK
    geteigenvecs(d,simple_getLiWWtLit);
#else
    printf("Error (precond_preprocess): ARPACK not linked in, so pc=23 not valid.\n");
    exit(0);
#endif
    for(k = 1; k <= d->numblk; k++) EASYDSCAL(d->Gnnz[k], sqrt(2.0), d->L[k]);
    EASYDSCAL(d->U_ncol*d->nr, sqrt(2.0*d->sigma), d->U);
    //EASYDCOPY(d->U_ncol*d->nr, d->U, d->U_save);
    SMW_preprocess(d);

    goto END_PC_PREPROCESS;

  }

END_PC_PREPROCESS:

  if(d->pc == PC_S_DIAG + PC_DIAG || d->pc == PC_S_CHOL + PC_DIAG)
    MYFREE(diagWWt);

  MYFREE(diagS);

  if(d->pc == PC_S_DIAG + PC_LARGEVECS || d->pc == PC_S_CHOL + PC_LARGEVECS)
    MYFREE(ARnorm);

  if(d->pc == PC_S_FACT + PC_DIAG || d->pc == PC_S_FACT + PC_LARGEVECS)
    MYFREE(tempvec);

  if(d->pc == PC_S_FACT + PC_LARGEVECS)
    MYFREE(normsave);


  return 0;
}

void initialize_pc(problemdata *data)
{
  int i, k;
  double tv;

  if(data->pc == PC_S_DIAG + PC_DIAG      ||
     data->pc == PC_S_DIAG + PC_LARGEVECS ||
     data->pc == PC_S_DIAG + PC_EIGENVECS ||
     data->pc == PC_S_FACT + PC_DIAG)
    MYCALLOC(data->D, double, data->nr + 1);

  if(data->pc > PC_S_CHOL) {

    MYCALLOC(data->Ecolptr, int*, data->numblk + 1);
    MYCALLOC(data->Ecolind, int*, data->numblk + 1);
    MYCALLOC(data->Econversion, int*, data->numblk + 1);
    for(k = 1; k <= data->numblk; k++) {
      if(data->blktype[k] == SDPBLK && data->XS_blksto[k] == SPARSE) {
        MYCALLOC(data->Ecolptr[k], int, data->blksz[k] + 2);
        MYCALLOC(data->Ecolind[k], int, data->XS_blkptr[k+1] - data->XS_blkptr[k] + 1);
        MYCALLOC(data->Econversion[k], int, data->XS_blkptr[k+1] - data->XS_blkptr[k] + 1);
        convert_sparse_format(data->Ecolptr[k], data->Ecolind[k], data->Econversion[k], data->XS_colptr[k], data->XS_rowind[k], data->blksz[k]);
      }
    }

    MYCALLOC(data->Gcolptr, int*, data->numblk + 1);
    MYCALLOC(data->Gcolind, int*, data->numblk + 1);
    MYCALLOC(data->Gnnz, int, data->numblk + 1);
    MYCALLOC(data->EinG, int*, data->numblk + 1);
    makeG(data, data->gdens);

    MYCALLOC(data->L, double*, data->numblk + 1);
    for(k = 1; k <= data->numblk; k++)
      MYCALLOC(data->L[k], double, data->Gnnz[k] + 1);

  }

  if(data->pc%10 == PC_LARGEVECS || data->pc%10 == PC_EIGENVECS) {
//   if(data->pc == PC_S_DIAG + PC_LARGEVECS ||
//      data->pc == PC_S_DIAG + PC_EIGENVECS ||
//      data->pc == PC_S_CHOL + PC_LARGEVECS ||
//      data->pc == PC_S_CHOL + PC_EIGENVECS ||
//      data->pc == PC_S_FACT + PC_LARGEVECS ||
//      data->pc == PC_S_FACT + PC_EIGENVECS) {
//   if(data->pc == 2 || data->pc == 3 || data->pc == 12 || data->pc == 13 || data->pc == 22 || data->pc == 23) {

    tv = (double)data->dthresh_dim;
    tv = 1000.0*pow(tv, 2.0)/(double)data->nr;

    data->U_ncol = mymax(1, mymin( (int)tv, (int)sqrt((double)data->m) ));

    MYCALLOC(data->U, double, data->U_ncol*data->nr + 1);
    MYCALLOC(data->U_vals, double, data->U_ncol + 1);
    MYCALLOC(data->U_symm, double, data->U_ncol*data->U_ncol + 1);
    //MYCALLOC(data->U_save, double, data->U_ncol*data->nr + 1);

    MYCALLOC(data->I, int, data->m + 1);
    for(i = 1; i <= data->m; i++) data->I[i] = i;
    data->szI = data->U_ncol;

  }

  if(data->pc == PC_LBFGS) {

    MYCALLOC(data->pcvecs, lbfgsvec, data->numbfgsvecs + 1);
    for(i = 1; i <= data->numbfgsvecs; i++) {
      MYCALLOC((data->pcvecs+i)->s, double, data->nr + 1);
      MYCALLOC((data->pcvecs+i)->y, double, data->nr + 1);
    }

    MYCALLOC(data->pcvecsnew, lbfgsvec, data->numbfgsvecs + 1);
    for(i = 1; i <= data->numbfgsvecs; i++) {
      MYCALLOC((data->pcvecsnew+i)->s, double, data->nr + 1);
      MYCALLOC((data->pcvecsnew+i)->y, double, data->nr + 1);
    }

  }

}

void deinitialize_pc(problemdata* data)
{
  int k;

//   if(data->pc == 1 || data->pc == 2 || data->pc == 3 || data->pc == 21)
  if(data->pc == PC_S_DIAG + PC_DIAG      ||
     data->pc == PC_S_DIAG + PC_LARGEVECS ||
     data->pc == PC_S_DIAG + PC_EIGENVECS ||
     data->pc == PC_S_FACT + PC_DIAG)
    MYFREE(data->D);

  // Deallocate structure involved with preconditioning
  if(data->pc > PC_S_CHOL) {
//   if(data->pc >= 11) {

    for(k = 1; k <= data->numblk; k++) {
      if(data->blktype[k] == 's' && data->XS_blksto[k] == 's') {

        MYFREE(data->Ecolptr[k]);
        MYFREE(data->Ecolind[k]);
        MYFREE(data->Econversion[k]);

        MYFREE(data->Gcolptr[k]);
        MYFREE(data->Gcolind[k]);
        MYFREE(data->EinG[k]);

      }
      MYFREE(data->L[k]);
    }
    MYFREE(data->Ecolptr);
    MYFREE(data->Ecolind);
    MYFREE(data->Econversion);

    MYFREE(data->Gcolptr);
    MYFREE(data->Gcolind);
    MYFREE(data->Gnnz);
    MYFREE(data->EinG);

    MYFREE(data->L);

  }

  if(data->pc%10 == PC_LARGEVECS || data->pc%10 == PC_EIGENVECS) {
//   if(data->pc == 2 || data->pc == 3 || data->pc == 12 || data->pc == 13 || data->pc == 22 || data->pc == 23) {
    MYFREE(data->U);
    //MYFREE(data->U_save);
    MYFREE(data->U_vals);
    MYFREE(data->U_symm);
    MYFREE(data->I);
  }

  if(data->pc == PC_LBFGS) {

    for(k = 1; k <= data->numbfgsvecs; k++) {
      MYFREE((data->pcvecs+k)->s);
      MYFREE((data->pcvecs+k)->y);
      MYFREE((data->pcvecsnew+k)->s);
      MYFREE((data->pcvecsnew+k)->y);
    }
    MYFREE(data->pcvecsnew);
    MYFREE(data->pcvecs);

  }
  
}


int getdiagS(problemdata* d, double* diagS)
{
  int base, i, j, k;
  int *colptr, *rowind;
  double *S;
  base = 0;

  ZEROVEC(diagS, d->n);
  for(k = 1; k <= d->numblk; k++) {
    S = d->S + d->XS_blkptr[k] - 1;
    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == SPARSE) {
      colptr = d->XS_colptr[k];
      rowind = d->XS_rowind[k];
      for(j = 1; j <= d->blksz[k]; j++)
        for(i = colptr[j]; i <= colptr[j+1]-1; i++) 
          if(rowind[i] == j) {
            diagS[base + j] = S[i];
            i = colptr[j+1]; // early quit
          }
    }
    else if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE) {
      for(j = 1; j <= d->blksz[k]; j++)
        diagS[base + j] = S[SMATIND(j,j,d->blksz[k])];
    }
    else if(d->blktype[k] == DIAGBLK) {
      for(j = 1; j <= d->blksz[k]; j++)
        diagS[base + j] = S[j];
    }
    base += d->blksz[k];
  }

  return 0;
}


int getdiagWWt(problemdata* d, double* diagWWt)
{
  int base, h, i, k, ell;

  ZEROVEC(diagWWt, d->nr);
  for(i = 1; i <= d->m; i++) {
    base = 0;
    for(k = 1; k <= d->numblk; k++) {
      for(h = 1; h <= d->ARind[i][k][0]; h++)
        for(ell = 1; ell <= d->rank[k]; ell++)
          diagWWt[base + (ell-1)*d->blksz[k] + d->ARind[i][k][h]] += pow(d->AR[i][k][(ell-1)*d->ARind[i][k][0] + h], 2.0);
      base += d->blksz[k]*d->rank[k];
    }
  }

  return 0;
}

int getL(problemdata* d, double* diagS)
{
  int base, k, ell;
  double total_diagelt, diagelt;
  int info, jobz, dimwork, maxblksz;
  char uplo='u';
  double *work=NULL, *evals=NULL, mineval; 

  // For all remaining pc'ing schemes >= 11,
  // do incomplete Chol factorization of S [for sparse blocks]
  // (based only on SDP blocks, not diagonal ones)

  // On the other hand, for dense blocks, do eval calculation
  // and resulting Cholesky factorization

  // Compute approximation LL^T for S based on G structure
  // (for all but one pc scheme)

  // Allocate storage that will be needed for calculation of
  // eigenvalues for S for dense blocks. Nothing is allocated
  // if there are no dense blocks.
  maxblksz = 0;
  for(k = 1; k <= d->numblk; k++)
    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE)
      maxblksz = mymax(maxblksz, d->blksz[k]);
  if(maxblksz > 0) {
    dimwork = mymax(1, 3*maxblksz - 1);
    MYCALLOC(work, double, dimwork + 1);
    MYCALLOC(evals, double, maxblksz + 1);
  }

  base = 0;

  for(k = 1; k <= d->numblk; k++) {

    if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == SPARSE) {

      // Approximate S
      total_diagelt = 1.0e10;
      for(ell = 1; ell <= d->blksz[k]; ell++)
        total_diagelt = mymin(total_diagelt, diagS[base + ell]);
      total_diagelt = fabs(mymin(0.0, total_diagelt));
      diagelt = 0.0;
      do {
        total_diagelt += fabs(diagelt) + 1.0; // (double)1.0/d->blksz[k];   //(double)1.0/d->n;
        ZEROVEC(d->L[k], d->Gnnz[k]);
        copy_sparse(d->L[k], d->S + d->XS_blkptr[k] - 1, d->Ecolptr[k], d->Econversion[k], d->EinG[k], d->blksz[k]);
        for(ell = 1; ell <= d->blksz[k]; ell++) d->L[k][ell] = total_diagelt + diagS[base + ell];
        diagelt = incomplete_Cholesky(d->L[k], d->Gcolptr[k], d->Gcolind[k], d->blksz[k]);
      } while(diagelt <= 0.0);

    }
    else if(d->blktype[k] == SDPBLK && d->XS_blksto[k] == DENSE) {

      // Compute minimum eigenvalue of S
      jobz = 'n'; uplo = 'l';
      EASYDCOPY(d->blksz[k]*d->blksz[k], d->S + d->XS_blkptr[k] - 1, d->L[k]);
      dsyev_(&jobz, &uplo, &(d->blksz[k]), d->L[k] + 1, &(d->blksz[k]),
             evals + 1, work + 1, &dimwork, &info);
      if(info != 0) {
        printf("Eigenvalue computation failed.\n");
        exit(0);
      }
      mineval = evals[1];


      // Alter the diag of S by |mineval| + eps
      // Note: We are writing S into L
      // Not sure about choice of eps
      EASYDCOPY(d->blksz[k]*d->blksz[k], d->S + d->XS_blkptr[k] - 1, d->L[k]);
      if(mineval <= 0.0) {
        for(ell = 1; ell <= d->blksz[k]; ell++)
          d->L[k][ SMATIND(ell, ell, d->blksz[k]) ] += fabs(mineval) + 1.0;
      }

      // Do Cholesky factorization
      uplo = 'l';
      dpotrf_(&uplo, &(d->blksz[k]), d->L[k] + 1, &(d->blksz[k]), &info);
      if(info != 0) {
        printf("Cholesky factorization failed.\n");
        exit(0);
      }

    }
    else if(d->blktype[k] == DIAGBLK) {

      for(ell = 1; ell <= d->blksz[k]; ell++) {
        // d->L[k][ell] = sqrt(mymax(d->S[k]->diag[ell], (double)1.0/d->n));
        // above line was causing trouble with inaccuracy/convergence // note old data structures
        // for pc=1 or pc=2
        d->L[k][ell] = 1.0;
      }

    }

    base += d->blksz[k];

  }

  if(maxblksz > 0) {
    MYFREE(work);
    MYFREE(evals);
  }

  return 0;

}

#ifdef __ARPACK

int geteigenvecs(problemdata* d, int (*matvec)(double*, double*))
{
  int j, ncv, maxitr, printlevel, nconv, nummatvec, ret;
  char which[3];

  global_data = d;
  ncv = mymin(mymax(2*d->U_ncol,10),d->nr);
  strcpy(which, "LA\0");
  maxitr = 300;
  printlevel = 1;

  // When finished, evals[1] gives first eigenvalue
  // and evecs[1] is start of eigenvectors
  ret = my_arpack(matvec, d->nr, d->U_ncol, ncv, which, maxitr, printlevel, d->U_vals+1, d->U+1, &nconv, &nummatvec);

  // Scale evecs by apppropriate evals
  for(j = 1; j <= d->U_ncol; j++) {
    if(fabs(d->U_vals[j]) < 1.0e-15) d->U_vals[j] = 0.0;
    //tempval = 2.0*d->sigma*d->U_vals[j]; // not needed at this moment // also deprecated b/c of change to sigma
    //tempval = sqrt(tempval/(1.0 + tempval)); // not needed at this moment
    EASYDSCAL(d->nr, sqrt(d->U_vals[j]), d->U + (j-1)*d->nr);
  }

  // Save W
  //EASYDCOPY(d->U_ncol*d->nr, d->U, d->U_save);

  return 0;
}

#endif

int getlargevecs(problemdata* d, double* ARnorm)
{
  int base1, base2, h, i, k, ell;
  double tempval;

  // Calculate the norms of AR (believe I have this right)
  for(i = 1; i <= d->m; i++) {

    tempval = 0.0;
    for(k = 1; k <= d->numblk; k++)
      tempval += pow(EASYDNRM2(d->ARind[i][k][0]*d->rank[k], d->AR[i][k]), 2.0);

    ARnorm[i] = sqrt(tempval);

  }

  // Sort in order of decreasing norm
  for(i = 1; i <= d->m; i++) d->I[i] = i;
  quicksort2d(ARnorm, d->I, 1, d->m);

  ZEROVEC(d->U, d->U_ncol*d->nr);
  base1 = 0;
  for(i = 1; i <= d->U_ncol; i++) {
    base2 = 0;
    for(k = 1; k <= d->numblk; k++) {
      for(h = 1; h <= d->ARind[ d->I[i] ][k][0]; h++)
        for(ell = 1; ell <= d->rank[k]; ell++)
          d->U[base1 + base2 + (ell-1)*d->blksz[k] + d->ARind[ d->I[i] ][k][h]] =
            d->AR[ d->I[i] ][k][(ell-1)*d->ARind[ d->I[i] ][k][0] + h];
      base2 += d->blksz[k]*d->rank[k];
    }
    base1 += d->nr;
  }

  // Save U
  //EASYDCOPY(d->U_ncol*d->nr, d->U, d->U_save);

  return 0;

}

int getLiAiR(problemdata* d, double* tempvec, int i)
{
  int base, h, k, ell;

  ZEROVEC(tempvec, d->nr);
  base = 0;
  for(k = 1; k <= d->numblk; k++) {
    for(h = 1; h <= d->ARind[i][k][0]; h++)
      for(ell = 1; ell <= d->rank[k]; ell++)
        tempvec[base + (ell-1)*d->blksz[k] + d->ARind[i][k][h]] = d->AR[i][k][(ell-1)*d->ARind[i][k][0] + h];
    base += d->blksz[k]*d->rank[k];
	}

  block_solve_L(d, tempvec);

  return 0;
}


int SMW_preprocess(problemdata* d)
{
  int i, j, ret, inc=1, inc2;
  double one=1.0, *tvec;
  char uplo='u', trans='t', diag='n';

  // As part of SMW, form I + W^T W
  uplo = 'l';
  ZEROVEC(d->U_symm, d->U_ncol*d->U_ncol);
  for(j = 1; j <= d->U_ncol; j++) d->U_symm[(j-1)*d->U_ncol + j] = 1.0;
  dsyrk_(&uplo, &trans, &(d->U_ncol), &(d->nr), &one, d->U+1, &(d->nr), &one, d->U_symm+1, &(d->U_ncol));
  for(j = 1; j <= d->U_ncol; j++)
    for(i = j+1; i <= d->U_ncol; i++)
      d->U_symm[(i-1)*d->U_ncol + j] = d->U_symm[(j-1)*d->U_ncol + i];

  //for(j = 1; j <= d->U_ncol; j++)
  //  for(i = 1; i <= d->U_ncol; i++)
  //    printf("symm(%d,%d) = %0.15f;\n", j, i, d->U_symm[(j-1)*d->U_ncol + i]);

  // Factor I + W^T W
  uplo = 'l';
  dpotrf_(&uplo, &(d->U_ncol), d->U_symm+1, &(d->U_ncol), &ret);

  if(ret != 0) {
    printf("Problem with factorization (dpotrf_)!!!\n");
    exit(0);
  }

  // Form auxilliary matrix, which will allow SMW,
  // stored in U
  MYCALLOC(tvec, double, d->U_ncol+1);
  inc2 = d->nr;
  for(j = 1; j <= d->nr; j++) {
    mydcopy(d->U_ncol, d->U + j, d->nr, tvec+1, 1);
    uplo = 'l'; trans = 'n'; diag = 'n';
    dtrsv_(&uplo, &trans, &diag, &(d->U_ncol), d->U_symm+1, &(d->U_ncol), tvec+1, &inc);
    mydcopy(d->U_ncol, tvec+1, 1, d->U + j, d->nr);
  }
  MYFREE(tvec);

  return 0;

}


int dopc(problemdata *d, double* z, double* tempvec)
{
  int j;

  if(d->pc == PC_S_DIAG + PC_DIAG) {

    for(j = 1; j <= d->nr; j++) z[j] /= d->D[j];

  }

  if(d->pc == PC_S_DIAG + PC_LARGEVECS || d->pc == PC_S_DIAG + PC_EIGENVECS) {

    for(j = 1; j <= d->nr; j++) z[j] /= d->D[j];
    do_smw(d,z,tempvec);
    for(j = 1; j <= d->nr; j++) z[j] /= d->D[j];

  }

  if(d->pc == PC_S_CHOL + PC_DIAG) {

    block_solve_L(d,z);
    block_solve_LT(d,z);

  }

  if(d->pc == PC_S_CHOL + PC_LARGEVECS ||
     d->pc == PC_S_CHOL + PC_EIGENVECS ||
     d->pc == PC_S_FACT + PC_LARGEVECS ||
     d->pc == PC_S_FACT + PC_EIGENVECS) {
//   if(d->pc == 12 || d->pc == 13 || d->pc == 22 || d->pc == 23) {

    block_solve_L(d,z);
    do_smw(d,z,tempvec);
    block_solve_LT(d,z);

  }

  if(d->pc == PC_S_FACT + PC_DIAG) {

    block_solve_L(d,z);
    for(j = 1; j <= d->nr; j++) z[j] /= d->D[j];
    block_solve_LT(d,z);

  }

  if(d->pc == PC_LBFGS) {
    dirlbfgs(d,d->pcvecs,z,z,global_oldest,d->numbfgsvecs,0);
  }

  if(d->pc == PC_LBFGS_NEW) {
    dirlbfgs(d,global_vecs,z,z,global_oldest,d->numbfgsvecs,0);
  }

  return 0;

}

int do_smw(problemdata* d, double* z, double* tempvec)
{
  int one=1;
  double one_d=1.0, minusone=-1.0, zero=0.0;
  char trans;

  trans = 't'; 
  dgemv_(&trans, &(d->nr), &(d->U_ncol), &one_d, d->U + 1,
         &(d->nr), z + 1, &one, &zero, tempvec + 1, &one);
  trans = 'n';
  dgemv_(&trans, &(d->nr), &(d->U_ncol), &minusone, d->U + 1,
         &(d->nr), tempvec + 1, &one, &one_d, z + 1, &one);

  return 0;

}



