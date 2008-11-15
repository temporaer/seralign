#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

// Parameters and macros just for this file
#define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
#define MYFREE(VAR) free(VAR)
#define mymin(A, B) ((A) < (B) ? (A) : (B))
#define DATABLOCKIND(DATA,BLOCK,NUMBLOCK) ((DATA+1)-1)*NUMBLOCK + BLOCK - 1

#define RANDOM 1

// Declaration of subroutines called in this file
int generate_params(void);
int getparams(char*, int*, double*, double*, int*, double*, int*, int*, int*, double*,
              int*, double*, int*, double*, int*, double*, int*, int*, double*, int*, int*);
int getstorage(int, int, int*, char*, double*, int*, int*, int*, char*, int*, int*, int*);
int writedata_raw(char*, int, int, int*, char*, double*, double*, int*, int*, int*, int*, char*, char*);
int writedata_sdplr(char*, int, int, int*, char*, double*, double*, int*, int*, int*, int*, char*, char*);
int readdata_raw(char*, int*, int*, int**, char**, double**, double**, int**, int**, int**, int**, char**, char**);
int readdata_sdplr(char*, int*, int*, int**, char**, double**, double**, int**, int**, int**, int**, char**, char**);
int readdata_sdpa(char*, int*, int*, int**, char**, double**, double**, int**, int**, int**, int**, char**, char**);
int sdplrlib(int, int, int*, char*, double*, double*, int*, int*, int*, int*, char*, char*,
             int, int, double, double, double, double, int, int, double, int, double, int, int, int, double, int, double, int, int,
             double*, double*, int*, int*, double*);
int readin(int, int, int*, char*, double*, double*, int*, int*, double*, FILE*);
int writeout(int, int, int*, char*, double*, double*, int*, int*, double*, FILE*);

int main(int argc, char *argv[])
{
  int h, k, ret, n, nr, *maxranks, *ranks;
  double *R, *lambda, pieces[8];
  FILE *fid;

  int m, numblk, *blksz;
  char *blktype;
  double *b, *CAent;
  int *CArow, *CAcol, *CAinfo_entptr, *CAinfo_rowcolptr;
  char *CAinfo_type, *CAinfo_storage;

  int method;
  int numbfgsvecs;
  double mixedtol;
  double rho_f;
  double rho_c;
  double sigmafac;
  int reorder;
  int precond;
  double gdens;
  int rankreduce;
  double gaptol;
  int checkbd;
  int typebd;
  int dthresh_dim;
  double dthresh_dens;
  int timelim;
  double rankredtol;
  int inputtype;
  int doAR;
  int printlevel;

  char *filein;
  char *fileout;

  // Begin: Perform some simple checks
  if(argc != 2 && argc != 3 && argc != 4 && argc != 5) {
    printf("Usage #1: %s <input_file> [params_file] [soln_in] [soln_out]\n", argv[0]);
    printf("Usage #2: %s gen_params\n", argv[0]);
    exit(0);
  }
  // End: Perform some simple checks

  if(argc == 2 && strcmp("gen_params", argv[1]) == 0) {

    generate_params();

    return 0;

  }

  if(argc == 2) {

    ret = getparams(NULL, &inputtype, &rho_f, &rho_c, &method,
                    &sigmafac, &rankreduce, &timelim, &dthresh_dim,
                    &dthresh_dens, &numbfgsvecs, &mixedtol, &precond,
                    &gdens, &reorder, &gaptol, &checkbd, &typebd, &rankredtol, &doAR, &printlevel);

  }

  if(argc != 2) {

    ret = getparams(argv[2], &inputtype, &rho_f, &rho_c, &method,
                    &sigmafac, &rankreduce, &timelim, &dthresh_dim,
                    &dthresh_dens, &numbfgsvecs, &mixedtol, &precond,
                    &gdens, &reorder, &gaptol, &checkbd, &typebd, &rankredtol, &doAR, &printlevel);

  }

  if(ret == 0) printf("Warning (main): Some problems getting parameters.\n");
  if(ret == -1) {
    printf("Error (main): Problem getting parameters.\n");
    exit(0);
  }


       if(argc == 4) { filein = argv[3]; fileout = NULL;    }
  else if(argc == 5) { filein = argv[3]; fileout = argv[4]; }
  else               { filein = NULL;    fileout = NULL;    }

  if(inputtype == 1) {

    readdata_sdpa(argv[1], &m, &numblk, &blksz, &blktype, &b, &CAent, &CArow, &CAcol,
                  &CAinfo_entptr, &CAinfo_rowcolptr, &CAinfo_type, &CAinfo_storage);

  }
  if(inputtype == 2) {

    readdata_sdplr(argv[1], &m, &numblk, &blksz, &blktype, &b, &CAent, &CArow, &CAcol,
                  &CAinfo_entptr, &CAinfo_rowcolptr, &CAinfo_type, &CAinfo_storage);

  }
  if(inputtype == 1000) {

     readdata_raw(argv[1], &m, &numblk, &blksz, &blktype, &b, &CAent, &CArow, &CAcol,
                  &CAinfo_entptr, &CAinfo_rowcolptr, &CAinfo_type, &CAinfo_storage);

  }

  // Quick error check

  for(k = 0; k < numblk; k++)
    if(blksz[k] == 0) {
      printf("Error (main): Block %d has size 0.\n", k);
      exit(0);
    }


  // Allocate space and setup info

  // First get dimensions and ranks
  MYCALLOC(maxranks, int, numblk);
  MYCALLOC(ranks, int, numblk);
  getstorage(m, numblk, blksz, blktype, CAent, CArow, CAcol, CAinfo_entptr, CAinfo_type, &n, &nr, maxranks);
  for(k = 0; k < numblk; k++) ranks[k] = maxranks[k];

  // Allocate space
  MYCALLOC(R, double, nr);
  MYCALLOC(lambda, double, m);

  // Setup initial point; either default or from file
  if(filein != NULL && (fid = fopen(filein, "r")) != NULL) {
    readin(m,numblk,blksz,blktype,R,lambda,maxranks,ranks,pieces,fid);
    fclose(fid);
  }
  else {
    if(RANDOM) srand( (unsigned)time( NULL ) );
    else       srand(925);
    for(h = 0; h < nr; h++) R[h] = (double)rand()/RAND_MAX - (double)rand()/RAND_MAX;
    pieces[0] = (double)0;
    pieces[1] = (double)0;
    pieces[2] = (double)0;
    pieces[3] = (double)0;
    pieces[4] = (double)0;
    pieces[5] = (double)0.0;
    pieces[6] = (double)1.0/n;
    pieces[7] = (double)1.0;
  }

//   writedata_sdpa("temp.dat-s", m, numblk, blksz, blktype, b, CAent, CArow, CAcol, CAinfo_entptr, CAinfo_rowcolptr, CAinfo_type, CAinfo_storage);
//   exit(0);
  
  sdplrlib(m, numblk, blksz, blktype, b, CAent, CArow, CAcol,
           CAinfo_entptr, CAinfo_rowcolptr, CAinfo_type, CAinfo_storage,
           method, numbfgsvecs, mixedtol, rho_f, rho_c, sigmafac, reorder, precond, gdens, rankreduce,
           gaptol, checkbd, typebd, dthresh_dim, dthresh_dens, timelim, rankredtol, doAR, printlevel,
           R-1, lambda-1, maxranks, ranks, pieces);

  if(fileout != NULL && (fid = fopen(fileout, "w")) != NULL) {
    writeout(m,numblk,blksz,blktype,R,lambda,maxranks,ranks,pieces,fid);
    fclose(fid);
  }

  MYFREE(R);
  MYFREE(lambda);
  MYFREE(maxranks);
  MYFREE(ranks);

  if(inputtype == 1 || inputtype == 2) {
    MYFREE(blksz);
    MYFREE(blktype);
    MYFREE(b);
    MYFREE(CAent);
    MYFREE(CArow);
    MYFREE(CAcol);
    MYFREE(CAinfo_entptr);
    MYFREE(CAinfo_rowcolptr);
    MYFREE(CAinfo_type);
    MYFREE(CAinfo_storage);
  }

  return 0;
}


int getstorage(int m, int numblk, int* blksz, char* blktype, double* CAent, int* CArow, int* CAcol,
               int* CAinfo_entptr, char* CAinfo_type, int* passedn, int* passednr, int* maxranks)
{
  int i, k, *cons, n, nr, ind;

  MYCALLOC(cons, int, m+1);

  n = 0;
  nr = 0;

  for(k = 1; k <= numblk; k++) {

    if(blktype[k-1] == 's') {

      for(i = 1; i <= m; i++) {

        ind = DATABLOCKIND(i,k,numblk);
        if(CAinfo_entptr[ind+1] > CAinfo_entptr[ind]) cons[i] = 1;
        else cons[i] = 0;

      }

      cons[0] = 0;
      for(i = 1; i <= m; i++) cons[0] += cons[i];
      // Note: this formula MUST match exactly the one found in initialize.c
      maxranks[k-1] = mymin((int)sqrt(2*cons[0]) + 1, blksz[k-1]);
      nr += blksz[k-1]*maxranks[k-1];
      n += blksz[k-1];

    }
    else if(blktype[k-1] == 'd') {
      maxranks[k-1] = 1;
      nr += blksz[k-1];
      n += blksz[k-1];
    }
  }

  *passedn  = n;
  *passednr = nr;

  MYFREE(cons);

  return 0;

}

int readin(int m, int numblk, int* blksz, char* blktype, double* R, double* lambda, int* maxranks, int* ranks, double* pieces, FILE* fid)
{
  int h, k, base, ti1, ti2, ti3;
  char tc1;

  fscanf(fid, "dual variable %d\n", &ti1);

  if(ti1 != m) {
    printf("Error (readin): Input solution and problem file do not match.\n");
    exit(0);
  }

  for(h = 0; h < m; h++) fscanf(fid, "%lf\n", &(lambda[h]));

  base = 0;
  for(k = 0; k < numblk; k++) {

    fscanf(fid, "primal variable %d %c %d %d %d\n", &ti1, &tc1, &ti2, &ti3, &(ranks[k])); ti1--;

    if(ti1 != k || tc1 != blktype[k] || ti2 != blksz[k] || ti3 != maxranks[k]) {
      printf("Error (readin): Input solution and problem file do not match.\n");
      exit(0);
    }

    for(h = 0; h < blksz[k]*ranks[k]; h++) fscanf(fid, "%lf\n", &(R[base + h]));
    base += blksz[k]*ranks[k];
  }

  fscanf(fid, "special majiter ");
  fscanf(fid, "%lf\n", &(pieces[0]));
  fscanf(fid, "special iter ");
  fscanf(fid, "%lf\n", &(pieces[1]));
  fscanf(fid, "special lambdaupdate ");
  fscanf(fid, "%lf\n", &(pieces[2]));
  fscanf(fid, "special CG ");
  fscanf(fid, "%lf\n", &(pieces[3]));
  fscanf(fid, "special curr_CG ");
  fscanf(fid, "%lf\n", &(pieces[4]));
  fscanf(fid, "special totaltime ");
  fscanf(fid, "%lf\n", &(pieces[5]));
  fscanf(fid, "special sigma ");
  fscanf(fid, "%lf\n", &(pieces[6]));
  fscanf(fid, "special scale ");
  fscanf(fid, "%lf\n", &(pieces[7]));

  return 0;
}

int writeout(int m, int numblk, int* blksz, char* blktype, double* R, double* lambda, int* maxranks, int* ranks, double* pieces, FILE* fid)
{
  int h, k, base;

  fprintf(fid, "dual variable %d\n", m);
  for(h = 0; h < m; h++) fprintf(fid, "%.16e\n", lambda[h]*pieces[7]);
  base = 0;
  for(k = 0; k < numblk; k++) {
    fprintf(fid, "primal variable %d %c %d %d %d\n", k+1, blktype[k], blksz[k], maxranks[k], ranks[k]);
    for(h = 0; h < blksz[k]*ranks[k]; h++)
      fprintf(fid, "%.16e\n", R[base + h]);
    base += blksz[k]*ranks[k];
  }
  fprintf(fid, "special majiter ");
  fprintf(fid, "%d\n", (int)pieces[0]);
  fprintf(fid, "special iter ");
  fprintf(fid, "%d\n", (int)pieces[1]);
  fprintf(fid, "special lambdaupdate ");
  fprintf(fid, "%d\n", (int)pieces[2]);
  fprintf(fid, "special CG ");
  fprintf(fid, "%d\n", (int)pieces[3]);
  fprintf(fid, "special curr_CG ");
  fprintf(fid, "%d\n", (int)pieces[4]);
  fprintf(fid, "special totaltime ");
  fprintf(fid, "%.16e\n", (double)pieces[5]);
  fprintf(fid, "special sigma ");
  fprintf(fid, "%.16e\n", (double)pieces[6]);
  fprintf(fid, "special scale ");
  fprintf(fid, "%.16e\n", (double)pieces[7]);

  return 0;
}


