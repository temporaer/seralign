typedef struct {
  double *s;
  double *y;
  double rho;
  double a;
} lbfgsvec;

typedef struct {
  double*  d;
  double*  ent;
  int      nrow;
  int      ncol;
} lowrankmat;

typedef struct {
  int*    row;
  int*    col;
  int     nnz;
  double* ent;
  int*    XS_in;
} sparsesymmmat;

typedef struct {
  int*    ind;
  int     nnz;
  double* ent;
  int*    XS_in;
} diagmat;

typedef struct {
  lowrankmat*    lr;
  sparsesymmmat* sp;
  diagmat*       diag;
  char           type;
  double*        multval;
  char*          label;
} datamat;

typedef struct {

  // user options
  int        method;
  int        numbfgsvecs;
  double     mixedtol;
  double     rho_f;
  double     rho_c;
  double     sigmafac;
  int        reorder;
  int        pc;
  double     gdens;
  int        rankreduce;
  double     gaptol;
  int        checkbd;
  int        typebd;
  int        dthresh_dim;
  double     dthresh_dens;
  int        timelim;
  double     rankredtol;
  int        doAR;
  int        printlevel;

  // very basic data (both SDPA and lowrank formats)
  int        m;
  int        numblk;
  int*       blksz;
  char*      blktype;
  datamat*** A;
  double*    b;
  datamat**  C;

  // auxiliary data, applies to lowrank format only
  double     evaladj;

  //
  // algorithm data, doesn't change once it is setup
  //
  // global dimension
  int        n;
  int        nr;
  // number of limited memory vecs to store
  int        M;
  // pointers to different types of matrices
  int**      lrind;
  int*       nlrind;
  int**      usind;
  int*       nusind;
  // rank information
  int*       rank;
  int*       maxrank;
  // norm of C for internal stopping criterion
  double     normC;
  // preconditioning data structures
  int** Ecolptr;
  int** Ecolind;
  int** Econversion;
  int** EinG;
  int** Gcolptr;
  int** Gcolind;
  int*  Gnnz;
  int   U_ncol;


  //
  // algorithm data, changes per iteration
  //
  // basic structures
  double*      lambda;
  double       sigma;
  double*      vio;
  double*      G;
  // CG structures
  double CG_tol;
  int    CG_maxiter;
  int    ***ARind;
  double ***AR;
  // preconditioning structures
  double *U;
  double *U_vals;
  double *U_symm;
  int *I;
  int szI;
  //double *U_save;
  double *D;
  double **L;
  int method3_switch;

  // timing data structures
  double totaltime;
  clock_t start;
  clock_t finish;

  // experimental
  //double *AAT;
  lbfgsvec *pcvecs;
  lbfgsvec *pcvecsnew;

  // user data structures
#ifdef __USER
  userdata* user;
#endif

  double   *S;
  double   *X;
  double   *y;
  int      *XS_blkptr;
  char     *XS_blksto;
  int     **XS_colptr;
  int     **XS_rowind;
  int      *AA_rowptr;
  int      *AA_colind;
  double   *AA_colval_one;
  double   *AA_colval_two;
  int      *lr_mat;
  int      *lr_blk;
  int       lr_num;

} problemdata;



