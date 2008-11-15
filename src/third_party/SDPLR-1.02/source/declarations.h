// copytstructures.c
int copystructures(problemdata*, int, int, int*, char*, double*, double*, int*, int*, int*, int*, char*, char*);

// dataoper.c
double function(problemdata* data, double* R);
int Aoper(problemdata* d, double* U, double* V, double* UVt, int same, int obj, double* results);
int Aoper_formUVt(problemdata* data, double* passedUVt, double* U, double* V, int same);
int gradient(problemdata* data, double* R);
int StimesR(problemdata *data, double *S, double *y, double *R, double *result, int save);
int Stimesmat(problemdata *data, double *S, double *y, double* vec, double* result, int n, int m, int k, int save);
int AToper(problemdata* data, double* y, double* S, int obj);
int sparsesymmmat_timesmat(problemdata *data, sparsesymmmat *S, double* vec, double* result, int n, int m);
double C_normdatamat(problemdata* data);
double normdatamat(problemdata* data, int matnum);
int essential_calcs(problemdata* data, double* R, double normC, double normb, double* val, double* rho_c_val, double* rho_f_val);

// eigen.c
int my_arpack(int (*)(double*, double*), int, int, int, char*, int, int, double*, double*, int*, int*);
int simple_getWWt(double*, double*);
int simple_Stimesvec_block(double*, double*);
int Smineval(problemdata*, double*);
int Sblockmineval(problemdata*, double*);
double Hmineval(problemdata*);
double Hmaxeval(problemdata*);
double dualbound(problemdata*);

// initialize.c
int initialize(problemdata*, int*);
int deinitialize(problemdata*);

// lbfgs.c
int dirlbfgs(problemdata*, lbfgsvec*, double*, double*, int, int, int);
int updatelbfgs1(problemdata*, lbfgsvec*, double*, int);
int updatelbfgs2(problemdata*, lbfgsvec*, double*, double*, double, int*, int);

// linesearch.c
double linesearch(problemdata*, double*, double*, double, double*, int);

// sdplrlib.c
int sdplrlib(int, int, int*, char*, double*, double*, int*, int*, int*, int*, char*, char*,
             int, int, double, double, double, double, int, int, double, int, double, int, int, int, double, int, double, int, int,
             double*, double*, int*, int*, double*);

// misc.c
int print_notes(void);
int printparams(problemdata*);
int print_dimacs_errors(problemdata*, double*);
int printstats(problemdata*);
int printheading(int);

// pcg.c
int getPCG(problemdata*, double*, double*, double*);
int AR(problemdata*, double*);
int getMv(problemdata *d, double* R, double* vec, double *S, double* Mv, double* temp);
int getWWt(problemdata*, double*, double*, double*);
int CG_alloc(problemdata*);
int CG_dealloc(problemdata*);

// precond.c
int solve_L_vec(double*, double*, int*, int*, int);
int solve_LT_vec(double*, double*, int*, int*, int);
int block_solve_L(problemdata*, double*);
int block_solve_LT(problemdata*, double*);
int multrhs_solve_L_vec(double*, double*, int*, int*, int, int);
int multrhs_solve_LT_vec(double*, double*, int*, int*, int, int);
int simple_getLiWWtLit(double* , double*);
int convert_sparse_format(int* Ecolptr, int* Ecolind, int* Econversion, int* cptr, int* rind, int n);
int copy_sparse();
double incomplete_Cholesky();
int makeG(problemdata*, double);
int symbfactor(int, int*, int*, int*, int*, int*, int);
int smbfct_2(int, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
int quicksort2(int*, int*, int, int);
int partition2(int*, int*, int, int);
int quicksort2d(double*, int*, int, int);
int partition2d(double*, int*, int, int);
int precond(problemdata*, lowrankmat**, lowrankmat**);
int precond_preprocess(problemdata*);
void initialize_pc(problemdata*);
void deinitialize_pc(problemdata*);
int getdiagS(problemdata*, double*);
int getdiagWWt(problemdata*, double*);
int getL(problemdata*, double*);
int geteigenvecs(problemdata*, int (*)(double*, double*));
int getlargevecs(problemdata*, double*);
int getLiAiR(problemdata*, double*, int);
int SMW_preprocess(problemdata*);
int dopc(problemdata*, double*, double*);
int do_smw(problemdata*, double*, double*);

// rankreduce.c
int dorankreduce(problemdata*, double*);

// timefuncs.c
#ifdef __PC
  double current_time(clock_t);
#else
  double current_time(double);
#endif

// util.c
int copyscaledvectovec(double*, double, double*, int);
int move_in_dir(double*, double*, double, double*, int);
int mydaxpy(int, double, double*, int, double*, int);
int mydcopy(int, double*, int, double*, int);
double myddot(int, double*, int, double*, int);
double mydnrm2(int, double*, int);
int mydscal(int, double, double*, int);
int createlowrankmat(problemdata*, lowrankmat**, int, int);
int destroylowrankmat(lowrankmat*);
int createsparsesymmmat(sparsesymmmat**, int);
int destroysparsesymmmat(sparsesymmmat*);
int createdatamat(problemdata*, datamat**, char, int, int, char*);
int destroydatamat(datamat*);

// external
extern void dgemm_();
extern void dgemv_();
extern void dgeqp3_();
extern void dsymv_();
extern void dsymm_();
extern void dsyr_();
extern void dsyrk_();
extern void dsyr2k_();
extern void dpotrf_();
extern void dtrsv_();
extern void dtrsm_();
extern void dsaupd_();
extern void dseupd_();
extern void dsyev_();
extern void dtrmm_();
extern int  idamax_();
extern void METIS_NodeND();
extern int gsl_poly_solve_cubic(double, double, double, double*, double*, double*);
extern double gsl_poly_eval(double*, int, double);

