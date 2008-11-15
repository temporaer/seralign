#include "sdplrlib.h"
#include "globalsmain.h"

#define PRINTFREQ      60.0

#define LAMBDAUPDATECT 1
#define SIGMASTRATEGY  1 // If 1, then LAMBDAUPDATECT should equal 1.

#define BEGIN 1
#define END   0

#define EASY   'e'
#define MEDIUM 'm'
#define HARD   'h'

void myprint(int majiter, int iter, double val, double rho_f_val, int CG, double totaltime);
int do_scaling(problemdata *data, double value, double *norm);

int sdplrlib (int m, int numblk, int *blksz, char *blktype, double *b,
              double *CAent, int *CArow, int *CAcol, int *CAinfo_entptr,
              int *CAinfo_rowcolptr, char *CAinfo_type, char
              *CAinfo_storage, int method, int numbfgsvecs, double
              mixedtol, double rho_f, double rho_c, double sigmafac,
              int reorder, int pc, double gdens, int rankreduce, double
              gaptol, int checkbd, int typebd, int dthresh_dim, double
              dthresh_dens, int timelim, double rankredtol, int doAR,
              int printlevel, double *R, double *lambda, int* maxranks,
              int *ranks, double *pieces)
{
  // Paramters that are stored in 'pieces'
  int    majiter, iter, CG, curr_CG, lambdaupdate;
  double timeoffset;

  // Algorithm declarations
  int          lambdaupdatect, oldest;
  double       val, rho_c_val, rho_f_val;
  double       alpha, rho_c_tol, normb, normC;
  double      *D;
  double       bestbd;
  lbfgsvec    *vecs=NULL;
  problemdata *data;

  // Various misc declarations
  int  i, k, ret, localiter=0, recalc=5, recalcct=10000000;
  char line[1000], difficulty='h';
  double bestinfeas=1.0e10;
  double sc, overallsc, tv; // scaling related
  double lastval; // monitor slowing down of optimization
  double origval;

  // Timing declarations and initialization
  double timeprintfreq = 0.0;
#ifdef __PC
  clock_t             timeorig;
  timeorig = clock();
#else
  double              timeorig;
  struct tms          timearray;
  times(&timearray);
  timeorig = timearray.tms_utime;
#endif
  timeoffset = 0.0;

  // Allocate space for main data
  MYCALLOC(data, problemdata, 1);

  // Copy input into data
  copystructures(data, m, numblk, blksz, blktype, b, CAent, CArow, CAcol,
                 CAinfo_entptr, CAinfo_rowcolptr, CAinfo_type, CAinfo_storage);

  // Set passed variables to data and global
  data->lambda = lambda;
  global_copyR = R;

  // Copy user params into data
  data->method       = method;
  data->numbfgsvecs  = numbfgsvecs;
  data->mixedtol     = mixedtol;
  data->rho_f        = rho_f;
  data->rho_c        = rho_c;
  data->sigmafac     = sigmafac;
  data->reorder      = reorder;
  data->pc           = pc;
  data->gdens        = gdens;
  data->rankreduce   = rankreduce;
  data->gaptol       = gaptol;
  data->checkbd      = checkbd;
  data->typebd       = typebd;
  data->dthresh_dim  = dthresh_dim;
  data->dthresh_dens = dthresh_dens;
  data->timelim      = timelim;
  data->rankredtol   = rankredtol;
  data->doAR         = doAR;
  data->printlevel   = printlevel;

  // Initialize all the data structures
  initialize(data, maxranks);


  /*** BEGIN: Create data structures needed in this function. (Note: Others created in initialize) ***/

  // Direction
  MYCALLOC (D, double, data->nr + 1);

  // Exact linesearch
  MYCALLOC(UVt, double, data->XS_blkptr[data->numblk+1]-1 + 1);
  MYCALLOC(ARD, double, data->m + 1);
  MYCALLOC(ADD, double, data->m + 1);

  // LBFGS
  if (data->method == 1 || data->method == 3 || (data->method == 2 && data->pc == PC_LBFGS_NEW) ) {
    MYCALLOC (vecs, lbfgsvec, data->numbfgsvecs + 1);
    for (i = 1; i <= data->numbfgsvecs; i++) {
      MYCALLOC ((vecs + i)->s, double, data->nr + 1);
      MYCALLOC ((vecs + i)->y, double, data->nr + 1);
    }
  }

  /*** END: Create data structures needed in this function. (Note: Others created in initialize) ***/


  // Setup algorithm parameters, quantities
  i = 1; normb          = fabs (data->b[idamax_ (&(data->m), data->b + 1, &i)]);
  normC                 = C_normdatamat (data);
  lambdaupdatect        = LAMBDAUPDATECT;
  oldest                = 1;
  data->method3_switch  = 0;
  bestbd                = -1.0e+20;

#ifdef __USER
  // This does not generalize well
  createlowrankmat (data, (lowrankmat **) & (data->user->QAP_SR), data->rank[1], data->blksz[1]);
#endif


  /*** Now really get started ***/
  /*** Now really get started ***/

  // printparams(data);
  if(data->printlevel > 0) printheading(BEGIN);

  // Setup ranks as read in from calling routine
  data->nr = 0;
  for (k = 1; k <= data->numblk; k++) {
    data->rank[k] = ranks[k - 1];
    data->nr += data->blksz[k] * data->rank[k];
  }

  // Get algorithm status parameters from calling routine
  majiter      = (int)    pieces[0];
  iter         = (int)    pieces[1];
  lambdaupdate = (int)    pieces[2];
  CG           = (int)    pieces[3];
  curr_CG      = (int)    pieces[4];
  timeoffset   = (double) pieces[5];
  data->sigma  = (double) pieces[6];
  overallsc    = (double) pieces[7];

  // Do scaling and inititalize total time correctly
  if(SCALE_OBJ && normC > 1.0e-10) do_scaling(data,overallsc,&normC);
  data->totaltime = timeoffset;

  // Now can setup rho_c_tol
  rho_c_tol = data->rho_c / data->sigma;

  // Calculate val, grad, etc.
  essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);

  // Save original value
  origval = overallsc*val;

  // If method=3, decide whether switch should be enforced
  if (data->method == 3 && data->method3_switch == 0 && rho_f_val < data->mixedtol)
    data->method3_switch = 1;

  /*** ITERATE ITERATE ITERATE ***/
  /*** ITERATE ITERATE ITERATE ***/
  /*** ITERATE ITERATE ITERATE ***/

  while (majiter++ < 100000) {

    while ( (!SIGMASTRATEGY && lambdaupdate < lambdaupdatect) || (SIGMASTRATEGY && difficulty != EASY) ) {

      // Increment lambda counter, reset local iter count and lastval
      lambdaupdate++; localiter = 0; lastval = 1.0e10;

      // Break if already meeting stopping criterion
      if (rho_c_val <= rho_c_tol) break;

      while(rho_c_val > rho_c_tol) { 

        if(data->method == 2 || (data->method == 3 && !data->method3_switch))
          if(localiter >= 11 && fabs(lastval - val)/(0.5*(fabs(lastval) + fabs(val))) <= mymin(data->rho_f*data->rho_c,data->rho_f*data->rho_c/data->sigma)) 
            break;

        // Increment both iter and localiter counts
        iter++; localiter++;

        // Direction calculation
        if(data->method == 1 && data->numbfgsvecs == 0) {
          copyscaledvectovec (D, -1.0, data->G, data->nr);
        }
        else if (data->method == 1 || (data->method == 3 && !data->method3_switch)) {
          dirlbfgs(data, vecs, D, data->G, oldest, data->numbfgsvecs, 1);
          updatelbfgs1(data, vecs, data->G, oldest);
        }
        else if (data->method == 2 || (data->method == 3 && data->method3_switch)) {
          global_generic_value = val; // What is this? Can't remember...
          if(data->pc == PC_LBFGS_NEW) {
            global_vecs = vecs;
            global_oldest = oldest;
          }
            ret = getPCG (data, R, data->G, D);
            if (ret <= 0) { CG += -ret; curr_CG += -ret; }
            else          { CG +=  ret; curr_CG +=  ret; }
          if(data->pc == PC_LBFGS_NEW) updatelbfgs1(data, vecs, data->G, oldest); // save LBFGS info
        }

        // If somehow we don't have descent, revert to steepest descent
        if (EASYDDOT (data->nr, D, data->G) >= 0.0)
          copyscaledvectovec (D, -1.0, data->G, data->nr);

        // Linesearch plus variable update
        lastval = val;
        alpha = linesearch (data, R, D, 1.0, &val, 1);
        EASYDAXPY (data->nr, alpha, D, R);

        // Refresh all the essentials
        if(recalc == 0) {
          essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
          recalc = recalcct;
        }
        else {
          gradient(data, R);
          rho_c_val = EASYDNRM2(data->nr, data->G)/(1.0 + normC);
          rho_f_val = EASYDNRM2(data->m, data->vio)/(1.0 + normb);
          recalc--;
        }

        // Direction calculation (necessary post-processing)
        if( (data->method == 1 && data->numbfgsvecs > 0)     ||
            (data->method == 3 && data->method3_switch == 0) ||
            (data->method == 2 && data->pc == PC_LBFGS_NEW) )
          updatelbfgs2 (data, vecs, D, data->G, alpha, &oldest, data->numbfgsvecs);

        // If PRINTFREQ seconds have passed since last major iteration, print an update
        timeprintfreq += current_time(timeorig) + timeoffset - data->totaltime;
        if(timeprintfreq > PRINTFREQ) {
          myprint(-1, iter, overallsc*val, rho_f_val, CG, data->totaltime);
          timeprintfreq -= PRINTFREQ;
        }

        // Update totaltime
        data->totaltime = current_time(timeorig) + timeoffset;

        // Possibly terminate
        if (data->totaltime >= data->timelim || rho_f_val <= data->rho_f || iter >= 1000000 || CG >= 1000000) {
            EASYDAXPY (data->m, -data->sigma, data->vio, data->lambda);
            essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
            goto END_CURR_MAJOR_ITERATION;
          }
        
        //printf("iter = %d  CGS = %4d  val = %f   ss = %.1e   \n", iter, CG, val, alpha);

        bestinfeas = mymin(rho_f_val, bestinfeas);
        if(data->method == 2 || (data->method == 3 && !data->method3_switch))
          if(rho_f_val/bestinfeas > 10.0)
            break;

      }

#ifdef __MEX
        printf(""); fflush(stdout);
        mexEvalString("drawnow;");
#endif


      // Update Lagrange multipliers and recalculate essentials
      EASYDAXPY (data->m, -data->sigma, data->vio, data->lambda);

      tv = EASYDNRM2(data->m,data->lambda);
      if(SCALE_OBJ && normC > 1.0e-10 && majiter >= 2 && (tv > 10.0 || tv < 0.1)) {
        if(tv > 10.0) sc = ( 1.0 - 0.9*pow(10.0/tv,0.1) )*tv;
        else          sc = ( 9.0*pow(10.0*tv,0.1) + 1.0 )*tv;
        overallsc *= sc;
        do_scaling(data,sc,&normC);
      }

      essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);

      // Decide how difficult this particular subproblem was
      if(SIGMASTRATEGY) {
          if(localiter <= 10)                           difficulty = EASY;
          else if(10+1 <= localiter && localiter <= 50) difficulty = MEDIUM;
          else if(50+1 <= localiter)                    difficulty = HARD;
      }

      if(data->method == 2 || (data->method == 3 && !data->method3_switch))
        if(rho_f_val/bestinfeas > 10.0) {
          bestinfeas *= sqrt(10.0);
          break;
        }

    } // end single major iteration

    // If user wants, calculate dual bound and save best so far.
    // Assumes only one block in SDP and dual shift is by identity
    // if (data->checkbd == 1) bestbd = mymax (bestbd, dualbound (data));
    
    // Check found grossly unbound value (indicative of infeasibility)
    if(overallsc*val > 1.0e10*fabs(origval)) {
      printf("Cannot reduce infeasibility any further.\n");
      goto END_MAJOR_ITERATIONS;
    }

  END_CURR_MAJOR_ITERATION:

    if(data->printlevel > 0) {
      if (data->checkbd == 1) {
        sprintf (line, "%3d %6d % .7e %.1e % .7e %5d    [ %5d  %5d ]\n", majiter, iter, val, rho_f_val, bestbd, (int) data->totaltime, curr_CG, CG);
        printf (line); fflush (stdout);
      }
      if (data->checkbd == 0 || data->checkbd == -1) myprint(majiter, iter, overallsc*val, rho_f_val, CG, data->totaltime);
    }

#ifdef __MEX
    mexEvalString("drawnow;");
#endif

    if (data->totaltime >= data->timelim || rho_f_val <= data->rho_f || iter >= 1000000 || CG >= 1000000)
      goto END_MAJOR_ITERATIONS;

    if (_isnan (val)) { printf ("Error(sdplrlib): Got NaN.\n"); exit (0); }

    /*** Get ready for next major iteration ***/

    // Adjust rank (hopefully reduce)
    if (data->rankreduce == 1) dorankreduce (data, R);

    // Update sigma
    do {
      data->sigma *= data->sigmafac;
      essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
      rho_c_tol = data->rho_c / data->sigma;
    } while (rho_c_val <= rho_c_tol);

    // Refresh some parameters
    lambdaupdate = 0;
    curr_CG = 0;
    if(SIGMASTRATEGY) difficulty = 'h';
    timeprintfreq = 0.0;

    // Decide whether to switch to CG method
    if (data->method == 3 && rho_f_val < data->mixedtol && data->method3_switch == 0)
      data->method3_switch = 1;

    // Clear BFGS vectors
    if (data->method == 1 || (data->method == 3 && data->method3_switch == 0) ||
        (data->method == 2 && data->pc == PC_LBFGS_NEW) )
      for (i = 1; i <= data->numbfgsvecs; i++) {
        ZEROVEC ((vecs + i)->s, data->nr);
        ZEROVEC ((vecs + i)->y, data->nr);
      }

    // Do some preconditioning operations
    if (data->pc == PC_LBFGS) {
      for (i = 1; i <= data->numbfgsvecs; i++) {
        ZEROVEC ((data->pcvecs + i)->s, data->nr);
        ZEROVEC ((data->pcvecs + i)->y, data->nr);
        ZEROVEC ((data->pcvecsnew + i)->s, data->nr);
        ZEROVEC ((data->pcvecsnew + i)->y, data->nr);
      }
    }

  } // end major iterations


END_MAJOR_ITERATIONS:

  // Final Lanczos
//   if (data->checkbd == 0) {
//     bestbd = mymax (bestbd, dualbound (data));
//     data->totaltime = current_time (timeorig) + timeoffset;
//     printf ("%.7e  %d\n", bestbd, (int) data->totaltime);
//   }

#ifdef __USER
  // This does not generalize well
  destroylowrankmat ((lowrankmat *) (data->user->QAP_SR));
  data->user->QAP_SR = NULL;
#endif

  if(data->printlevel > 0) {
    printheading(END);
    essential_calcs (data, R, normC, normb, &val, &rho_c_val, &rho_f_val);
    print_dimacs_errors (data, R);
  }

  // Return info to calling program
  pieces[0] = (double) majiter;
  pieces[1] = (double) iter;
  pieces[2] = (double) lambdaupdate;
  pieces[3] = (double) CG;
  pieces[4] = (double) curr_CG;
  pieces[5] = (double) data->totaltime;
  pieces[6] = (double) data->sigma;
  pieces[7] = (double) overallsc;
  for (k = 1; k <= data->numblk; k++)
    ranks[k - 1] = data->rank[k];


  // Clean up items allocated in this file
  MYFREE (UVt);
  MYFREE (ARD);
  MYFREE (ADD);
  if (data->method == 1 || data->method == 3 || (data->method == 2 && data->pc == PC_LBFGS_NEW)) {
    for (i = 1; i <= data->numbfgsvecs; i++) {
      MYFREE ((vecs + i)->s);
      MYFREE ((vecs + i)->y);
    }
    MYFREE (vecs);
  }
  MYFREE (D);
  deinitialize (data);
  MYFREE (data);

  return 0;
}

void myprint(int majiter, int iter, double val, double rho_f_val, int CG, double totaltime)
{
  char line[1000];
  if(majiter < 0) sprintf(line, "       %7d  % .8e  %.1e  %8d  %6d",          iter, val, rho_f_val, CG, (int)totaltime);
  else            sprintf(line, "  %3d  %7d  % .8e  %.1e  %8d  %6d", majiter, iter, val, rho_f_val, CG, (int)totaltime);
  printf(line);
  printf("\n");
  fflush(stdout);
}

int do_scaling(problemdata *data, double value, double *norm)
{
  int j, k;

  for(k = 1; k <= data->numblk; k++) {
    if(data->blktype[k] == 's') {
      if(data->C[k]->type == 's') {
        for(j = 1; j <= data->C[k]->sp->nnz; j++) {
          data->C[k]->sp->ent[j] /= value;
        }
      }
      else if(data->C[k]->type == 'l') 
        EASYDSCAL(data->C[k]->lr->ncol, 1.0/value, data->C[k]->lr->d);
    }
    else if(data->blktype[k] == 'd')
      for(j = 1; j <= data->C[k]->diag->nnz; j++) 
        data->C[k]->diag->ent[j] /= value;
  }

  for(j = data->AA_rowptr[0]; j <= data->AA_rowptr[1]-1; j++) {
    data->AA_colval_one[j] /= value;
    data->AA_colval_two[j] /= value;
  }

  *norm = C_normdatamat (data);

  EASYDSCAL(data->m, 1.0/value, data->lambda);

  return 0;

}
