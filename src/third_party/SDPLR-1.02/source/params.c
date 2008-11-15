#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "version.h"

#define NUMPARAMS 20

// Default params (still hard coded in generate_params)
#define INPUTTYPE     1
#define RHO_F         1.0e-5
#define RHO_C         1.0e-1
#define METHOD        1
#define SIGMAFAC      2.0
#define RANKREDUCE    0
#define TIMELIM       3600
#define DTHRESH_DIM   10
#define DTHRESH_DENS  0.75
#define NUMBFGSVECS   4
#define MIXEDTOL      1.0e-5
#define PRECOND       0
#define GDENS         0.0
#define REORDER       0
#define GAPTOL        1.0e-3
#define CHECKBD       -1
#define TYPEBD        1
#define RANKREDTOL    2.2204460492503131e-16
#define DOAR          1
#define PRINTLEVEL    1

int generate_params_info(int);
int getparams_maxlinelength(FILE*);
int getparams_getline(FILE*, char*, int);
int getparams_tolower(char*, int);
int getparams_isnumber(char*);

// Macros just for this file
#define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
#define MYFREE(VAR) free(VAR)

int getparams(
char*  paramfile,
int*    inputtype,
double* rho_f,
double* rho_c,
int*    method,
double* sigmafac,
int*    rankreduce,
int*    timelim,
int*    dthresh_dim,
double* dthresh_dens,
int*    numbfgsvecs,
double* mixedtol,
int*    precond,
double* gdens,
int*    reorder,
double* gaptol,
int*    checkbd,
int*    typebd,
double* rankredtol,
int*    doAR,
int*    printlevel)
{
  int i, buffsz, ret, numparams=NUMPARAMS;
  int assigned[NUMPARAMS];
  double values[NUMPARAMS];
  char *buff, *match;
  FILE *fid;

  char paramstr[NUMPARAMS][50] = {
    "input type",
    "feasibility tolerance",
    "centrality tolerance",
    "step direction calculation",
    "penalty factor increase",
    "rank reduction",
    "time limit",
    "dimension threshold for dense matrices",
    "density threshold for dense matrices",
    "number of lbfgs vectors",
    "feasibility level for switch",
    "preconditioning method",
    "incomplete cholesky density",
    "reorder data to reduce fill-in",
    "duality gap tolerance",
    "check dual bound",
    "dual bound type",
    "rank tolerance",
    "auxiliary structures to facilitate tn",
    "print level"
  };

  // Assign default params just to be sure.
  *inputtype    = INPUTTYPE;         values[0]  = (double)INPUTTYPE; 
  *rho_f        = RHO_F;             values[1]  = (double)RHO_F;     
  *rho_c        = RHO_C;             values[2]  = (double)RHO_C;     
  *method       = METHOD;            values[3]  = (double)METHOD;    
  *sigmafac     = SIGMAFAC;          values[4]  = (double)SIGMAFAC;  
  *rankreduce   = RANKREDUCE;        values[5]  = (double)RANKREDUCE;
  *timelim      = TIMELIM;           values[6]  = (double)TIMELIM;   
  *dthresh_dim  = DTHRESH_DIM;       values[7]  = (double)DTHRESH_DIM;
  *dthresh_dens = DTHRESH_DENS;      values[8]  = (double)DTHRESH_DENS;
  *numbfgsvecs  = NUMBFGSVECS;       values[9]  = (double)NUMBFGSVECS;
  *mixedtol     = MIXEDTOL;          values[10] = (double)MIXEDTOL;  
  *precond      = PRECOND;           values[11] = (double)PRECOND;
  *gdens        = GDENS;             values[12] = (double)GDENS;
  *reorder      = REORDER;           values[13] = (double)REORDER;
  *gaptol       = GAPTOL;            values[14] = (double)GAPTOL; 
  *checkbd      = CHECKBD;           values[15] = (double)CHECKBD;
  *typebd       = TYPEBD;            values[16] = (double)TYPEBD;
  *rankredtol   = RANKREDTOL;        values[17] = (double)RANKREDTOL;
  *doAR         = DOAR;              values[18] = (double)DOAR;
  *printlevel   = PRINTLEVEL;        values[19] = (double)PRINTLEVEL;


  // Simple case (param file = NULL)
  if(paramfile == NULL) return 1;

  // Keep track of what has been assigned by file.
  for(i = 0; i < numparams; i++) assigned[i] = 0;

  // Open file for reading
  fid = fopen(paramfile, "r");

  if(fid == NULL) {
    printf("Warning (getparams): File %s not found. Using default parameters.\n", paramfile);
    return 0;
  }
  else if(fid != NULL) {

    buffsz = getparams_maxlinelength(fid) + 10;
    fclose(fid);

    fid = fopen(paramfile, "r");

    MYCALLOC(buff, char, buffsz);

    do {
      ret = getparams_getline(fid, buff, buffsz);
      getparams_tolower(buff, buffsz);

      for(i = 0; i < numparams; i++) {
        match = strstr(buff, paramstr[i]);
        if(match != NULL) {
          if(assigned[i] == 0) {
            match = strchr(buff, ':');
            if(match == NULL) {
              printf("Error (getparams): Parameter file has wrong format.\n");
              return -1;
            }
            match++;
            if(getparams_isnumber(match) != 1) {
              printf("Error (getparams): Parameter file has wrong format.\n");
              return -1;
            }
            values[i] = atof(match);
            assigned[i] = 1;
          }
          else if(assigned[i] == 1)
            printf("Warning (getparams): Attempt to assign parameter '%s' twice.\n", paramstr[i]);
        }
      }

    } while(ret != 0);

    MYFREE(buff);

    fclose(fid);

    // If some values have not been assigned, print warning
    for(i = 0; i < numparams; i++)
      if(assigned[i] == 0)
        printf("Warning (getparams): Some parameters not assigned. Using default values.\n");

    // assign values read from file
    // *INDENT-OFF*
    *inputtype    = (int)values[0];
    *rho_f        =      values[1];
    *rho_c        =      values[2];
    *method       = (int)values[3];
    *sigmafac     =      values[4];
    *rankreduce   = (int)values[5];
    *timelim      = (int)values[6];
    *dthresh_dim  = (int)values[7];
    *dthresh_dens =      values[8];
    *numbfgsvecs  = (int)values[9];
    *mixedtol     =      values[10];
    *precond      = (int)values[11];
    *gdens        =      values[12];
    *reorder      = (int)values[13];
    *gaptol       =      values[14];
    *checkbd      = (int)values[15];
    *typebd       = (int)values[16];
    *rankredtol   =      values[17];
    *doAR         = (int)values[18];
    *printlevel   = (int)values[19];
    // *INDENT-ON*

    // Perform some simple checks

    // Singleton checks

    if(*inputtype != 1 && *inputtype != 2 && *inputtype != 1000) {
      printf("Error (params): Parameter '%s' must be 1 or 2.\n", paramstr[0]);
      return -1;
    }

    if(*rho_f <= 0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[1]);
      return -1;
    }

    if(*rho_c <= 0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[2]);
      return -1;
    }


    if(*method != 1 && *method != 2 && *method != 3) {
      printf("Error (params): Parameter '%s' must be 1, 2, or 3.\n", paramstr[3]);
      return -1;
    }

    if(*sigmafac <= 1.0) {
      printf("Error (params): Parameter '%s' must be greater than 1.0.\n", paramstr[4]);
      return -1;
    }

    if(*rankreduce != 0 && *rankreduce != 1) {
      printf("Error (params): Parameter '%s' must be 0 or 1.\n", paramstr[5]);
      return -1;
    }

    if(*timelim <= 0) {
      printf("Parameter '%s' must be positive.\n", paramstr[6]);
      return -1;
    }

    if(*dthresh_dim < 0) {
      printf("Parameter '%s' must be nonnegative.\n", paramstr[7]);
      return -1;
    }

    if(*dthresh_dens < 0.0 || *dthresh_dens > 1.0) {
      printf("Parameter '%s' must be in [0,1].\n", paramstr[8]);
      return -1;
    }

    if(*numbfgsvecs < 0) {
      printf("Error (params): Parameter '%s' must be a non-negative integer.\n", paramstr[9]);
      return -1;
    }

    if(*mixedtol <= 0.0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[10]);
      return -1;
    }

    if(*precond != 0  &&
       *precond != 1  && *precond != 2  && *precond !=  3 &&
       *precond != 11 && *precond != 12 && *precond != 13 &&
       *precond != 21 && *precond != 22 && *precond != 23 &&
       *precond != 100 && *precond != 200) {
      printf("Error (params): Parameter '%s' has an illegal value.\n", paramstr[11]);
      return -1;
    }

#ifndef __ARPACK
    if(*precond == 3 || *precond == 13 || *precond == 23) {
      printf("Error (params): ARPACK not linked in, so paramater '%s' cannot be 3, 13, or 23.\n", paramstr[11]);
      return -1;
    }
#endif

    if(*gdens < 0.0 || *gdens > 1.0) {
      printf("Parameter '%s' must be in [0,1].\n", paramstr[12]);
      return -1;
    }

    if(*reorder != 0 && *reorder != 1) {
      printf("Error (params): Parameter '%s' must be 0 or 1.\n", paramstr[13]);
      return -1;
    }

#ifndef __METIS
    if(*reorder == 1) {
      printf("Error (params): METIS not linked in, so parameter '%s' must be 0.\n", paramstr[13]);
      return -1;
    }
#endif

    if(*gaptol <= 0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[14]);
      return -1;
    }

    // Hopefully temporary
    if(*checkbd != -1) {
      printf("Error (params): At this time, parameter '%s' must be -1.\n", paramstr[15]);
      return -1;
    }

    if(*checkbd != -1 && *checkbd != 0 && *checkbd != 1) {
      printf("Error (params): Parameter '%s' must be -1, 0, or 1.\n", paramstr[15]);
      return -1;
    }

    if(*typebd != 1) {
      printf("Error (params): Currently, parameter '%s' must equal 1.\n", paramstr[16]);
      return -1;
    }

    if(*rankredtol <= 0.0) {
      printf("Error (params): Parameter '%s' must be positive.\n", paramstr[17]);
      return -1;
    }

    if(*doAR != 0 && *doAR != 1) {
      printf("Error (params): Parameter '%s' must be 0 or 1.\n", paramstr[18]);
      return -1;
    }

    if(*printlevel != 0 && *printlevel != 1) {
      printf("Error (params): Parameter '%s' must be 0 or 1.\n", paramstr[19]);
      return -1;
    }

    // Pairwise checks

    if(*doAR == 0 && (*precond == 1 || *precond == 11)) {
      printf("Error (params): Incompatible '%s' and '%s' parameters.\n", paramstr[18], paramstr[11]);
      return -1;
    }


    return 1;

  }

  return -1;

}

int getparams_maxlinelength(FILE *datafile)
{
  int maxlen=0, k, c;
  
  do {

    k = 0;

    do {
      c = getc(datafile); 
      k++;
    } while(c != '\n' && c != EOF);

    if(k > maxlen) maxlen=k;

  } while(c != EOF);
  
  return(maxlen);
}

int getparams_getline(FILE *datafile, char *buffer, int bufsiz)
{
  int k;
  char c;
  
  k = 0;

  do {
    if(k >= bufsiz) {
      printf("Error (getparams_getline): Line too long!  Adjust bufsiz.\n");
      return(-1);
    }
    c = getc(datafile);
    buffer[k] = c;
    k++;
  } while(c != '\n' && c != EOF);

  if(c != '\n') {
    buffer[k] = '\n';
    buffer[k+1] = 0;
  }
  else buffer[k] = 0;

  if(c == EOF) return 0;
  else return 1;
}

int getparams_tolower(char* buff, int buffsz)
{
  int i;

  for(i = 0; i < buffsz; i++)
    buff[i] = tolower(buff[i]);

  return 1;
}

int getparams_isnumber(char* str)
{
  int i, length;

  length = strlen(str);

  for(i = 0; i < length; i++)
    if(isdigit(str[i]) == 0 && str[i] != '.' &&
       str[i] != '-' && str[i] != 'e' && isspace(str[i]) == 0 &&
       str[i] != '\n' && str[i] != '\0' && str[i] != EOF && str[i] != '+') {
      return 0;
    }

  return 1;
}

int generate_params(void)
{
  /*
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
  */

  int i;
  char in[NUMPARAMS][110], *ret, filename[] = "sdplr.params", usrfilename[100];
  FILE *fid;
  char paramstr[NUMPARAMS][100] = {
    "Input type (1=SDPA, 2=SDPLR)",
    "Feasibility tolerance",
    "Centrality tolerance",
    "Step direction calculation (1=LBFGS, 2=TN, 3=switch)",
    "Penalty factor increase",
    "Rank reduction (1=yes, 0=no)",
    "Time limit (in seconds)",
    "Print level (0=nothing, 1=all)",
    "Dimension threshold for dense matrices",
    "Density threshold for dense matrices (in [0,1])",
    "Number of LBFGS vectors",
    "Feasibility level for switch",
    "Auxiliary structures to facilitate TN (1=yes, 0=no)",
    "Rank tolerance (e.g., machine precision)",
    "Preconditioning method (0,1,2,3,11,12,13,21,22,23,100)",
    "Incomplete Cholesky density (in [0,1])",
    "Reorder data to reduce fill-in (1=yes, 0=no)",
    "Duality gap tolerance",
    "Check dual bound (-1, 0, 1)",
    "Dual bound type (1=+/-1, 2=0/1)"
  };
  char def[NUMPARAMS][50] = {
    "1",
    "1.0e-5",
    "1.0e-1",
    "1",
    "2.0",
    "0",
    "3600",
    "1",
    "10",
    "0.75",
    "4",
    "1.0e-5",
    "1",
    "2.2204460492503131e-16",
    "0",
    "0.0",
    "0",
    "1.0e-3",
    "-1",
    "1"
  };

  printf("\nSDPLR %.2f  --  Automatic Paramater File Generation\n\n", VERSION);

  do {
    printf("\n");
    printf("Parameter file name [%s]: ", filename); fflush(stdout);
    ret = fgets(usrfilename, 100, stdin);
    if(ret == NULL) {
      printf("Error\n");
      exit(0);
    }
    usrfilename[strlen(usrfilename)-1] = '\0';
    if(strlen(usrfilename) == 0) strcpy(usrfilename, filename);
    fid = fopen(usrfilename, "w");
  } while(fid == NULL);

  printf("\n\nPress 'i' for information at any time.\n");
  printf("Press 'i' for information at any time.\n");
  printf("Press 'i' for information at any time.\n\n"); fflush(stdout);

  for(i = 0; i < NUMPARAMS; i++) {
    do {
      printf("\n");
      printf("%s [%s]: ", paramstr[i], def[i]); fflush(stdout);
      ret = fgets(in[i], 100, stdin);
      if(ret == NULL) {
        printf("Error\n");
        exit(0);
      }
      in[i][strlen(in[i])-1] = '\0';
      if(strlen(in[i]) == 0) strcpy(in[i], def[i]);
      if(in[i][0] == 'i' || in[i][0] == 'I') generate_params_info(i);
    } while(getparams_isnumber(in[i]) == 0);
  }

  fprintf(fid, "SDPLR %.2f paramter file (automatically generated)\n\n", VERSION);
  fprintf(fid, "--> Basic parameters <--\n\n");
  for(i = 0; i < 10; i++)
    fprintf(fid, "%s : %s\n", paramstr[i], in[i]);
  fprintf(fid, "\n--> Conditional parameters <--\n\n");
  for(i = 10; i < 14; i++)
    fprintf(fid, "%s : %s\n", paramstr[i], in[i]);
  fprintf(fid, "\n--> In-development parameters <--\n\n");
  for(i = 14; i < NUMPARAMS; i++)
    fprintf(fid, "%s : %s\n", paramstr[i], in[i]);

  fclose(fid);

  printf("\n");

  return 0;

}

int generate_params_info(int num)
{
  switch(num)
  {
  case 0:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies which format your datafile is in, either\n");
      printf("SDPA or SDPLR format. SDPA is the most common format.\n\n");
      printf("Information on the SDPA format can be found at\n");
      printf("http://www.nmt.edu/~sdplib/FORMAT. Information on the SDPLR format\n");
      printf("can be found in the SDPLR v1.0 User's Guide.\n");
      printf("\n");
      break;
    }
  case 1:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies how accurately you would like to satisfy\n");
      printf("the primal constraints. SDPLR is a primal infeasible method that\n");
      printf("works towards feasibility as the algorithm progresses. Primal\n");
      printf("infeasibility is measured in relative terms as\n");
      printf("\n");
      printf("            || b - A(RR^T) || / (1 + |b_max|)\n");
      printf("\n");
      printf("where R is the decision variable, A represents the constraints, b\n");
      printf("is the right-hand side, and b_max is the largest entry of b in\n");
      printf("absolute value.\n");
      printf("\n");
      printf("Smaller values will cause SDPLR to do more work. On most problems,\n");
      printf("reasonable values are in the range 1.0e-5 to 1.0e-8.\n");
      printf("\n");
      break;
    }
  case 2:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies how accurately you would like to solve\n");
      printf("each augmented Lagrangian subproblem. SDPLR uses the augmented\n");
      printf("Lagrangian approach for nonlinear problems. Optimality for a\n");
      printf("subproblem is measured in relative terms as\n");
      printf("\n");
      printf("            || gradient || = || 2SR ||_F / (1 + |C_max|)\n");
      printf("\n");
      printf("where R is the decision variable, S is the dual variable estimate,\n");
      printf("C is the objective cost matrix, and C_max is the largest entry of\n");
      printf("C in absolute value. Based on the above formula, this parameter\n");
      printf("can be interpreted as enforcing complementary slackness between\n");
      printf("the primal and dual problems.\n");
      printf("\n");
      printf("Smaller values will cause SDPLR to do more work. On most problems,\n");
      printf("reasonable values are in the range 1.0e-1 and 1.0e-3.\n");
      printf("\n");
      break;
    }
  case 3:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies how the step direction will be calculated\n");
      printf("in each iteration. The LBFGS (limited memory BFGS) direction is\n");
      printf("quicker but worse quality. The TN (truncated Newton) direction is\n");
      printf("slower but better quality. On most problems, LBFGS performs better\n");
      printf("overall, but TN may be helpful on your problem.\n");
      printf("\n");
      printf("An additional option, Switch, is available to do LBFGS early on\n");
      printf("and then to switch to TN later.\n");
      printf("\n");
      printf("If LBFGS or Switch is chosen, then the parameter 'number of LBFGS\n");
      printf("vectors' must be specified.\n");
      printf("\n");
      printf("If TN or Switch is chosen, then the parameter 'auxiliary\n");
      printf("structures to facilitate TN' must be specified.\n");
      printf("\n");
      printf("If Switch is chosen, then the parameter 'feasibility level for\n");
      printf("switch' must be specified.\n");
      printf("\n");
      break;
    }
  case 4:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the factor by which the penalty parameter\n");
      printf("is increased every major iteration. SDPLR uses the penalty\n");
      printf("parameter to enforce feasibility in its augmented Lagrangian\n");
      printf("approach. This parameter should be greater than 1.0.\n");
      printf("\n");
      printf("Smaller values are considered more conservative; higher values are\n");
      printf("more aggressive. Reasonable values are between 2.0 and 10.0.\n");
      printf("\n");
      break;
    }
  case 5:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies whether or not to perform the rank\n");
      printf("reduction procedure, which is a dynamic way to reduce the\n");
      printf("dimensionality of the problem, thereby speeding up the algorithm.\n");
      printf("On many problems, rank reduction is very effective, but on some\n");
      printf("problems it can actually increase the overall time required.\n");
      printf("\n");
      break;
    }
  case 6:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the overall time limit for the algorithm\n");
      printf("(in seconds). The algorithm is terminated immediately after the\n");
      printf("completion of the first minor iteration past the time limit.\n");
      printf("\n");
      break;
    }
  case 7:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the amount of information printed to the\n");
      printf("screen by the algorithm: 0 prints nothing, and 1 prints\n");
      printf("everything.\n");
      printf("\n");
      break;
    }
  case 8:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the threshold N such that all matrices\n");
      printf("with both dimensions <= N are stored as dense matrices. The idea\n");
      printf("is that, even if a small matrix is sparse, it is often quicker to\n");
      printf("treat it as dense.\n");
      printf("\n");
      printf("In the current implementation of SDPLR, this parameter applies\n");
      printf("only to storage of the dual variable S, and not the data matrices\n");
      printf("C, A_i, which are always stored as sparse. The primal variable R\n");
      printf("is always stored as dense.\n");
      printf("\n");
      break;
    }
  case 9:
    {
      printf("\n");
      printf("\n");
      printf("This parameter specifies the threshold d such that all matrices\n");
      printf("with density >= d are stored as dense matrices. Here, d is in\n");
      printf("[0,1].\n");
      printf("\n");
      printf("In the current implementation of SDPLR, this parameter applies\n");
      printf("only to storage of the dual variable S, and not the data matrices\n");
      printf("C, A_i, which are always stored as sparse. The primal variable R\n");
      printf("is always stored as dense.\n");
      printf("\n");
      break;
    }
  case 10:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* Applicable only if parameter 'step direction calculation'\n");
      printf("is set to LBFSG=1 or to Switch=3.\n");
      printf("\n");
      printf("This parameter specifies how many limited memory BFGS (LBFGS)\n");
      printf("pairs to store. A higher number produces better LBFGS directions\n");
      printf("but requires more computation. Reasonable values are between 3 and\n");
      printf("10.\n");
      printf("\n");
      printf("In-development: This parameter also affects the algorithm when\n");
      printf("'preconditioning method' is set to 100. See the technical report\n");
      printf("\"Computational Enhancements in Low-Rank Semidefinite Programming\"\n");
      printf("available at http://dollar.biz.uiowa.edu/sburer or contact\n");
      printf("samuel-burer@uiowa.edu for more information.\n");
      printf("\n");
      break;
    }
  case 11:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* Applicable only if parameter 'step direction calculation'\n");
      printf("is set to Switch=3.\n");
      printf("\n");
      printf("This parameter specifies the feasibility level at which to switch\n");
      printf("from LBFGS to TN. More precisely, the switch occurs after the\n");
      printf("first major iteration for which the primal infeasibility passes\n");
      printf("below the parameter. Primal infeasibility is measured in relative\n");
      printf("terms as\n");
      printf("\n");
      printf("            || b - A(RR^T) || / (1 + |b_max|)\n");
      printf("\n");
      printf("where R is the decision variable, A represents the constraints, b\n");
      printf("is the right-hand side, and b_max is the largest entry of b in\n");
      printf("absolute value.\n");
      printf("\n");
      printf("This parameter should be set relative to your choice for the\n");
      printf("'feasibility tolerance' parameter.\n");
      printf("\n");
      break;
    }
  case 12:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* Applicable only if parameter 'step direction calculation'\n");
      printf("is set to TN=2 or to Switch=3.\n");
      printf("\n");
      printf("This parameter specifies whether or not to utilize extra storage\n");
      printf("space to aid the TN method. The main benefit is increased\n");
      printf("numerical accuracy. Usually there is no noticeable change in\n");
      printf("speed. On a few problems, the extra storage requirements may be\n");
      printf("prohibitive.\n");
      printf("\n");
      break;
    }
  case 13:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* Only applicable if parameter 'rank reduction' is set to\n");
      printf("yes=1.\n");
      printf("\n");
      printf("This parameter controls the determination of numerical rank of the\n");
      printf("primal variable R in the rank reduction procedure. Let the\n");
      printf("parameter be called TOL. In essence, if a certain submatrix M (k\n");
      printf("columns) of R (r columns) satisfies\n");
      printf("\n");
      printf("           || M ||_F <= TOL * || R ||_F,\n");
      printf("\n");
      printf("then the numerical rank of R is r-k. A good choice for TOL is the\n");
      printf("machine precision. If TOL is set too high, then the rank reduction\n");
      printf("procedure may become unreliable.\n");
      printf("\n");
      break;
    }
  case 14:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development! See the technical report \"Computational Enhancements\n");
      printf("in Low-Rank Semidefinite Programming\" available at\n");
      printf("http://dollar.biz.uiowa.edu/sburer or contact\n");
      printf("samuel-burer@uiowa.edu for more information.\n");
      printf("\n");
      printf("*NOTE* Applicable only if parameter 'step direction calculation'\n");
      printf("is set to TN=2 or to Switch=3.\n");
      printf("\n");
      printf("This parameter specifies which preconditioning method to use for\n");
      printf("the TN method.\n");
      printf("\n");
      break;
    }
  case 15:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development! See the technical report \"Computational Enhancements\n");
      printf("in Low-Rank Semidefinite Programming\" available at\n");
      printf("http://dollar.biz.uiowa.edu/sburer or contact\n");
      printf("samuel-burer@uiowa.edu for more information.\n");
      printf("\n");
      printf("*NOTE* Applicable only if parameter 'step direction calculation'\n");
      printf("is set to TN=2 or to Switch=3 and if parameter 'preconditioning\n");
      printf("method' is set to 11, 12, 13, 21, 22, or 23.\n");
      printf("\n");
      printf("This parameter specifies the percentage of fill-in with which the\n");
      printf("incomplete Cholesky factorization of S is calculated. Valid values\n");
      printf("are in [0,1]. A value of 0.0 means \"no fill-in,\" and a value of\n");
      printf("1.0 means \"all fill-in induced by S.\"\n");
      printf("\n");
      break;
    }
  case 16:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development! See the technical report \"Computational Enhancements\n");
      printf("in Low-Rank Semidefinite Programming\" available at\n");
      printf("http://dollar.biz.uiowa.edu/sburer or contact\n");
      printf("samuel-burer@uiowa.edu for more information.\n");
      printf("\n");
      printf("*NOTE* Applicable only if parameter 'step direction calculation'\n");
      printf("is set to TN=2 or to Switch=3 and if parameter 'preconditioning\n");
      printf("method' is set to 11, 12, 13, 21, 22, or 23.\n");
      printf("\n");
      printf("*NOTE* Applicable only if the library METIS has been linked in\n");
      printf("with SDPLR.\n");
      printf("\n");
      printf("This parameter specifies whether or not to reorder the data so as\n");
      printf("to achieve less fill-in in the incomplete Cholesky factorization\n");
      printf("of S.\n");
      printf("\n");
      break;
    }
  case 17:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development but currently inactive! It has to do with calculating\n");
      printf("dual bounds with SDPLR. Contact samuel-burer@uiowa.edu for more\n");
      printf("information.\n");
      printf("\n");
      printf("This parameter specifies the target duality gap. SDPLR terminates\n");
      printf("once both the target feasibility and gap are met.\n");
      printf("\n");
      break;
    }
  case 18:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development but currently inactive! It has to do with calculating\n");
      printf("dual bounds with SDPLR. Contact samuel-burer@uiowa.edu for more\n");
      printf("information.\n");
      printf("\n");
      printf("This parameter specifies when to calculate a dual bound. The value\n");
      printf("-1 means never, 0 means only at the end of the algorithm, 1 means\n");
      printf("after every major iteration.\n");
      printf("\n");
      break;
    }
  case 19:
    {
      printf("\n");
      printf("\n");
      printf("*NOTE* This parameter refers to a feature of SDPLR that is in\n");
      printf("development but currently inactive! It has to do with calculating\n");
      printf("dual bounds with SDPLR. Contact samuel-burer@uiowa.edu for more\n");
      printf("information.\n");
      printf("\n");
      printf("This parameter specifies what type of dual bound to calculate. The\n");
      printf("value 1 refers to +/-1 combinatorial bounds, and the value 2\n");
      printf("refers to 0/1 combinatorial bounds.\n");
      printf("\n");
      break;
    }
  default:
    {
      printf("default\n");
      break;
    }
  }

  fflush(stdout);

  return 0;

}



