#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define user_MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE));
#define user_MYFREE(VAR) free(VAR);

typedef struct {
  double** col;
  int      ncol;
} user_lowrankmat;

typedef struct {
  int      QAP_n;
  double **QAP_A;
  int     *QAP_Ai;
  int     *QAP_Aj;
  int      QAP_ctA;
  double  *QAP_Aent;
  double **QAP_B;
  int     *QAP_Bi;
  int     *QAP_Bj;
  int      QAP_ctB;
  double  *QAP_Bent;
  double **QAP_temp1;
  double **QAP_temp2;
  int      QAP_flag;
  user_lowrankmat *QAP_SR;
} userdata;


int user_readdata(userdata*, FILE*, char*);
int user_deinitialize(userdata*);
double user_AdotRRt(userdata*, char*, user_lowrankmat*);
int user_Stimesvec(userdata*, char*, double*, double*, int);
double user_normdatamat(userdata*, char*);
int mydcopy(int, double*, int, double*, int);
double myddot(int, double*, int, double*, int);
#define EASYDDOT(DIM,X,Y)         myddot(DIM,X+1,1,Y+1,1)

