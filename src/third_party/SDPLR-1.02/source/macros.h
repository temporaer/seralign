#ifndef __PC
#define _isnan isnan
#endif

#ifdef __MEX
#define printf mexPrintf
#endif

#define MYCALLOC(VAR,TYPE,SIZE) VAR = (TYPE*)calloc(SIZE, sizeof(TYPE))
#define MYFREE(VAR) free(VAR);

#define mymax(A, B) ((A) > (B) ? (A) : (B))
#define mymin(A, B) ((A) < (B) ? (A) : (B))

#define SMATIND(I,J,DIM) (J-1)*DIM + I
#define LMATIND(I,J,DIM) DIM*(DIM+1)/2 -(DIM-J+1)*(DIM-J+2)/2 + I - J + 1

#define EASYDAXPY(DIM,SCAL,X,Y)   mydaxpy(DIM,SCAL,X+1,1,Y+1,1)
#define EASYDCOPY(DIM,X,Y)        mydcopy(DIM,X+1,1,Y+1,1)
#define EASYDDOT(DIM,X,Y)         myddot(DIM,X+1,1,Y+1,1)
#define EASYDNRM2(DIM,X)          mydnrm2(DIM,X+1,1)
#define EASYDSCAL(DIM,SCAL,X)     mydscal(DIM,SCAL,X+1,1)
#define ZEROVEC(X,DIM)            mydscal(DIM,0.0,X+1,1)

/*
#define __OPT_MYBLAS        0     // Not including dot product or norm
#define __OPT_ALIAS         0
#define __OPT_BACKLOOP      0
#define __OPT_UNROLLEDLOOP  0
#define __OPT_CONST         0
*/

// Some parameters for testing purposes
#define RANDOM 1 // default 1
#define SCALE_OBJ 1

#define SDPBLK  's'
#define DIAGBLK 'd'
#define SPARSE  's'
#define DENSE   'd'

#define PC_S_DIAG         0  // S approximated by diagonal
#define PC_S_CHOL        10  // S approximated by incomplete Cholesky
#define PC_S_FACT        20  // S approximated by incomplete Cholesky, then LL + WW^T = L(I + UU^T)L^T is approximated
#define PC_DIAG           1  // WW^T is approximated by diagonal
#define PC_LARGEVECS      2  // WW^T is approxiamted by large members of AiR
#define PC_EIGENVECS      3  // WW^T is approximated by top eigenvectors of AiR
#define PC_LBFGS        100  // L-BFGS preconditioner of Morales-Nocedal
#define PC_LBFGS_NEW    200  // My new attempt (let's hope it works!)
