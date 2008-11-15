#include "user.h"

int user_readdata(userdata* user, FILE* datafile, char* label)
{
  int i, j, h, symma=1, symmb=1, cta=0, ctb=0;

  user->QAP_n = 0;

  if(strcmp(label, "QAP1") == 0) {

    fscanf(datafile, "%d\n", &(user->QAP_n));

    user_MYCALLOC(user->QAP_A, double*, user->QAP_n + 1);
    for(i = 1; i <= user->QAP_n; i++)
      user_MYCALLOC(user->QAP_A[i], double, user->QAP_n + 1);

    user_MYCALLOC(user->QAP_B, double*, user->QAP_n + 1);
    for(i = 1; i <= user->QAP_n; i++)
      user_MYCALLOC(user->QAP_B[i], double, user->QAP_n + 1);

    user_MYCALLOC(user->QAP_temp1, double*, user->QAP_n + 1);
    for(i = 1; i <= user->QAP_n; i++)
      user_MYCALLOC(user->QAP_temp1[i], double, user->QAP_n + 1);

    user_MYCALLOC(user->QAP_temp2, double*, user->QAP_n + 1);
    for(i = 1; i <= user->QAP_n; i++)
      user_MYCALLOC(user->QAP_temp2[i], double, user->QAP_n + 1);

    for(i = 1; i <= user->QAP_n; i++) {
      for(j = 1; j <= user->QAP_n; j++)
        fscanf(datafile, "%lf ", &(user->QAP_A[i][j]));
      fscanf(datafile, "\n");
    }

    for(i = 1; i <= user->QAP_n; i++) {
      for(j = 1; j <= user->QAP_n; j++)
        fscanf(datafile, "%lf ", &(user->QAP_B[i][j]));
      fscanf(datafile, "\n");
    }

    for(i = 1; i <= user->QAP_n; i++) {
      for(j = 1; j <= i-1; j++) {
        if(user->QAP_A[i][j] != user->QAP_A[j][i]) symma = 0;
        if(user->QAP_B[i][j] != user->QAP_B[j][i]) symmb = 0;
      }
    }

    if(symma == 0 && symmb == 0) {
      printf("Cannot easily symmetrize...exiting\n");
      exit(0);
    }

    for(i = 1; i <= user->QAP_n; i++) {
      for(j = 1; j <= i-1; j++) {
        user->QAP_A[i][j] = user->QAP_A[j][i] = 0.5*(user->QAP_A[i][j] + user->QAP_A[j][i]);
        user->QAP_B[i][j] = user->QAP_B[j][i] = 0.5*(user->QAP_B[i][j] + user->QAP_B[j][i]);
      }
    }

    for(i = 1; i <= user->QAP_n; i++) for(j = 1; j <= i; j++) {
      if(user->QAP_A[i][j] != 0.0) cta++;
      if(user->QAP_B[i][j] != 0.0) ctb++;
    }

    user_MYCALLOC(user->QAP_Ai, int, cta + 1);
    user_MYCALLOC(user->QAP_Aj, int, cta + 1);
    user_MYCALLOC(user->QAP_Aent, double, cta + 1);
    user_MYCALLOC(user->QAP_Bi, int, ctb + 1);
    user_MYCALLOC(user->QAP_Bj, int, ctb + 1);
    user_MYCALLOC(user->QAP_Bent, double, ctb + 1);

    h=1;
    for(i = 1; i <= user->QAP_n; i++) for(j = 1; j <= i; j++) {
      if(user->QAP_A[i][j] != 0.0) {
        user->QAP_Ai[h] = i;
        user->QAP_Aj[h] = j;
        user->QAP_Aent[h++] = user->QAP_A[i][j];
      }
    }
    if(h != cta + 1) { printf("problem!!!!\n"); }
    user->QAP_ctA = cta;

    h=1;
    for(i = 1; i <= user->QAP_n; i++) for(j = 1; j <= i; j++) {
      if(user->QAP_B[i][j] != 0.0) {
        user->QAP_Bi[h] = i;
        user->QAP_Bj[h] = j;
        user->QAP_Bent[h++] = user->QAP_B[i][j];
      }
    }
    if(h != ctb + 1) { printf("problem!!!!\n"); }
    user->QAP_ctB = ctb;

  }

  return 1;
}

int user_deinitialize(userdata* user)
{
  int i;

  for(i = 1; i <= user->QAP_n; i++)
    if(user->QAP_A[i] != NULL) user_MYFREE(user->QAP_A[i]);
  if(user->QAP_A != NULL) user_MYFREE(user->QAP_A);

  for(i = 1; i <= user->QAP_n; i++)
    if(user->QAP_B[i] != NULL) user_MYFREE(user->QAP_B[i]);
  if(user->QAP_B != NULL) user_MYFREE(user->QAP_B);

  for(i = 1; i <= user->QAP_n; i++)
    if(user->QAP_temp1[i] != NULL) user_MYFREE(user->QAP_temp1[i]);
  if(user->QAP_temp1 != NULL) user_MYFREE(user->QAP_temp1);

  for(i = 1; i <= user->QAP_n; i++)
    if(user->QAP_temp2[i] != NULL) user_MYFREE(user->QAP_temp2[i]);
  if(user->QAP_temp2 != NULL) user_MYFREE(user->QAP_temp2);

  if(user->QAP_SR != NULL) user_MYFREE(user->QAP_SR);

  if(user->QAP_Ai != NULL) user_MYFREE(user->QAP_Ai);
  if(user->QAP_Aj != NULL) user_MYFREE(user->QAP_Aj);
  if(user->QAP_Aent != NULL) user_MYFREE(user->QAP_Aent);
  if(user->QAP_Bi != NULL) user_MYFREE(user->QAP_Bi);
  if(user->QAP_Bj != NULL) user_MYFREE(user->QAP_Bj);
  if(user->QAP_Bent != NULL) user_MYFREE(user->QAP_Bent);

  return 1;
}


double user_AdotRRt(userdata* user, char* label, user_lowrankmat* R)
{
  int j;
  double sum=0.0;
  /* int i, j, k, n, r, tempintai, tempintaj, tempintbi, tempintbj;
  double sum=0.0, tempval, tempvala, tempvalb;

  n = user->QAP_n;
  r = R->ncol;

  if(strcmp(label, "QAP1") == 0) {
    for(i = 2; i <= n*n+1; i++)
      for(j = 2; j <= i; j++) {

        tempval = 0.0;
        for(k = 1; k <= r; k++)
          tempval += R->col[k][i]*R->col[k][j];

        if((i-1)%n == 0) { tempintai = (i-1)/n;     tempintbi = n;       }
        else             { tempintai = (i-1)/n + 1; tempintbi = (i-1)%n; }
        if((j-1)%n == 0) { tempintaj = (j-1)/n;     tempintbj = n;       }
        else             { tempintaj = (j-1)/n + 1; tempintbj = (j-1)%n; }

        tempvala = user->QAP_A[tempintai][tempintaj];
        tempvalb = user->QAP_B[tempintbi][tempintbj];

        if(i == j) sum +=     tempvala*tempvalb*tempval;
        else       sum += 2.0*tempvala*tempvalb*tempval;

      }
  } */

  if(strcmp(label, "QAP1") == 0) {
    for(j = 1; j <= R->ncol; j++) {
      user_Stimesvec(user, label, R->col[j], user->QAP_SR->col[j], 0);
      sum += EASYDDOT(user->QAP_n*user->QAP_n + 1, R->col[j], user->QAP_SR->col[j]);
    }
  }

  return sum;
}

int user_Stimesvec(userdata* user, char* label, double* in, double* result, int flag)
{
  int i, j, k, n;
  double *vecin, *vecout;

  if(strcmp(label, "QAP1") == 0) {

    if(flag == 0) {

      n = user->QAP_n;
      vecin = in + 1;
      result[1] = 0;
      vecout = result + 1;

      for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
        user->QAP_temp1[i][j] = vecin[n*(j-1) + i];

      for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
        user->QAP_temp2[i][j] = 0.0;

      //for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
      //  for(k = 1; k <= n; k++)
      //    user->QAP_temp2[i][j] += user->QAP_temp1[i][k]*user->QAP_A[k][j];

      for(k = 1; k <= user->QAP_ctB; k++) {
        for(i = 1; i <= n; i++)
          user->QAP_temp2[i][user->QAP_Bj[k]] += user->QAP_temp1[i][user->QAP_Bi[k]]*user->QAP_Bent[k];
        if(user->QAP_Bi[k] != user->QAP_Bj[k])
          for(i = 1; i <= n; i++)
            user->QAP_temp2[i][user->QAP_Bi[k]] += user->QAP_temp1[i][user->QAP_Bj[k]]*user->QAP_Bent[k];
      }

      for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
        user->QAP_temp1[i][j] = 0.0;

      //for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
      //  for(k = 1; k <= n; k++)
      //    user->QAP_temp1[i][j] += user->QAP_A[i][k]*user->QAP_temp2[k][j];

      for(k = 1; k <= user->QAP_ctA; k++) {
        for(j = 1; j <= n; j++)
          user->QAP_temp1[user->QAP_Ai[k]][j] += user->QAP_Aent[k]*user->QAP_temp2[user->QAP_Aj[k]][j];
        if(user->QAP_Ai[k] != user->QAP_Aj[k])
          for(j = 1; j <= n; j++)
            user->QAP_temp1[user->QAP_Aj[k]][j] += user->QAP_Aent[k]*user->QAP_temp2[user->QAP_Ai[k]][j];

      }

      for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
        vecout[n*(j-1) + i] = user->QAP_temp1[i][j];
    }

    else {
      n = user->QAP_n;
      mydcopy(n*n + 1, user->QAP_SR->col[flag] + 1, 1, result + 1, 1);
    }

  }

  return 1;
}

double user_normdatamat(userdata* user, char* label)
{
  int i, j, n;
  double sum=0.0, tempval=0.0, tempval2=0.0;

  if(strcmp(label, "QAP1") == 0) {
    n = user->QAP_n;
    for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
      tempval += pow(user->QAP_B[i][j], 2.0);
    for(i = 1; i <= n; i++) for(j = 1; j <= n; j++)
      tempval2 += pow(user->QAP_A[i][j], 2.0);
    sum = sqrt(tempval2*tempval);
  }

  return sum;
}

